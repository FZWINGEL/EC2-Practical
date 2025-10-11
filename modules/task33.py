import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from .config import (
    DATA_DIR,
    GC_AREA_CM2,
    NU_CMS2,
    N_ELECTRONS,
    R_GAS,
    F_CONST,
    T_K,
    C_BULK_MOL_PER_CM3,
)
from .diffusion import get_calculated_diffusion, add_diffusion_coefficient
from .io_utils import read_biologic_table, write_csv
from .pipeline import TaskReport
from .plotting import beautify_axes, safe_save
from .task32_eis import fit_with_impedance_py
from .task32_eis import ECFFit
from .utils import extract_rpm_from_name, get_column, label_from_filename


@dataclass
class KLPoint:
    E_V: float
    inv_sqrt_omega_s05: float
    inv_j_cm2_A_inv: float
    jk_A_cm2: float


def _interpolate_to_common_grid(curves: Dict[int, pd.DataFrame]) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    Emins, Emaxs = [], []
    for df in curves.values():
        e = pd.to_numeric(df["Ewe/V"], errors="coerce").to_numpy()
        Emins.append(np.nanmin(e))
        Emaxs.append(np.nanmax(e))
    Emin, Emax = max(Emins), min(Emaxs)
    if not np.isfinite(Emin) or not np.isfinite(Emax) or Emin >= Emax:
        raise ValueError("LSV potentials have no overlap across rpms.")
    E_grid = np.linspace(Emin, Emax, 200)

    j_by_rpm: Dict[int, np.ndarray] = {}
    for rpm, df in curves.items():
        e = pd.to_numeric(df["Ewe/V"], errors="coerce").to_numpy()
        current_candidates = [col for col in ["<I>/mA", "I/mA"] if col in df.columns]
        if not current_candidates:
            raise KeyError("No current column '<I>/mA' or 'I/mA' found in LSV data.")
        i_mA = pd.to_numeric(df[current_candidates[0]], errors="coerce").to_numpy()
        j_A_cm2 = (i_mA / 1000.0) / GC_AREA_CM2
        order = np.argsort(e)
        e_sorted = e[order]
        j_sorted = j_A_cm2[order]
        j_interp = np.interp(E_grid, e_sorted, j_sorted)
        j_by_rpm[rpm] = j_interp
    return E_grid, j_by_rpm


def analyze_task33_lsv_and_eis() -> TaskReport:
    report = TaskReport(name="Task 3.3 - LSV/EIS")
    lsv_dir = DATA_DIR / "Task 3.3 LSV"
    if not lsv_dir.exists():
        message = "[Task 3.3] Directory not found; skipping."
        print(message)
        report.add_message(message)
        return report

    lsv_files = sorted([p for p in lsv_dir.glob("*_LSV_C01.txt") if "PEIS" not in p.name])
    curves: Dict[int, pd.DataFrame] = {}

    if lsv_files:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        plot_data = []  # Store plot data for sorting
        
        for path in lsv_files:
            try:
                df = read_biologic_table(path)
            except Exception as exc:
                message = f"[LSV] Failed to read {path.name}: {exc}"
                print(message)
                report.add_warning(message)
                continue

            e_col = get_column(df, ["Ewe/V"]) or "Ewe/V"
            i_col = get_column(df, ["<I>/mA", "I/mA"]) or "<I>/mA"
            if not ({e_col, i_col} <= set(df.columns)):
                message = f"[LSV] Missing columns in {path.name}; skipping."
                print(message)
                report.add_warning(message)
                continue

            e_v = pd.to_numeric(df[e_col], errors="coerce").to_numpy()
            j_mA_cm2 = pd.to_numeric(df[i_col], errors="coerce").to_numpy() / GC_AREA_CM2
            rpm = extract_rpm_from_name(path)
            label = f"{rpm} rpm" if rpm is not None else label_from_filename(path)
            
            # Store plot data with rpm for sorting
            plot_data.append((rpm if rpm is not None else 0, e_v, j_mA_cm2, label))

            if rpm is not None:
                df_use = df[[e_col, i_col]].copy()
                df_use.columns = ["Ewe/V", i_col]
                if i_col != "<I>/mA":
                    df_use["<I>/mA"] = df_use[i_col]
                curves[rpm] = df_use

        # Sort by rpm and plot
        plot_data.sort(key=lambda x: x[0])
        for _, e_v, j_mA_cm2, label in plot_data:
            ax.plot(e_v, j_mA_cm2, lw=1.0, label=label)

        ax.set_xlabel("E [V vs Ag/AgCl]")
        ax.set_ylabel("j [mA cm$^{-2}$]")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        report.record_figure(safe_save(fig, "T3.3_LSV_overlay.png"))
    else:
        message = "[LSV] No LSV files found; skipping LSV plot."
        print(message)
        report.add_message(message)

    if len(curves) >= 3:
        try:
            E_grid, j_by_rpm = _interpolate_to_common_grid(curves)
        except Exception as exc:
            message = f"[KL] Could not interpolate to common grid: {exc}"
            print(message)
            report.add_warning(message)
            E_grid, j_by_rpm = None, None

        if E_grid is not None and j_by_rpm is not None:
            rpms = np.array(sorted(j_by_rpm.keys()))
            omega = 2 * np.pi * (rpms / 60.0)
            inv_sqrt_omega = 1.0 / np.sqrt(omega)

            E_targets = np.linspace(E_grid.min(), E_grid.max(), 50)  # Use all potentials
            kl_rows: List[KLPoint] = []
            jk_vs_E: List[Tuple[float, float]] = []
            kl_summ_rows: List[Dict[str, float]] = []
            accepted_E: List[float] = []

            for E0 in E_targets:
                j_at_E = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E0))] for rpm in rpms])
                mask = np.isfinite(j_at_E) & (np.abs(j_at_E) > 1e-6)
                if mask.sum() < 3:
                    continue
                y = 1.0 / np.abs(j_at_E[mask])
                x = inv_sqrt_omega[mask]

                # Do KL analysis on all potentials
                if len(x) >= 3:
                    A = np.vstack([x, np.ones_like(x)]).T
                    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                    y_pred = slope * x + intercept
                    ss_res = float(np.sum((y - y_pred) ** 2))
                    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
                    r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                    if (r2 >= 0.90) and np.isfinite(slope) and (slope > 0) and np.isfinite(intercept) and (intercept != 0):
                        jk = 1.0 / intercept
                        kl_rows.extend(
                            KLPoint(E_V=float(E0), inv_sqrt_omega_s05=float(xi), inv_j_cm2_A_inv=float(yi), jk_A_cm2=float(jk))
                            for xi, yi in zip(x, y)
                        )
                        jk_vs_E.append((E0, jk))
                        accepted_E.append(float(E0))

                # KL slope S = 1/B; B = 0.62 n F D^(2/3) nu^(-1/6) C*
                # Work on current density, so area cancels
                if 'slope' in locals() and np.isfinite(slope) and (slope > 0):
                    B_val = 1.0 / float(slope)
                    denom = 0.62 * N_ELECTRONS * F_CONST * C_BULK_MOL_PER_CM3
                    if denom > 0 and (NU_CMS2 > 0):
                        D_pow23 = (B_val * (NU_CMS2 ** (1.0 / 6.0))) / denom
                        D_val = float(D_pow23 ** 1.5) if D_pow23 > 0 else float("nan")
                    else:
                        D_val = float("nan")
                    kl_summ_rows.append({
                        "E_V": float(E0),
                        "slope_1_over_B": float(slope),
                        "B_A_cm2_s05": float(B_val),
                        "D_KL_cm2_s": float(D_val),
                        "R2": float(r2),
                    })

            if jk_vs_E:
                df_kl = pd.DataFrame([point.__dict__ for point in kl_rows])
                report.record_table(write_csv(df_kl, "T3.3_KouteckyLevich_points.csv"))
                report.add_message(f"[KL] Logged {len(kl_rows)} Koutecky-Levich points across {len(jk_vs_E)} potentials.")

                df_jk = pd.DataFrame({"E_V": [p[0] for p in jk_vs_E], "j_k_A_cm2": [p[1] for p in jk_vs_E]}).sort_values("E_V")
                report.record_table(write_csv(df_jk, "T3.3_jk_vs_E.csv"))
                report.add_message(f"[KL] Derived {len(df_jk)} kinetic-current values.")

                # Persist KL-derived diffusion coefficients per potential and summarize
                if kl_summ_rows:
                    df_klD = pd.DataFrame(kl_summ_rows).sort_values("E_V")
                    report.record_table(write_csv(df_klD, "T3.3_KL_diffusion.csv"))
                    D_vals = df_klD["D_KL_cm2_s"].to_numpy()
                    D_vals = D_vals[np.isfinite(D_vals) & (D_vals > 0)]
                    if D_vals.size:
                        D_med = float(np.nanmedian(D_vals))
                        add_diffusion_coefficient(D_med, "Koutecky-Levich (Task 3.3)")
                        report.add_message(f"[KL] D_KL~{D_med:.2e} cm^2/s (median across E)")

                mid = len(accepted_E) // 2
                if mid < len(accepted_E) and len(accepted_E) > 0:
                    E_mid = accepted_E[mid]
                    j_mid = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E_mid))] for rpm in rpms])
                    mask = np.isfinite(j_mid) & (np.abs(j_mid) > 1e-6)
                    x = inv_sqrt_omega[mask]
                    y = 1.0 / np.abs(j_mid[mask])
                    A = np.vstack([x, np.ones_like(x)]).T
                    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                    figkl, axkl = plt.subplots(figsize=(5.4, 4.0))
                    axkl.plot(x, y, "o")
                    axkl.plot(x, slope * x + intercept, "-")
                    axkl.set_xlabel(r"$1/\sqrt{\omega}$ [s$^{1/2}$]")
                    axkl.set_ylabel(r"$1/j$ [cm$^2$ A$^{-1}$]")
                    axkl.set_title(f"Koutecky-Levich at E = {E_mid:.3f} V")
                    beautify_axes(axkl)
                    report.record_figure(safe_save(figkl, "T3.3_KL_example.png"))

                zero_Es = []
                for rpm, j in j_by_rpm.items():
                    sign = np.sign(j)
                    zcross = np.where(np.diff(sign))[0]
                    if len(zcross):
                        idx = zcross[0]
                        E0 = float(np.interp(0.0, [j[idx], j[idx + 1]], [E_grid[idx], E_grid[idx + 1]]))
                        zero_Es.append(E0)
                if zero_Es:
                    E_eq = float(np.median(zero_Es))
                else:
                    idx = int(np.argmin(np.abs(df_jk["j_k_A_cm2"])) )
                    E_eq = float(df_jk["E_V"].iloc[idx])

                jk_data = df_jk.copy()
                jk_data["eta_V"] = jk_data["E_V"] - E_eq
                jk_data = jk_data[np.isfinite(jk_data["eta_V"]) & np.isfinite(jk_data["j_k_A_cm2"])].copy()
                jk_data["sign"] = np.sign(jk_data["j_k_A_cm2"])
                sign_to_use = 1.0 if (jk_data["sign"] > 0).sum() >= (jk_data["sign"] < 0).sum() else -1.0
                branch = jk_data[(jk_data["sign"] == sign_to_use) & (np.abs(jk_data["j_k_A_cm2"]) > 1e-8)].copy()

                if len(branch) >= 3:
                    branch["log10_jk"] = np.log10(np.abs(branch["j_k_A_cm2"]))
                    
                    # Select points in the log10(|jk|) range -2.5 to -1.35
                    selected = (branch["log10_jk"] >= -2.5) & (branch["log10_jk"] <= -1.35)
                    
                    if selected.sum() >= 3:
                        x = branch.loc[selected, "log10_jk"].to_numpy()
                        y = branch.loc[selected, "eta_V"].to_numpy()
                        A = np.vstack([x, np.ones_like(x)]).T
                        slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                        y_pred = slope * x + intercept
                        ss_res = np.sum((y - y_pred) ** 2)
                        ss_tot = np.sum((y - np.mean(y)) ** 2)
                        r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                        best = {
                            "log10_jk_min": -2.5,
                            "log10_jk_max": -1.35,
                            "slope": float(slope),
                            "intercept": float(intercept),
                            "r2": float(r2),
                            "x": x,
                            "y": y,
                            "n_points": len(x),
                        }
                    else:
                        # Fallback: use all points if not enough in target range
                        x = branch["log10_jk"].to_numpy()
                        y = branch["eta_V"].to_numpy()
                        A = np.vstack([x, np.ones_like(x)]).T
                        slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                        y_pred = slope * x + intercept
                        ss_res = np.sum((y - y_pred) ** 2)
                        ss_tot = np.sum((y - np.mean(y)) ** 2)
                        r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                        best = {
                            "log10_jk_min": float(branch["log10_jk"].min()),
                            "log10_jk_max": float(branch["log10_jk"].max()),
                            "slope": float(slope),
                            "intercept": float(intercept),
                            "r2": float(r2),
                            "x": x,
                            "y": y,
                            "n_points": len(x),
                        }

                    if best is not None:
                        log10_j0 = -best["intercept"] / best["slope"]
                        j0 = 10 ** log10_j0

                        figtf, axtf = plt.subplots(figsize=(5.4, 4.0))
                        
                        # Plot all points in gray
                        axtf.plot(branch["log10_jk"], branch["eta_V"], "o", color="lightgray", alpha=0.6, label=f"All Data ({len(branch)} points)")
                        
                        # Highlight selected points used for fitting
                        axtf.plot(best["x"], best["y"], "o", color="blue", label=f"Selected for Fit ({best['n_points']} points)")
                        
                        # Extend fit line across the whole plot, at least to j0 intercept
                        xlim = axtf.get_xlim()
                        x_intersect = -best["intercept"] / best["slope"]  # When η=0
                        x_min = min(xlim[0], x_intersect - 0.1)  # Extend left to include j0
                        xfit = np.linspace(x_min, xlim[1], 100)
                        axtf.plot(xfit, best["slope"] * xfit + best["intercept"], "-", color="red", label=f"Fit (R²={best['r2']:.3f})")
                        
                        # Calculate intersection with y-axis (η=0)
                        x_intersect = -best["intercept"] / best["slope"]  # When η=0
                        axtf.plot(x_intersect, 0, "ro", markersize=8, label=f"j₀: {j0:.2e} A/cm²")
                        
                        axtf.set_xlabel(r"log$_{10}$(|j$_k$|) [A cm$^{-2}$]")
                        axtf.set_ylabel(r"$\eta$ [V]")
                        beautify_axes(axtf)
                        axtf.legend(loc="upper right", fontsize=8)
                        report.record_figure(safe_save(figtf, "T3.3_Tafel.png"))

                        # Also estimate k0 from Tafel j0 if PEIS is not available later (fallback)
                        k0_app = float(j0 / (N_ELECTRONS * F_CONST * C_BULK_MOL_PER_CM3)) if (N_ELECTRONS > 0 and F_CONST > 0 and C_BULK_MOL_PER_CM3 > 0) else float("nan")
                        report.record_table(
                            write_csv(
                                pd.DataFrame(
                                    {
                                        "slope_V_per_decade": [best["slope"]],
                                        "intercept_V": [best["intercept"]],
                                        "E_eq_est_V": [E_eq],
                                        "j0_A_cm2": [float(j0)],
                                        "k0_app_cm_s": [float(k0_app)],
                                        "R2": [best["r2"]],
                                        "branch_sign": [int(sign_to_use)],
                                        "log10_jk_min": [best["log10_jk_min"]],
                                        "log10_jk_max": [best["log10_jk_max"]],
                                        "n_points_used": [best["n_points"]],
                                        "n_points_total": [len(branch)],
                                    }
                                ),
                                "T3.3_Tafel_fit.csv",
                            )
                        )
                        report.add_message(f"[Tafel] j0~{j0:.2e} A/cm^2; slope={best['slope']:.3f} V/dec")
                    else:
                        message = "[Tafel] Not enough points in any reasonable region to fit."
                        print(message)
                        report.add_warning(message)
                else:
                    message = "[Tafel] Not enough points on a single branch to fit."
                    print(message)
                    report.add_warning(message)
    else:
        message = "[KL] Need LSV data at >=3 rotation rates for Koutecky-Levich and Tafel."
        print(message)
        report.add_warning(message)

    peis_files = sorted([p for p in lsv_dir.glob("*_PEIS_C01.txt")])
    if peis_files:
        D_calc = get_calculated_diffusion()
        if not np.isfinite(D_calc):
            message = "[PEIS] No calculated D available yet (RS/Sand/Warburg). delta will be NaN. Run Task 3.2 first in the same run."
            print(message)
            report.add_warning(message)
        fig, ax = plt.subplots(figsize=(6.0, 4.8))
        peis_rows = []
        for path in peis_files:
            try:
                df = read_biologic_table(path)
            except Exception as exc:
                message = f"[PEIS] Failed to read {path.name}: {exc}"
                print(message)
                report.add_warning(message)
                continue

            re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
            im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
            if not ({re_col, im_col} <= set(df.columns)):
                message = f"[PEIS] Missing impedance columns in {path.name}; skipping."
                print(message)
                report.add_warning(message)
                continue

            z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
            z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
            rpm = extract_rpm_from_name(path)
            label = f"{rpm} rpm" if rpm is not None else label_from_filename(path)
            ax.plot(z_re, z_im, marker="o", ms=2.5, lw=0.8, label=label)

            f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
            f_all = pd.to_numeric(df[f_col], errors="coerce").to_numpy()
            mask = np.isfinite(f_all) & (f_all > 0)
            omega_all = 2 * np.pi * f_all[mask]

            try:
                fit = fit_with_impedance_py(df, circuit_str='R0-p(R1,CPE1)-p(R2,CPE2)')
            except Exception as exc:
                message = f"[PEIS] Fit failed for {path.name}: {exc}"
                print(message)
                report.add_warning(message)
                fit = ECFFit(model="", Rs=float("nan"), Rct=float("nan"), Cdl=None, sigma=None, Q=None, alpha=None, rchi2=float("nan"), n_points=0)

            if fit.Rs is not None and hasattr(fit, 'parameters_') and len(fit.parameters_) >= 7:
                Rs_area = fit.parameters_[0] * GC_AREA_CM2
                R1 = float(fit.parameters_[1])
                Q1 = float(fit.parameters_[2])
                alpha1 = float(fit.parameters_[3])
                R2 = float(fit.parameters_[4])
                # Q2 = float(fit.parameters_[5])  # not used currently
                # alpha2 = float(fit.parameters_[6])
                Rct_area = (R1 + R2) * GC_AREA_CM2  # R1 + R2
                # Effective Cdl from first branch at its peak frequency: C_eff = Q^(1/alpha) / R^(1 - 1/alpha)
                if np.isfinite(R1) and R1 > 0 and np.isfinite(Q1) and Q1 > 0 and np.isfinite(alpha1) and (0.5 <= alpha1 < 1.0):
                    C_eff = (Q1 ** (1.0 / alpha1)) / (R1 ** (1.0 - 1.0 / alpha1))
                    Cdl_per_cm2 = C_eff / GC_AREA_CM2
                else:
                    Cdl_per_cm2 = float("nan")
            else:
                Rs_area = float("nan")
                Rct_area = float("nan")
                Cdl_per_cm2 = float("nan")

            f_hz = pd.to_numeric(df[f_col], errors="coerce").to_numpy()
            idx_peak = int(np.nanargmax(z_im))
            f_peak = float(f_hz[idx_peak]) if np.isfinite(f_hz[idx_peak]) else float("nan")

            if rpm is not None and rpm > 0 and np.isfinite(D_calc):
                omega = 2 * np.pi * (rpm / 60.0)
                delta_cm = 1.61 * (D_calc ** (1 / 3)) * (NU_CMS2 ** (1 / 6)) * (omega ** (-0.5))
            else:
                delta_cm = float("nan")

            # k0 from PEIS (primary): k0 = RT / (n^2 F^2 C* Rct)
            if np.isfinite(Rct_area) and (Rct_area > 0) and (C_BULK_MOL_PER_CM3 > 0):
                k0_cm_s = float((R_GAS * T_K) / ((N_ELECTRONS ** 2) * (F_CONST ** 2) * C_BULK_MOL_PER_CM3 * Rct_area))
            else:
                k0_cm_s = float("nan")

            peis_rows.append(
                {
                    "label": label,
                    "rpm": rpm,
                    "Rs_area_Ohm_cm2": Rs_area,
                    "Rct_area_Ohm_cm2": Rct_area,
                    "Cdl_F_cm2": Cdl_per_cm2,
                    "f_peak_Hz": f_peak,
                    "delta_cm": delta_cm,
                    "D_used_cm2_s": D_calc if np.isfinite(D_calc) else float("nan"),
                    "k0_cm_s": k0_cm_s,
                }
            )

            z_fit_plot = None
            if fit.fitted_impedance is not None:
                z_model = np.asarray(fit.fitted_impedance)
                mask_model = np.isfinite(z_model.real) & np.isfinite(z_model.imag)
                if mask_model.any():
                    z_fit_plot = z_model[mask_model]

            fig_debug, ax_debug = plt.subplots(figsize=(5.0, 4.0))
            ax_debug.plot(z_re, z_im, 'o', ms=3, label='data')
            if z_fit_plot is not None:
                ax_debug.plot(z_fit_plot.real, -z_fit_plot.imag, '-', lw=1.5, label='fit')

            ax_debug.set_xlabel("Z' (ohm cm^2)")
            ax_debug.set_ylabel("-Z'' (ohm cm^2)")
            ax_debug.set_title(f"Debug PEIS Fit: {label}")
            ax_debug.text(0.05, 0.90, f'rchi2 = {fit.rchi2:.2e}\nN = {fit.n_points}', transform=ax_debug.transAxes, fontsize=10)
            ax_debug.legend()
            beautify_axes(ax_debug)
            ax_debug.set_aspect("equal", adjustable="box")

            debug_dir = Path("results/debugplots")
            debug_dir.mkdir(exist_ok=True, parents=True)
            debug_path = debug_dir / f"debug_PEIS_fit_{label_from_filename(path)}.png"
            fig_debug.savefig(debug_path, dpi=300, bbox_inches="tight", format='png', facecolor="white")
            plt.close(fig_debug)
            report.record_figure(debug_path)

        ax.set_xlabel(r"Z' [Ω cm$^2$]")
        ax.set_ylabel(r"-Z'' [Ω cm$^2$]")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        ax.set_aspect("equal", adjustable="box")
        report.record_figure(safe_save(fig, "T3.3_EIS_Nyquist.png"))

        if peis_rows:
            df_peis = pd.DataFrame(peis_rows)
            report.record_table(write_csv(df_peis, "T3.3_PEIS_summary.csv"))
            report.add_message(f"[PEIS] Processed {len(peis_rows)} PEIS spectra.")
            if "k0_cm_s" in df_peis.columns and df_peis["k0_cm_s"].notna().any():
                k0_med = float(np.nanmedian(df_peis["k0_cm_s"]))
                report.add_message(f"[PEIS] k0~{k0_med:.2e} cm/s (median across rotations)")
    else:
        message = "[PEIS] No PEIS files found; skipping PEIS plot."
        print(message)
        report.add_message(message)

    # --- Comparison vs Interface 1 (Task 3.2) ---
    try:
        # Interface 2 D from KL output (do not use global aggregator to avoid mixing with I1)
        D_I2 = float("nan")
        kl_path = Path("results") / "T3.3_KL_diffusion.csv"
        if kl_path.exists():
            df_kl = pd.read_csv(kl_path)
            if "D_KL_cm2_s" in df_kl.columns:
                D_I2 = float(np.nanmedian(pd.to_numeric(df_kl["D_KL_cm2_s"], errors="coerce")))
        k0_I2 = float("nan")
        peis_summary_path = Path("results") / "T3.3_PEIS_summary.csv"
        if peis_summary_path.exists():
            df_peis2 = pd.read_csv(peis_summary_path)
            if "k0_cm_s" in df_peis2.columns:
                k0_I2 = float(np.nanmedian(pd.to_numeric(df_peis2["k0_cm_s"], errors="coerce")))

        # Interface 1 values from saved results if available
        D_I1 = float("nan")
        k0_I1 = float("nan")
        rs_path = Path("results") / "T3.2_CV_randles_sevcik.csv"
        eis_path = Path("results") / "T3.2_EIS_summary.csv"
        if rs_path.exists():
            df_rs = pd.read_csv(rs_path)
            if "D_RandlesSevcik_cm2_s" in df_rs.columns:
                D_I1 = float(np.nanmedian(pd.to_numeric(df_rs["D_RandlesSevcik_cm2_s"], errors="coerce")))
        if eis_path.exists():
            df_e1 = pd.read_csv(eis_path)
            # Prefer EIS Warburg D if present
            if "D_Warburg_cm2_s" in df_e1.columns and pd.to_numeric(df_e1["D_Warburg_cm2_s"], errors="coerce").notna().any():
                D_I1 = float(np.nanmedian(pd.to_numeric(df_e1["D_Warburg_cm2_s"], errors="coerce"))) if not np.isfinite(D_I1) else D_I1
            if "k0_cm_s" in df_e1.columns:
                k0_I1 = float(np.nanmedian(pd.to_numeric(df_e1["k0_cm_s"], errors="coerce")))

        rows_compare: List[Dict[str, object]] = []
        if np.isfinite(D_I1):
            rows_compare.append({"metric": "D", "interface": 1, "method": "I1_RS/EIS", "value": D_I1, "units": "cm^2/s"})
        if np.isfinite(D_I2):
            rows_compare.append({"metric": "D", "interface": 2, "method": "I2_KL", "value": D_I2, "units": "cm^2/s"})
        if np.isfinite(k0_I1):
            rows_compare.append({"metric": "k0", "interface": 1, "method": "I1_PEIS", "value": k0_I1, "units": "cm/s"})
        if np.isfinite(k0_I2):
            rows_compare.append({"metric": "k0", "interface": 2, "method": "I2_PEIS", "value": k0_I2, "units": "cm/s"})

        if rows_compare:
            df_cmp = pd.DataFrame(rows_compare)
            report.record_table(write_csv(df_cmp, "T3.3_compare_interface1_vs_interface2.csv"))
            if np.isfinite(D_I1) and np.isfinite(D_I2) and (D_I1 > 0):
                diffD = (D_I2 - D_I1) / D_I1 * 100.0
                report.add_message(f"[Compare] D_I2 vs D_I1: {diffD:+.1f}%")
            if np.isfinite(k0_I1) and np.isfinite(k0_I2) and (k0_I1 > 0):
                diffk = (k0_I2 - k0_I1) / k0_I1 * 100.0
                report.add_message(f"[Compare] k0_I2 vs k0_I1: {diffk:+.1f}%")
    except Exception as exc:
        report.add_warning(f"[Compare] Failed to generate comparison: {exc}")

    return report
