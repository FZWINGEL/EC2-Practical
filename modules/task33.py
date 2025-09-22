import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from .config import DATA_DIR, GC_AREA_CM2, NU_CMS2, N_ELECTRONS
from .diffusion import get_calculated_diffusion
from .io_utils import read_biologic_table, write_csv
from .plotting import beautify_axes, safe_save
from .task32_eis import _extract_eis_parameters_both
from .task32_eis import fit_with_impedance_py
from .task32_eis import ECFFit
from .utils import extract_rpm_from_name, get_column, label_from_filename


@dataclass
class KLPoint:
    E_V: float
    inv_sqrt_omega_s05: float
    inv_j_cm2_A_inv: float


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


def analyze_task33_lsv_and_eis() -> None:
    lsv_dir = DATA_DIR / "Task 3.3 LSV"
    if not lsv_dir.exists():
        print("[Task 3.3] Directory not found; skipping.")
        return

    lsv_files = sorted([p for p in lsv_dir.glob("*_LSV_C01.txt") if "PEIS" not in p.name])
    curves: Dict[int, pd.DataFrame] = {}

    if lsv_files:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        for path in lsv_files:
            try:
                df = read_biologic_table(path)
            except Exception as exc:
                print(f"[LSV] Failed to read {path.name}: {exc}")
                continue

            e_col = get_column(df, ["Ewe/V"]) or "Ewe/V"
            i_col = get_column(df, ["<I>/mA", "I/mA"]) or "<I>/mA"
            if not ({e_col, i_col} <= set(df.columns)):
                print(f"[LSV] Missing columns in {path.name}; skipping.")
                continue

            e_v = pd.to_numeric(df[e_col], errors="coerce").to_numpy()
            j_mA_cm2 = pd.to_numeric(df[i_col], errors="coerce").to_numpy() / GC_AREA_CM2
            rpm = extract_rpm_from_name(path)
            label = f"{rpm} rpm" if rpm is not None else label_from_filename(path)
            ax.plot(e_v, j_mA_cm2, lw=1.0, label=label)

            if rpm is not None:
                df_use = df[[e_col, i_col]].copy()
                df_use.columns = ["Ewe/V", i_col]
                if i_col != "<I>/mA":
                    df_use["<I>/mA"] = df_use[i_col]
                curves[rpm] = df_use

        ax.set_xlabel("E (V vs Ag/AgCl)")
        ax.set_ylabel("j (mA cm$^{-2}$)")
        ax.set_title("Task 3.3 - LSVs at various rotation rates")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        safe_save(fig, "T3.3_LSV_overlay.png")
    else:
        print("[LSV] No LSV files found; skipping LSV plot.")

    if len(curves) >= 3:
        try:
            E_grid, j_by_rpm = _interpolate_to_common_grid(curves)
        except Exception as exc:
            print(f"[KL] Could not interpolate to common grid: {exc}")
            E_grid, j_by_rpm = None, None

        if E_grid is not None and j_by_rpm is not None:
            rpms = np.array(sorted(j_by_rpm.keys()))
            omega = 2 * np.pi * (rpms / 60.0)
            inv_sqrt_omega = 1.0 / np.sqrt(omega)

            E_targets = np.linspace(E_grid.min() + 0.05 * np.ptp(E_grid), E_grid.max() - 0.05 * np.ptp(E_grid), 8)
            kl_rows: List[KLPoint] = []
            jk_vs_E: List[Tuple[float, float]] = []

            for E0 in E_targets:
                j_at_E = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E0))] for rpm in rpms])
                mask = np.isfinite(j_at_E) & (np.abs(j_at_E) > 1e-6)
                if mask.sum() < 3:
                    continue
                y = 1.0 / j_at_E[mask]
                x = inv_sqrt_omega[mask]
                A = np.vstack([x, np.ones_like(x)]).T
                slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                jk = 1.0 / intercept
                kl_rows.extend(
                    KLPoint(E_V=float(E0), inv_sqrt_omega_s05=float(xi), inv_j_cm2_A_inv=float(yi))
                    for xi, yi in zip(x, y)
                )
                jk_vs_E.append((E0, jk))

            if jk_vs_E:
                df_kl = pd.DataFrame([point.__dict__ for point in kl_rows])
                write_csv(df_kl, "T3.3_KouteckyLevich_points.csv")

                df_jk = pd.DataFrame({"E_V": [p[0] for p in jk_vs_E], "j_k_A_cm2": [p[1] for p in jk_vs_E]}).sort_values("E_V")
                write_csv(df_jk, "T3.3_jk_vs_E.csv")

                mid = len(jk_vs_E) // 2
                if mid < len(jk_vs_E):
                    E_mid = jk_vs_E[mid][0]
                    j_mid = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E_mid))] for rpm in rpms])
                    mask = np.isfinite(j_mid) & (np.abs(j_mid) > 1e-6)
                    x = inv_sqrt_omega[mask]
                    y = 1.0 / j_mid[mask]
                    A = np.vstack([x, np.ones_like(x)]).T
                    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                    figkl, axkl = plt.subplots(figsize=(5.4, 4.0))
                    axkl.plot(x, y, "o")
                    axkl.plot(x, slope * x + intercept, "-")
                    axkl.set_xlabel("1/sqrt(omega) (s^0.5)")
                    axkl.set_ylabel("1/j (cm^2 A^-1)")
                    axkl.set_title(f"Koutecky-Levich at E = {E_mid:.3f} V")
                    beautify_axes(axkl)
                    safe_save(figkl, "T3.3_KL_example.png")

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
                    candidates = [(0.03, 0.15), (0.02, 0.20), (0.02, 0.25), (0.015, 0.30), (0.01, 0.35)]
                    best = None
                    for eta_min, eta_max in candidates:
                        selected = (np.abs(branch["eta_V"]) >= eta_min) & (np.abs(branch["eta_V"]) <= eta_max)
                        if selected.sum() >= 3:
                            x = branch.loc[selected, "log10_jk"].to_numpy()
                            y = branch.loc[selected, "eta_V"].to_numpy()
                            A = np.vstack([x, np.ones_like(x)]).T
                            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                            y_pred = slope * x + intercept
                            ss_res = np.sum((y - y_pred) ** 2)
                            ss_tot = np.sum((y - np.mean(y)) ** 2)
                            r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                            if best is None or r2 > best["r2"]:
                                best = {
                                    "eta_min": eta_min,
                                    "eta_max": eta_max,
                                    "slope": float(slope),
                                    "intercept": float(intercept),
                                    "r2": float(r2),
                                    "x": x,
                                    "y": y,
                                }

                    if best is None:
                        selected = np.abs(branch["eta_V"]) >= 0.02
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
                                "eta_min": 0.02,
                                "eta_max": float(np.abs(branch["eta_V"]).max()),
                                "slope": float(slope),
                                "intercept": float(intercept),
                                "r2": float(r2),
                                "x": x,
                                "y": y,
                            }

                    if best is not None:
                        log10_j0 = -best["intercept"] / best["slope"]
                        j0 = 10 ** log10_j0

                        figtf, axtf = plt.subplots(figsize=(5.4, 4.0))
                        axtf.plot(best["x"], best["y"], "o", label="data")
                        xfit = np.linspace(np.min(best["x"]), np.max(best["x"]), 100)
                        axtf.plot(xfit, best["slope"] * xfit + best["intercept"], "-", label=f"fit; R^2={best['r2']:.2f}; j0~{j0:.2e} A/cm^2")
                        axtf.set_xlabel("log10(|j_k| / A cm^-2)")
                        axtf.set_ylabel("eta (V)")
                        axtf.set_title("Task 3.3 - Tafel from kinetic currents")
                        beautify_axes(axtf)
                        axtf.legend(title="Data and fit", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
                        safe_save(figtf, "T3.3_Tafel.png")

                        write_csv(
                            pd.DataFrame(
                                {
                                    "slope_V_per_decade": [best["slope"]],
                                    "intercept_V": [best["intercept"]],
                                    "E_eq_est_V": [E_eq],
                                    "j0_A_cm2": [float(j0)],
                                    "R2": [best["r2"]],
                                    "branch_sign": [int(sign_to_use)],
                                    "eta_min_V": [best["eta_min"]],
                                    "eta_max_V": [best["eta_max"]],
                                }
                            ),
                            "T3.3_Tafel_fit.csv",
                        )
                    else:
                        print("[Tafel] Not enough points in any reasonable region to fit.")
                else:
                    print("[Tafel] Not enough points on a single branch to fit.")
    else:
        print("[KL] Need LSV data at >=3 rotation rates for Koutecky-Levich and Tafel.")

    peis_files = sorted([p for p in lsv_dir.glob("*_PEIS_C01.txt")])
    if peis_files:
        D_calc = get_calculated_diffusion()
        if not np.isfinite(D_calc):
            print("[PEIS] No calculated D available yet (RS/Sand/Warburg). delta will be NaN. Run Task 3.2 first in the same run.")
        fig, ax = plt.subplots(figsize=(6.0, 4.8))
        peis_rows = []
        for path in peis_files:
            try:
                df = read_biologic_table(path)
            except Exception as exc:
                print(f"[PEIS] Failed to read {path.name}: {exc}")
                continue

            re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
            im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
            if not ({re_col, im_col} <= set(df.columns)):
                print(f"[PEIS] Missing impedance columns in {path.name}; skipping.")
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
                fit = fit_with_impedance_py(df, circuit_str='R0-p(R1,C1)-p(R2,C2)')
            except Exception as exc:
                print(f"[PEIS] Fit failed for {path.name}: {exc}")
                fit = ECFFit(model="", Rs=float("nan"), Rct=float("nan"), Cdl=None, sigma=None, Q=None, alpha=None, rchi2=float("nan"), n_points=0)

            if fit.Rs is not None and hasattr(fit, 'parameters_') and len(fit.parameters_) >= 5:
                Rs_area = fit.parameters_[0] * GC_AREA_CM2
                Rct_area = (fit.parameters_[1] + fit.parameters_[3]) * GC_AREA_CM2  # R1 + R2
                Cdl_per_cm2 = fit.parameters_[2] / GC_AREA_CM2  # Use first C as main Cdl
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
                }
            )

            z_fit_plot = None
            if fit.fitted_impedance is not None:
                z_model = np.asarray(fit.fitted_impedance)
                mask_model = np.isfinite(z_model.real) & np.isfinite(z_model.imag)
                if mask_model.any():
                    z_fit_plot = z_model[mask_model]
                    ax.plot(z_fit_plot.real, -z_fit_plot.imag, lw=1.2, alpha=0.9, label=f"{label} (fit)")

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
            fig_debug.savefig(debug_dir / f"debug_PEIS_fit_{label_from_filename(path)}.png", dpi=150, bbox_inches="tight")
            plt.close(fig_debug)

        ax.set_xlabel("Z' (ohm cm^2)")
        ax.set_ylabel("-Z'' (ohm cm^2)")
        ax.set_title("Task 3.3 - PEIS Nyquist (area-normalized)")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        ax.set_aspect("equal", adjustable="box")
        safe_save(fig, "T3.3_EIS_Nyquist.png")

        if peis_rows:
            write_csv(pd.DataFrame(peis_rows), "T3.3_PEIS_summary.csv")
    else:
        print("[PEIS] No PEIS files found; skipping PEIS plot.")

