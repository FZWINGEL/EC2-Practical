import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from modules.config import (
    DATA_DIR,
    GC_AREA_CM2,
    NU_CMS2,
    N_ELECTRONS,
    R_GAS,
    F_CONST,
    T_K,
    C_BULK_MOL_PER_CM3,
)
from modules.diffusion import get_calculated_diffusion, add_diffusion_coefficient
from modules.io_utils import read_biologic_table, write_csv
from modules.plotting import beautify_axes, safe_save
from modules.task32_eis import fit_with_impedance_py
from modules.task32_eis import ECFFit
from modules.utils import extract_rpm_from_name, get_column, label_from_filename


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


def analyze_task33_test():
    """Test version that uses ALL potentials for KL analysis"""
    print("=== Task 3.3 Test - Using ALL Potentials ===")
    
    lsv_dir = DATA_DIR / "Task 3.3 LSV"
    if not lsv_dir.exists():
        print("[Task 3.3] Directory not found; skipping.")
        return

    lsv_files = sorted([p for p in lsv_dir.glob("*_LSV_C01.txt") if "PEIS" not in p.name])
    curves: Dict[int, pd.DataFrame] = {}

    if lsv_files:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        plot_data = []  # Store plot data for sorting
        
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
        safe_save(fig, "T3.3_LSV_overlay_test.png")
        print(f"[LSV] Saved LSV overlay plot")
    else:
        print("[LSV] No LSV files found; skipping LSV plot.")

    if len(curves) >= 3:
        try:
            E_grid, j_by_rpm = _interpolate_to_common_grid(curves)
        except Exception as exc:
            print(f"[KL] Could not interpolate to common grid: {exc}")
            return

        rpms = np.array(sorted(j_by_rpm.keys()))
        omega = 2 * np.pi * (rpms / 60.0)
        inv_sqrt_omega = 1.0 / np.sqrt(omega)

        # MODIFIED: Use ALL potentials instead of filtering
        E_targets = np.linspace(E_grid.min(), E_grid.max(), 50)  # More points, full range
        kl_rows: List[KLPoint] = []
        jk_vs_E: List[Tuple[float, float]] = []
        kl_summ_rows: List[Dict[str, float]] = []
        accepted_E: List[float] = []

        print(f"[KL] Testing {len(E_targets)} potentials across full range...")

        for E0 in E_targets:
            j_at_E = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E0))] for rpm in rpms])
            mask = np.isfinite(j_at_E) & (np.abs(j_at_E) > 1e-6)
            if mask.sum() < 3:
                continue
            y = 1.0 / np.abs(j_at_E[mask])
            x = inv_sqrt_omega[mask]

            # MODIFIED: Skip plateau filtering, just check linearity
            if len(x) >= 3:
                A = np.vstack([x, np.ones_like(x)]).T
                slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                y_pred = slope * x + intercept
                ss_res = float(np.sum((y - y_pred) ** 2))
                ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
                r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                
                # MODIFIED: Lower R² threshold to include more points
                if (r2 >= 0.90) and np.isfinite(slope) and (slope > 0) and np.isfinite(intercept) and (intercept != 0):
                    jk = 1.0 / intercept
                    kl_rows.extend(
                        KLPoint(E_V=float(E0), inv_sqrt_omega_s05=float(xi), inv_j_cm2_A_inv=float(yi), jk_A_cm2=float(jk))
                        for xi, yi in zip(x, y)
                    )
                    jk_vs_E.append((E0, jk))
                    accepted_E.append(float(E0))

                    # Calculate diffusion coefficient
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

        print(f"[KL] Accepted {len(accepted_E)} potentials (vs 9 in original)")

        if jk_vs_E:
            df_kl = pd.DataFrame([point.__dict__ for point in kl_rows])
            write_csv(df_kl, "T3.3_KouteckyLevich_points_test.csv")
            print(f"[KL] Saved {len(kl_rows)} Koutecky-Levich points")

            df_jk = pd.DataFrame({"E_V": [p[0] for p in jk_vs_E], "j_k_A_cm2": [p[1] for p in jk_vs_E]}).sort_values("E_V")
            write_csv(df_jk, "T3.3_jk_vs_E_test.csv")
            print(f"[KL] Derived {len(df_jk)} kinetic-current values")

            # Save KL diffusion results
            if kl_summ_rows:
                df_klD = pd.DataFrame(kl_summ_rows).sort_values("E_V")
                write_csv(df_klD, "T3.3_KL_diffusion_test.csv")
                D_vals = df_klD["D_KL_cm2_s"].to_numpy()
                D_vals = D_vals[np.isfinite(D_vals) & (D_vals > 0)]
                if D_vals.size:
                    D_med = float(np.nanmedian(D_vals))
                    print(f"[KL] D_KL~{D_med:.2e} cm^2/s (median across E)")

            # Tafel analysis with ALL points
            print("\n=== Tafel Analysis (ALL Points) ===")
            
            # Estimate equilibrium potential
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
                idx = int(np.argmin(np.abs(df_jk["j_k_A_cm2"])))
                E_eq = float(df_jk["E_V"].iloc[idx])

            jk_data = df_jk.copy()
            jk_data["eta_V"] = jk_data["E_V"] - E_eq
            jk_data = jk_data[np.isfinite(jk_data["eta_V"]) & np.isfinite(jk_data["j_k_A_cm2"])].copy()
            jk_data["sign"] = np.sign(jk_data["j_k_A_cm2"])
            sign_to_use = 1.0 if (jk_data["sign"] > 0).sum() >= (jk_data["sign"] < 0).sum() else -1.0
            branch = jk_data[(jk_data["sign"] == sign_to_use) & (np.abs(jk_data["j_k_A_cm2"]) > 1e-8)].copy()

            print(f"[Tafel] Using {len(branch)} points for Tafel analysis")

            if len(branch) >= 3:
                branch["log10_jk"] = np.log10(np.abs(branch["j_k_A_cm2"]))
                
                # MODIFIED: Use ALL points for Tafel fit (no filtering)
                x = branch["log10_jk"].to_numpy()
                y = branch["eta_V"].to_numpy()
                A = np.vstack([x, np.ones_like(x)]).T
                slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                y_pred = slope * x + intercept
                ss_res = np.sum((y - y_pred) ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0

                log10_j0 = -intercept / slope
                j0 = 10 ** log10_j0

                figtf, axtf = plt.subplots(figsize=(5.4, 4.0))
                
                # Plot all points
                axtf.plot(branch["log10_jk"], branch["eta_V"], "o", color="blue", label=f"All {len(branch)} Points")
                
                # Fit line across the whole plot
                xlim = axtf.get_xlim()
                x_intersect = -intercept / slope  # When η=0
                x_min = min(xlim[0], x_intersect - 0.1)  # Extend left to include j0
                xfit = np.linspace(x_min, xlim[1], 100)
                axtf.plot(xfit, slope * xfit + intercept, "-", color="red", label="Fit")
                
                # Calculate intersection with y-axis (η=0)
                axtf.plot(x_intersect, 0, "ro", markersize=8, label=f"j₀: {j0:.2e} A/cm²")
                
                axtf.set_xlabel(r"log$_{10}$(|j$_k$|) [A cm$^{-2}$]")
                axtf.set_ylabel(r"$\eta$ [V]")
                beautify_axes(axtf)
                axtf.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
                safe_save(figtf, "T3.3_Tafel_test.png")

                # Save Tafel results
                k0_app = float(j0 / (N_ELECTRONS * F_CONST * C_BULK_MOL_PER_CM3)) if (N_ELECTRONS > 0 and F_CONST > 0 and C_BULK_MOL_PER_CM3 > 0) else float("nan")
                df_tafel = pd.DataFrame({
                    "slope_V_per_decade": [slope],
                    "intercept_V": [intercept],
                    "E_eq_est_V": [E_eq],
                    "j0_A_cm2": [float(j0)],
                    "k0_app_cm_s": [float(k0_app)],
                    "R2": [r2],
                    "branch_sign": [int(sign_to_use)],
                    "n_points": [len(branch)],
                })
                write_csv(df_tafel, "T3.3_Tafel_fit_test.csv")
                
                print(f"[Tafel] j0~{j0:.2e} A/cm^2; slope={slope:.3f} V/dec; R²={r2:.3f}")
                print(f"[Tafel] Used {len(branch)} points (vs 4 in original)")
            else:
                print("[Tafel] Not enough points for Tafel analysis")
    else:
        print("[KL] Need LSV data at >=3 rotation rates for Koutecky-Levich and Tafel.")


if __name__ == "__main__":
    analyze_task33_test()
