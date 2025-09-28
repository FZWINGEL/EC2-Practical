from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from modules.config import C_BULK_MOL_PER_CM3, GC_AREA_CM2, N_ELECTRONS
from modules.diffusion import add_diffusion_coefficient
from modules.io_utils import read_biologic_table, write_csv
from modules.plotting import beautify_axes, safe_save
from modules.utils import get_column, label_from_filename

from ..config import AnalysisConfig
from ..pipeline.reporting import TaskReport


@dataclass
class CVPeak:
    scan_rate_Vs: float
    ipa_A: float
    Epa_V: float
    ipc_A: float
    Epc_V: float
    delta_Ep_V: float
    ipa_ipc_ratio: float


def estimate_scan_rate_Vs(e_v: np.ndarray, t_s: Optional[np.ndarray]) -> Optional[float]:
    if t_s is None:
        return None
    e = np.asarray(e_v)
    t = np.asarray(t_s)
    if len(e) < 5 or len(t) < 5:
        return None
    dt = np.gradient(t)
    dt[dt == 0] = np.nan
    dEdt = np.gradient(e) / dt
    value = np.nanmedian(np.abs(dEdt))
    return float(value) if np.isfinite(value) else None


def _pick_last_cycle(df: pd.DataFrame) -> pd.DataFrame:
    cycle_col = get_column(df, ["cycle number", "cycle", "Cycle"])
    if cycle_col is None:
        return df
    cycle = pd.to_numeric(df[cycle_col], errors="coerce")
    if cycle.notna().any():
        last_cycle = int(cycle.max())
        return df[cycle == last_cycle]
    return df


def _find_cv_extrema(e_v: np.ndarray, i_A: np.ndarray) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    window = max(5, len(i_A) // 200)
    if window % 2 == 0:
        window += 1
    window = min(window, 51)
    smoothed = pd.Series(i_A).rolling(window, center=True, min_periods=1).mean().to_numpy()
    ipa = float(np.nanmax(smoothed))
    Epa = float(e_v[np.nanargmax(smoothed)])
    ipc = float(np.nanmin(smoothed))
    Epc = float(e_v[np.nanargmin(smoothed)])
    return (Epa, ipa), (Epc, ipc)


def _load_cv_table(path: Path, report: TaskReport) -> Optional[pd.DataFrame]:
    try:
        return read_biologic_table(path)
    except Exception as exc:  # pragma: no cover - IO guard
        message = f"[CV] Failed to read {path.name}: {exc}"
        print(message)
        report.add_warning(message)
        return None


def _fit_randles_sevcik(df_summary: pd.DataFrame, report: TaskReport) -> None:
    df_rs = df_summary.dropna(subset=["scan_rate_Vs"])
    if len(df_rs) < 2:
        message = "[CV] Need >=2 scan rates with time to estimate D via Randles-Sevcik."
        print(message)
        report.add_warning(message)
        return

    v = np.sqrt(df_rs["scan_rate_Vs"].to_numpy())
    ip = np.abs(df_rs[["ipa_A", "ipc_A"]]).max(axis=1).to_numpy()
    A_mat = np.vstack([v, np.ones_like(v)]).T
    slope, intercept = np.linalg.lstsq(A_mat, ip, rcond=None)[0]
    const = 2.69e5 * (N_ELECTRONS ** 1.5) * GC_AREA_CM2 * C_BULK_MOL_PER_CM3
    D_half = slope / const
    D_rs = (D_half ** 2)
    report.record_table(
        write_csv(
            pd.DataFrame(
                {
                    "slope_A_per_Vs_sqrt": [slope],
                    "intercept_A": [intercept],
                    "D_RandlesSevcik_cm2_s": [D_rs],
                }
            ),
            "T3.2_CV_randles_sevcik.csv",
        )
    )
    add_diffusion_coefficient(D_rs, "Randles-Sevcik (combined)")
    if np.isfinite(D_rs):
        report.add_message(f"[CV] Randles-Sevcik D~{D_rs:.2e} cm^2/s")

    v_sqrt = v
    ipa_A = np.abs(df_rs["ipa_A"].to_numpy())
    ipc_A = np.abs(df_rs["ipc_A"].to_numpy())

    A_fit = np.vstack([v_sqrt, np.ones_like(v_sqrt)]).T
    slope_a, intercept_a = np.linalg.lstsq(A_fit, ipa_A, rcond=None)[0]
    slope_c, intercept_c = np.linalg.lstsq(A_fit, ipc_A, rcond=None)[0]

    D_half_a = slope_a / const
    D_half_c = slope_c / const
    D_rs_a = float(D_half_a ** 2) if np.isfinite(D_half_a) else np.nan
    D_rs_c = float(D_half_c ** 2) if np.isfinite(D_half_c) else np.nan
    add_diffusion_coefficient(D_rs_a, "Randles-Sevcik (anodic)")
    add_diffusion_coefficient(D_rs_c, "Randles-Sevcik (cathodic)")

    report.record_table(
        write_csv(
            pd.DataFrame(
                {
                    "slope_anodic_A_per_Vs_sqrt": [float(slope_a)],
                    "intercept_anodic_A": [float(intercept_a)],
                    "D_RandlesSevcik_anodic_cm2_s": [D_rs_a],
                    "slope_cathodic_A_per_Vs_sqrt": [float(slope_c)],
                    "intercept_cathodic_A": [float(intercept_c)],
                    "D_RandlesSevcik_cathodic_cm2_s": [D_rs_c],
                }
            ),
            "T3.2_CV_randles_sevcik_fits.csv",
        )
    )

    fig_rs, ax_rs = plt.subplots(figsize=(6.0, 4.2))
    ax_rs.plot(v_sqrt, ipa_A * 1e3, "o", label="|i_pa| (anodic)")
    ax_rs.plot(v_sqrt, ipc_A * 1e3, "s", label="|i_pc| (cathodic)")

    xfit = np.linspace(np.nanmin(v_sqrt), np.nanmax(v_sqrt), 100)
    yfit_a_mA = (slope_a * xfit + intercept_a) * 1e3
    yfit_c_mA = (slope_c * xfit + intercept_c) * 1e3
    ax_rs.plot(xfit, yfit_a_mA, "-", label=f"fit anodic; D~{D_rs_a:.2e} cm^2/s")
    ax_rs.plot(xfit, yfit_c_mA, "--", label=f"fit cathodic; D~{D_rs_c:.2e} cm^2/s")

    ax_rs.set_xlabel(r"$\sqrt{v}$ (V$^{1/2}$ s$^{-1/2}$)")
    ax_rs.set_ylabel(r"$i_p$ (mA)")
    ax_rs.set_title("Task 3.2 - Randles-Sevcik (peak current vs sqrt(v))")
    ax_rs.legend(title="Peaks and fits", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax_rs)
    report.record_figure(safe_save(fig_rs, "T3.2_CV_RandlesSevcik.png"))


def run_cv_analysis(config: AnalysisConfig) -> TaskReport:
    report = TaskReport(name="Task 3.2 - CV")
    cv_dir = config.data_dir / "Task 3.2 CV"
    files = sorted(cv_dir.glob("*_CV_C01.txt"))
    if not files:
        message = "[CV] No CV files found; skipping."
        print(message)
        report.add_message(message)
        return report

    rows: List[CVPeak] = []
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for path in files:
        df = _load_cv_table(path, report)
        if df is None:
            continue

        df_use = _pick_last_cycle(df)
        e_col = get_column(df_use, ["Ewe/V"]) or "Ewe/V"
        i_col = get_column(df_use, ["<I>/mA", "I/mA"]) or "<I>/mA"
        t_col = get_column(df_use, ["time/s", "Time/s"])
        if e_col not in df_use or i_col not in df_use:
            message = f"[CV] Missing columns in {path.name}; skipping."
            print(message)
            report.add_warning(message)
            continue

        e_v = pd.to_numeric(df_use[e_col], errors="coerce").to_numpy()
        i_mA = pd.to_numeric(df_use[i_col], errors="coerce").to_numpy()
        i_A = i_mA / 1000.0
        j_mA_cm2 = i_mA / GC_AREA_CM2

        scan_rate = None
        if t_col:
            t_s = pd.to_numeric(df_use[t_col], errors="coerce").to_numpy()
            scan_rate = estimate_scan_rate_Vs(e_v, t_s)
        if scan_rate is not None and np.isfinite(scan_rate):
            legend_label = f"v~{scan_rate*1000:.0f} mV/s"
        else:
            legend_label = label_from_filename(path)

        ax.plot(e_v, j_mA_cm2, lw=1.1, label=legend_label)

        (Epa, ipa_A), (Epc, ipc_A) = _find_cv_extrema(e_v, i_A)
        delta_ep = Epa - Epc
        ratio = abs(ipa_A) / abs(ipc_A) if ipc_A != 0 else np.nan

        rows.append(
            CVPeak(
                scan_rate_Vs=scan_rate if scan_rate is not None else float("nan"),
                ipa_A=float(ipa_A),
                Epa_V=float(Epa),
                ipc_A=float(ipc_A),
                Epc_V=float(Epc),
                delta_Ep_V=float(delta_ep),
                ipa_ipc_ratio=float(ratio),
            )
        )

    ax.set_xlabel("E (V vs Ag/AgCl)")
    ax.set_ylabel("j (mA cm$^{-2}$)")
    ax.set_title("Task 3.2 - CVs (last cycle per file)")
    ax.legend(title="Scan rate (approx.)", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    report.record_figure(safe_save(fig, "T3.2_CV_overlay.png"))

    if not rows:
        report.add_warning("[CV] No valid CV curves processed.")
        return report

    df_summary = pd.DataFrame([peak.__dict__ for peak in rows])
    report.record_table(write_csv(df_summary, "T3.2_CV_peaks.csv"))

    _fit_randles_sevcik(df_summary, report)

    return report


__all__ = ["run_cv_analysis", "estimate_scan_rate_Vs", "CVPeak"]


