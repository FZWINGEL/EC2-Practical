from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from modules.config import C_BULK_mM, GC_AREA_CM2, N_ELECTRONS
from modules.diffusion import add_diffusion_coefficient
from modules.io_utils import read_biologic_table, write_csv
from modules.plotting import beautify_axes, safe_save
from modules.utils import get_column, label_from_filename

from ..config import AnalysisConfig
from ..pipeline.reporting import TaskReport


@dataclass
class CPResult:
    label: str
    current_mA: float
    tau_s: float
    E_tau_over_4_V: float
    D_Sand_cm2_s: float


def _segment_by_current_plateaus(t_s: np.ndarray, i_mA: Optional[np.ndarray]) -> List[Tuple[int, int]]:
    n = len(t_s)
    if i_mA is None or len(i_mA) != n or n < 30:
        return [(0, n - 1)]

    window = max(7, n // 300)
    if window % 2 == 0:
        window += 1
    window = min(window, 101)
    smoothed = pd.Series(i_mA).rolling(window, center=True, min_periods=1).median().to_numpy()

    diff = np.abs(np.diff(smoothed))
    denom = np.maximum(np.maximum(np.abs(smoothed[:-1]), np.abs(smoothed[1:])), 1e-6)
    rel = diff / denom

    median_abs = float(np.nanmedian(np.abs(smoothed))) if np.isfinite(np.nanmedian(np.abs(smoothed))) else 0.0
    abs_guard = max(0.0002, 0.1 * median_abs)
    change_mask = (rel > 0.20) | (diff > abs_guard)
    change_idxs = np.where(change_mask)[0]

    if len(change_idxs) == 0:
        return [(0, n - 1)]

    gap = max(5, n // 200)
    groups: List[List[int]] = []
    current: List[int] = [int(change_idxs[0])]
    for index in change_idxs[1:]:
        if int(index) - current[-1] <= gap:
            current.append(int(index))
        else:
            groups.append(current)
            current = [int(index)]
    groups.append(current)

    cuts = [grp[-1] + 1 for grp in groups]
    boundaries = [0] + [c for c in cuts if 0 < c < n - 1] + [n - 1]
    boundaries = sorted(set(boundaries))

    min_pts = max(20, n // 200)
    min_dt = 1.0
    segments: List[Tuple[int, int]] = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if end <= start:
            continue
        if (end - start + 1) < min_pts:
            continue
        if (float(t_s[end]) - float(t_s[start])) < min_dt:
            continue
        segments.append((int(start), int(end)))

    return segments if segments else [(0, n - 1)]


def _estimate_tau_s(t_s: np.ndarray, e_v: np.ndarray) -> float:
    n = len(t_s)
    if n < 5:
        return float("nan")

    window = max(11, n // 200)
    if window % 2 == 0:
        window += 1
    window = min(window, 101)
    smoothed = pd.Series(e_v).rolling(window, center=True, min_periods=1).mean().to_numpy()

    e0, e1 = float(smoothed[0]), float(smoothed[-1])
    emin, emax = float(np.nanmin(smoothed)), float(np.nanmax(smoothed))
    if not all(np.isfinite(val) for val in [e0, e1, emin, emax]):
        return float("nan")
    rising = e1 >= e0
    if emax - emin <= 1e-6:
        return float("nan")
    target = (emin + 0.9 * (emax - emin)) if rising else (emin + 0.1 * (emax - emin))

    lo = int(0.05 * n)
    idxs = np.arange(lo, n)
    if rising:
        crossed = idxs[smoothed[lo:] >= target]
    else:
        crossed = idxs[smoothed[lo:] <= target]
    if len(crossed):
        return float(t_s[int(crossed[0])])

    de_dt = np.gradient(smoothed, t_s)
    d2e_dt2 = np.gradient(de_dt, t_s)
    lo2 = int(0.40 * n)
    hi2 = int(0.98 * n)
    if hi2 <= lo2:
        lo2, hi2 = 0, n
    idx = np.argmax(np.abs(d2e_dt2[lo2:hi2])) + lo2
    tau = float(t_s[int(idx)])
    return tau if (np.isfinite(tau) and tau > 0) else float("nan")


def _fallback_tau(t_seg: np.ndarray, e_seg: np.ndarray) -> float:
    try:
        de_dt = np.gradient(pd.Series(e_seg).rolling(11, center=True, min_periods=1).mean().to_numpy(), t_seg)
        early = de_dt[: max(5, len(de_dt) // 10)]
        threshold = 3.0 * float(np.nanmedian(np.abs(early))) if np.isfinite(np.nanmedian(np.abs(early))) else float("nan")
        idx = int(np.argmax(de_dt > threshold)) if np.any(de_dt > threshold) else (len(t_seg) - 1)
        return float(t_seg[idx])
    except Exception:
        return float("nan")


def _load_cp_table(path: Path, report: TaskReport) -> Optional[pd.DataFrame]:
    try:
        return read_biologic_table(path)
    except Exception as exc:  # pragma: no cover - IO guard
        message = f"[CP] Failed to read {path.name}: {exc}"
        print(message)
        report.add_warning(message)
        return None


def run_cp_analysis(config: AnalysisConfig) -> TaskReport:
    report = TaskReport(name="Task 3.2 - CP")
    cp_dir = config.data_dir / "Task 3.2 CP"
    files = sorted(cp_dir.glob("*_CP_C01.txt"))
    if not files:
        message = "[CP] No CP files found; skipping."
        print(message)
        report.add_message(message)
        return report

    results: List[CPResult] = []
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for path in files:
        df = _load_cp_table(path, report)
        if df is None:
            continue

        t_col = get_column(df, ["time/s", "Time/s"]) or "time/s"
        e_col = get_column(df, ["Ewe/V"]) or "Ewe/V"
        i_col = get_column(df, ["<I>/mA", "I/mA"])
        if t_col not in df or e_col not in df:
            message = f"[CP] Missing time/E columns in {path.name}; skipping."
            print(message)
            report.add_warning(message)
            continue

        t_s = pd.to_numeric(df[t_col], errors="coerce").to_numpy()
        e_v = pd.to_numeric(df[e_col], errors="coerce").to_numpy()
        ax.plot(t_s, e_v, lw=1.0, label=label_from_filename(path))

        i_series = pd.to_numeric(df[i_col], errors="coerce").to_numpy() if (i_col and i_col in df) else None
        segments = _segment_by_current_plateaus(t_s, i_series)

        file_label = label_from_filename(path)
        for seg_idx, (start, end) in enumerate(segments, start=1):
            t_seg = t_s[start : end + 1] - float(t_s[start])
            e_seg = e_v[start : end + 1]
            if i_series is not None:
                i_seg_mA = pd.Series(i_series[start : end + 1]).rolling(11, center=True, min_periods=1).median().to_numpy()
                i_applied_mA = float(np.nanmedian(i_seg_mA))
            else:
                i_applied_mA = float("nan")

            if not np.isfinite(i_applied_mA) or abs(i_applied_mA) < 0.0002:
                continue

            try:
                tau = _estimate_tau_s(t_seg, e_seg)
            except Exception:
                tau = float("nan")

            if not np.isfinite(tau) or tau <= 0:
                tau = _fallback_tau(t_seg, e_seg)

            if np.isfinite(tau) and tau > 0:
                E_tau4 = float(np.interp(tau / 4.0, t_seg, e_seg))
            else:
                E_tau4 = float("nan")

            if np.isfinite(i_applied_mA) and np.isfinite(tau) and tau > 0:
                D_half = (abs(i_applied_mA) * math.sqrt(tau)) / (85.5 * N_ELECTRONS * GC_AREA_CM2 * C_BULK_mM)
                D_sand = float(D_half ** 2)
                add_diffusion_coefficient(D_sand, "Sand")
                if np.isfinite(D_sand):
                    report.add_message(f"[CP] Segment {seg_idx} D~{D_sand:.2e} cm^2/s")
            else:
                D_sand = float("nan")

            seg_label = f"{file_label} [seg {seg_idx}; I~{i_applied_mA * 1000:.0f} uA]"
            results.append(
                CPResult(
                    label=seg_label,
                    current_mA=i_applied_mA,
                    tau_s=tau,
                    E_tau_over_4_V=E_tau4,
                    D_Sand_cm2_s=D_sand,
                )
            )

    ax.set_xlabel("t (s)")
    ax.set_ylabel("E (V vs Ag/AgCl)")
    ax.set_title("Task 3.2 - Chronopotentiometry (E vs t)")
    ax.legend(title="Applied current", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    report.record_figure(safe_save(fig, "T3.2_CP_transients.png"))

    if results:
        df_results = pd.DataFrame([res.__dict__ for res in results])
        report.record_table(write_csv(df_results, "T3.2_CP_summary.csv"))
    else:
        report.add_warning("[CP] No usable Sand segments detected.")

    return report


__all__ = ["run_cp_analysis", "CPResult"]


