# CoinCells/plot.py
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from scipy.signal import savgol_filter
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# -----------------------------
# Configuration
# -----------------------------
# Theoretical capacity for graphite anode (adjust if needed)
THEORETICAL_CAPACITY_MAH_G: float = 372.0

# Mass of active material per coin cell in mg, derived from Skript.md specs.
# Using graphite anode areal loading 6.3 mg/cm^2 and anode disc Ø14 mm:
# mass ≈ 6.3 mg/cm^2 × π × (0.7 cm)^2 ≈ 9.70 mg per cell.
# Keys are inferred from filenames: DS, FZ, HE
MASS_MG: Dict[str, float] = {
    "DS": 9.70,
    "FZ": 9.70,
    "HE": 9.70,
}

# Plot styling
plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False


# -----------------------------
# Paths
# -----------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"
FIG_DIR.mkdir(parents=True, exist_ok=True)
RES_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------
# Utilities
# -----------------------------
def extract_sample_id(stem: str) -> str:
    """
    Extract short sample ID from filename (e.g., DS, FZ, HE).
    """
    m = re.search(r"_([A-Z]{2})_", stem)
    if m:
        return m.group(1)
    # Fallback: end-of-stem token
    parts = re.split(r"[_\s]", stem)
    for token in parts[::-1]:
        if len(token) == 2 and token.isalpha() and token.isupper():
            return token
    return stem


def read_coin_csv(path: Path) -> pd.DataFrame:
    """
    Robust CSV reader for various lab exports.
    - Auto-detects delimiter
    - Normalizes column names
    - Tries to locate time, voltage, current, capacity, cycle
    """
    def detect_data_header_line(p: Path) -> Optional[int]:
        """
        Scan file for the line that starts the actual data table.
        Looks for a header containing 'Time /s' and enough comma-separated fields.
        Returns 0-based line index or None if not found.
        """
        for enc in ("utf-8", "latin-1"):
            try:
                with open(p, "r", encoding=enc, errors="ignore") as f:
                    for idx, line in enumerate(f):
                        s = line.strip().lower()
                        # Use a robust condition: header line contains 'time /s' and many commas
                        if ("time /s" in s) and (s.count(",") >= 10):
                            return idx
            except Exception:
                continue
        return None

    def _read(p: Path, skiprows: Optional[int] = None, header: Optional[int] = "infer") -> pd.DataFrame:
        try:
            return pd.read_csv(
                p,
                sep=",",
                engine="python",
                encoding="utf-8",
                skiprows=skiprows,
                header=header,
                skipinitialspace=True,
                on_bad_lines="skip",
            )
        except Exception:
            return pd.read_csv(
                p,
                sep=",",
                engine="python",
                encoding="latin-1",
                skiprows=skiprows,
                header=header,
                skipinitialspace=True,
                on_bad_lines="skip",
            )

    header_line = detect_data_header_line(path)
    if header_line is not None:
        df = _read(path, skiprows=header_line, header=0)
    else:
        df = _read(path)

    # Normalize column names
    def norm(c: str) -> str:
        c = c.strip()
        c = c.replace("(", "[").replace(")", "]")
        c = re.sub(r"\s+", "_", c)
        return c.lower()

    df.columns = [norm(c) for c in df.columns]

    # Attempt to coerce numeric where possible (without deprecated errors='ignore')
    for c in df.columns:
        if df[c].dtype == object:
            series_str = df[c].astype(str).str.replace(",", ".", regex=False)
            series_num = pd.to_numeric(series_str, errors="coerce")
            # Keep numeric conversion if at least half the entries are numbers
            if series_num.notna().mean() >= 0.5:
                df[c] = series_num
            else:
                df[c] = series_str

    # Identify columns
    col_time = first_existing(df, [
        "time_s", "time_[s]", "time_/s", "time", "test_time[s]", "elapsed_time[s]", "t_[s]", "record_time_[s]"
    ])
    col_volt = first_existing(df, [
        "voltage", "voltage_[v]", "ewe/v", "e_[v]", "v", "cell_voltage_[v]", "u/[v]", "u_/v"
    ])
    col_curr_a = first_existing(df, [
        "current_[a]", "i_[a]", "current(a)", "current", "i(a)", "i_/a", "current_/a"
    ])
    col_curr_ma = first_existing(df, [
        "current_[ma]", "i_[ma]", "current(ma)", "i(ma)", "i_/ma", "current_/ma"
    ])
    col_cyc = first_existing(df, [
        "cycle", "cycle_index", "cycle_number", "cycle#", "loop", "loop_number"
    ])
    col_step = first_existing(df, [
        "step", "step_index", "step_number", "stage", "state", "mode", "step_name"
    ])
    col_type = first_existing(df, [
        "type", "technique", "method"
    ])
    col_task = first_existing(df, [
        "task"
    ])
    col_q_chg = first_existing(df, [
        "charge_capacity_[mah]", "q_charge_[mah]", "charge_capacity_mah", "charge_cap_[mah]",
        "q_charge_/ah", "charge_capacity_[ah]", "charge_capacity_ah"
    ])
    col_q_dchg = first_existing(df, [
        "discharge_capacity_[mah]", "q_discharge_[mah]", "discharge_capacity_mah", "discharge_cap_[mah]",
        "q_discharge_/ah", "discharge_capacity_[ah]", "discharge_capacity_ah"
    ])
    col_q_total = first_existing(df, [
        "capacity_[mah]", "q_[mah]", "capacity_mah", "cap_[mah]",
        "capacity_/ah", "q_/ah", "capacity_[ah]", "capacity_ah"
    ])

    # Construct canonical columns
    out = pd.DataFrame(index=df.index)

    # time in seconds
    if col_time is not None:
        out["time_s"] = df[col_time].astype(float)
        # If time seems monotonically increasing and in hours, convert if values are small?
        # We leave as-is, assuming seconds; user exports typically seconds.
    else:
        # fabricate monotonic time if missing
        out["time_s"] = np.arange(len(df), dtype=float)

    # voltage in V
    if col_volt is None:
        raise ValueError(f"Could not find a voltage column in {path.name}")
    out["voltage_v"] = df[col_volt].astype(float)

    # current in A
    if col_curr_a is not None:
        out["current_a"] = df[col_curr_a].astype(float)
    elif col_curr_ma is not None:
        out["current_a"] = df[col_curr_ma].astype(float) / 1000.0
    else:
        # Try to infer if any column looks like current by unit hints
        cand = [c for c in df.columns if "current" in c or c.startswith("i_")]
        if cand:
            # assume A
            out["current_a"] = pd.to_numeric(df[cand[0]], errors="coerce").fillna(0.0)
        else:
            raise ValueError(f"Could not find a current column in {path.name}")

    # step/cycle metadata (optional)
    if col_cyc is not None:
        out["cycle_index"] = pd.to_numeric(df[col_cyc], errors="coerce")
    else:
        out["cycle_index"] = np.nan

    if col_step is not None:
        out["step_label"] = df[col_step].astype(str)
    else:
        out["step_label"] = ""

    if col_type is not None:
        out["type_label"] = df[col_type].astype(str).str.lower()
    else:
        out["type_label"] = ""

    if col_task is not None:
        out["task_label"] = df[col_task].astype(str).str.lower()
    else:
        out["task_label"] = ""

    # capacity (mAh) if present
    if col_q_total is not None:
        series = pd.to_numeric(df[col_q_total], errors="coerce")
        if ("ah" in col_q_total) and ("mah" not in col_q_total):
            series = series * 1000.0
        out["capacity_mAh_raw"] = series
    else:
        out["capacity_mAh_raw"] = np.nan

    if col_q_chg is not None:
        series = pd.to_numeric(df[col_q_chg], errors="coerce")
        if ("ah" in col_q_chg) and ("mah" not in col_q_chg):
            series = series * 1000.0
        out["charge_capacity_mAh_raw"] = series
    else:
        out["charge_capacity_mAh_raw"] = np.nan

    if col_q_dchg is not None:
        series = pd.to_numeric(df[col_q_dchg], errors="coerce")
        if ("ah" in col_q_dchg) and ("mah" not in col_q_dchg):
            series = series * 1000.0
        out["discharge_capacity_mAh_raw"] = series
    else:
        out["discharge_capacity_mAh_raw"] = np.nan

    # Derived differentials
    out["dt_s"] = np.diff(out["time_s"], prepend=out["time_s"].iloc[0])
    out["dE_dt"] = np.gradient(out["voltage_v"].values, out["time_s"].values, edge_order=1)

    return out


def first_existing(df: pd.DataFrame, names: List[str]) -> Optional[str]:
    for n in names:
        if n in df.columns:
            return n
    return None


def smooth(y: np.ndarray, window: int = 21, poly: int = 3) -> np.ndarray:
    n = len(y)
    if n < 5:
        return y
    if HAS_SCIPY and window < n and window % 2 == 1:
        try:
            return savgol_filter(y, window_length=window, polyorder=min(poly, window - 2), mode="interp")
        except Exception:
            pass
    # fallback rolling mean
    s = pd.Series(y, dtype=float).rolling(window=min(11, max(3, n//10)), center=True, min_periods=1).mean().values
    return s


# -----------------------------
# Segmentation and metrics
# -----------------------------
def identify_cc_segments(df: pd.DataFrame, i_threshold_a: float = 1e-4) -> pd.DataFrame:
    """
    Identify constant-current like segments by sign of current.
    Adds columns:
      - seg_id: contiguous block id
      - is_discharge: bool
    """
    i = df["current_a"].values
    sign = np.where(np.abs(i) < i_threshold_a, 0, np.sign(i))
    seg_id = np.zeros(len(df), dtype=int)
    sid = 0
    for k in range(1, len(df)):
        if sign[k] != sign[k-1]:
            sid += 1
        seg_id[k] = sid
    out = df.copy()
    out["seg_id"] = seg_id
    out["is_discharge"] = (sign < 0)
    out["is_charge"] = (sign > 0)
    return out


def integrate_capacity_per_segment(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each constant-sign segment, integrate current over time to get capacity (mAh).
    Also compute per-segment normalized capacity (mAh/g) once mass is provided later.
    """
    out = df.copy()
    out["cap_mAh_seg"] = np.nan

    for sid, g in out.groupby("seg_id", sort=True):
        dq_mAh = np.cumsum(g["current_a"].values * g["dt_s"].values) / 3600.0 * 1000.0
        out.loc[g.index, "cap_mAh_seg"] = dq_mAh - (dq_mAh[0] if len(dq_mAh) else 0.0)

    return out


def estimate_cycles_from_segments(df: pd.DataFrame) -> pd.DataFrame:
    """
    Estimate cycle index by pairing charge and discharge segments.
    Assign cycle_idx starting from 1.
    """
    out = df.copy()
    seg_meta = []
    for sid, g in out.groupby("seg_id", sort=True):
        if len(g) == 0:
            continue
        sgn = np.sign(np.nansum(g["current_a"].values))
        seg_meta.append((sid, int(np.sign(sgn))))

    # Pair a charge (+1) followed by discharge (-1) as one cycle
    cycle_idx_by_seg: Dict[int, int] = {}
    cycle = 0
    k = 0
    while k < len(seg_meta):
        sid, sgn = seg_meta[k]
        if sgn > 0:  # charge
            # find next discharge
            if k + 1 < len(seg_meta) and seg_meta[k+1][1] < 0:
                cycle += 1
                cycle_idx_by_seg[sid] = cycle
                cycle_idx_by_seg[seg_meta[k+1][0]] = cycle
                k += 2
                continue
            else:
                cycle += 1
                cycle_idx_by_seg[sid] = cycle
        elif sgn < 0:
            cycle += 1
            cycle_idx_by_seg[sid] = cycle
        k += 1

    out["cycle_est"] = out["seg_id"].map(cycle_idx_by_seg).fillna(0).astype(int)
    return out


def compute_ce_per_cycle(df: pd.DataFrame, mass_g: float) -> pd.DataFrame:
    """
    Compute Coulombic Efficiency per estimated cycle using per-segment integrated capacity.
    Returns a dataframe with cycle, Qc, Qd, CE, and estimated C-rate.
    """
    # Summarize per segment
    per_seg = df.groupby("seg_id").agg(
        cycle=("cycle_est", "first"),
        is_discharge=("is_discharge", "first"),
        mean_I_A=("current_a", "mean"),
        cap_end_mAh=("cap_mAh_seg", lambda x: float(np.nanmax(np.abs(x)) if len(x) else np.nan)),
    ).reset_index()

    # Aggregate into cycles
    records = []
    for cyc, g in per_seg.groupby("cycle", sort=True):
        if cyc == 0:
            continue
        # Guard against all-NaN or empty segments
        Qc_vals = g.loc[g["is_discharge"] == False, "cap_end_mAh"].dropna()
        Qd_vals = g.loc[g["is_discharge"] == True, "cap_end_mAh"].dropna()
        if Qc_vals.empty or Qd_vals.empty:
            continue
        Qc_mAh = Qc_vals.sum()
        Qd_mAh = Qd_vals.sum()
        # Normalize to mAh/g (capacity already in mAh). Do not scale by 1000 again.
        Qc_mAh_g = Qc_mAh / mass_g
        Qd_mAh_g = Qd_mAh / mass_g
        CE = (Qd_mAh / Qc_mAh * 100.0) if Qc_mAh > 0 else np.nan
        I_mean = g["mean_I_A"].abs().mean()
        c_rate = estimate_c_rate(I_mean, mass_g, THEORETICAL_CAPACITY_MAH_G)
        records.append(dict(
            cycle=int(cyc),
            Qc_mAh_g=Qc_mAh_g,
            Qd_mAh_g=Qd_mAh_g,
            CE_percent=CE,
            C_rate=c_rate,
        ))
    return pd.DataFrame.from_records(records).sort_values("cycle")


def estimate_c_rate(I_A: float, mass_g: float, q_theor_mAh_g: float) -> float:
    """
    C-rate estimate: C = I / (mass * Q_theor_Ah/g) where Q_theor_Ah/g = q_theor_mAh_g/1000
    """
    denom = mass_g * (q_theor_mAh_g / 1000.0)
    return float(I_A / denom) if denom > 0 else np.nan


# -----------------------------
# CV detection and plotting
# -----------------------------
def try_extract_cv_cycles(df: pd.DataFrame, max_cycles: int = 3) -> List[pd.DataFrame]:
    """
    Heuristic extraction of CV cycles:
      - If step labels contain 'cv', use those rows.
      - Else, try to detect voltage scan reversals (sign flips in dE/dt) and build cycles.
    Returns a list of dataframes (each a CV cycle with columns E, I).
    """
    data = df.copy()

    # Quick sanity: if current is nearly zero everywhere, this is not CV
    if np.nanmax(np.abs(data["current_a"].values)) < 1e-9:
        return []

    # Prefer explicit 'cv' in labels
    mask_cv = (
        data["step_label"].str.lower().str.contains("cv", na=False)
        | data["type_label"].str.contains("cv", na=False)
        | data["task_label"].str.contains("cv", na=False)
    )
    cv_df = data.loc[mask_cv].copy()
    # Only accept CV if we have a meaningful number of rows tagged as CV
    if len(cv_df) < 50:
        return []

    # Identify reversal points from dE/dt sign changes
    de = cv_df["dE_dt"].values
    sign = np.sign(de)
    sign[np.isnan(sign)] = 0
    flips = np.where(np.diff(sign) != 0)[0] + 1
    if len(flips) < 3:
        # not convincing CV; return empty
        return []

    # Build cycles between alternating vertices; assume 2 scans per cycle
    # Use windows of 2 consecutive flips for one scan, 4 flips for one full cycle
    cycles = []
    start = 0
    # Each full cycle is approx 2 segments (forward+reverse), i.e., 2 flips
    # Use 2 flips per cycle for robustness
    for k in range(0, min(len(flips) - 1, 2 * max_cycles), 2):
        end = flips[k + 1]
        cyc = cv_df.iloc[start:end+1].copy()
        # Clean small NaNs
        cyc = cyc.dropna(subset=["voltage_v", "current_a"])
        if len(cyc) > 5:
            cycles.append(cyc[["voltage_v", "current_a"]].rename(columns={
                "voltage_v": "E_V", "current_a": "I_A"
            }))
        start = end
        if len(cycles) >= max_cycles:
            break

    return cycles


# -----------------------------
# Plotting
# -----------------------------
def plot_cv_overlay(sample_id: str, cycles: List[pd.DataFrame], outpath: Path) -> None:
    if not cycles:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No CV detected", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = plt.cm.tab10.colors
    for i, cyc in enumerate(cycles, start=1):
        ax.plot(cyc["E_V"].values, cyc["I_A"].values * 1000.0, lw=1.5, color=colors[(i - 1) % 10], label=f"Cycle {i}")
    ax.set_xlabel("E (V)")
    ax.set_ylabel("I (mA)")
    ax.set_title(f"CV overlay (first 3 cycles) - {sample_id}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_potential_vs_capacity(sample_id: str, df: pd.DataFrame, mass_g: float, outpath: Path) -> None:
    """
    Overlay charge/discharge segments across cycles in one plot: E vs Q (mAh/g).
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = plt.cm.tab10.colors

    # Ensure normalized capacity per segment exists
    df = df.copy()
    # Convert integrated capacity (mAh) to specific capacity (mAh/g)
    df["cap_mAh_g_seg"] = df["cap_mAh_seg"] / mass_g

    labels_done: set = set()
    labelled_cycles: set = set()

    # Choose a compact set of cycles to show in the legend (max 8)
    cycles_present = sorted([int(c) for c in df["cycle_est"].dropna().unique() if int(c) > 0])
    if cycles_present:
        num_to_label = min(8, len(cycles_present))
        # Evenly spaced selection including first and last
        sel_idx = np.unique(np.linspace(0, len(cycles_present) - 1, num=num_to_label, dtype=int))
        selected_cycles = {cycles_present[i] for i in sel_idx}
    else:
        selected_cycles = set()
    for cyc, gcyc in df.groupby("cycle_est", sort=True):
        if cyc == 0:
            continue
        color = colors[(cyc - 1) % len(colors)]
        # plot each segment in this cycle
        for sid, gseg in gcyc.groupby("seg_id", sort=True):
            if gseg["cap_mAh_g_seg"].abs().max() < 1e-6:
                continue
            is_discharge = bool(gseg["is_discharge"].iloc[0])
            # Label only once per cycle (use discharge as representative)
            if is_discharge and cyc in selected_cycles and cyc not in labelled_cycles:
                lbl = f"Cycle {cyc}"
                labelled_cycles.add(cyc)
            else:
                lbl = None
            ax.plot(
                gseg["cap_mAh_g_seg"].abs().values,
                gseg["voltage_v"].values,
                color=color,
                lw=1.4 if is_discharge else 1.2,
                ls="-" if is_discharge else "--",
                alpha=0.95 if is_discharge else 0.6,
                label=lbl,
            )
            labels_done.add(lbl or f"{sid}-{cyc}")

    # Y-axis: exactly min..max of plotted data (no extra margins)
    try:
        mask_plot = (
            (df["cycle_est"] > 0)
            & df["cap_mAh_g_seg"].abs().gt(1e-6)
            & df["voltage_v"].notna()
        )
        y_min = float(df.loc[mask_plot, "voltage_v"].min()) - 0.1
        y_max = float(df.loc[mask_plot, "voltage_v"].max()) + 0.1
        if np.isfinite(y_min) and np.isfinite(y_max) and y_max > y_min:
            ax.set_ylim(y_min, y_max)
            ax.margins(y=0)
    except Exception:
        pass

    ax.set_xlabel("Capacity (mAh/g)")
    ax.set_ylabel("Potential (V)")
    ax.set_title(f"Potential vs Capacity - {sample_id}")
    # Nice x-limits
    try:
        qmax = float(np.nanmax(np.abs(df["cap_mAh_g_seg"].values)))
        ax.set_xlim(0, min(max(200.0, qmax * 1.05), 380.0))
    except Exception:
        pass

    # Legend outside (right): cycles + style key
    handles, labels = ax.get_legend_handles_labels()
    style_handles = [
        Line2D([0], [0], color="0.2", lw=1.6, ls="-", label="Discharge"),
        Line2D([0], [0], color="0.2", lw=1.2, ls="--", label="Charge"),
    ]
    handles_all = handles + style_handles
    labels_all = labels + [h.get_label() for h in style_handles]
    if handles_all:
        ax.legend(
            handles_all,
            labels_all,
            title="Cycles",
            fontsize=8,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            framealpha=0.95,
        )
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def plot_dqdE_cathodic(sample_id: str, df: pd.DataFrame, mass_g: float, outpath: Path) -> None:
    """
    Compute and plot dQ/dE vs E for cathodic (discharge) segments.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = plt.cm.tab10.colors
    color_idx = 0

    # Normalized capacity within segment
    df = df.copy()
    # Specific capacity in mAh/g for discharge segments
    df["cap_mAh_g_seg"] = df["cap_mAh_seg"] / mass_g

    for cyc, gcyc in df.groupby("cycle_est", sort=True):
        if cyc == 0:
            continue
        for sid, gseg in gcyc.groupby("seg_id", sort=True):
            if not bool(gseg["is_discharge"].iloc[0]):
                continue
            E = gseg["voltage_v"].values
            Q = gseg["cap_mAh_g_seg"].abs().values
            if len(E) < 5:
                continue
            # Smooth E and Q with Savitzky–Golay (fallback if SciPy unavailable)
            if HAS_SCIPY:
                n_pts = len(E)
                win = min(11100, n_pts if n_pts % 2 == 1 else n_pts - 1)
                if win >= 5:
                    poly = 3 if win > 5 else 2
                    Es = savgol_filter(E, window_length=win, polyorder=min(poly, win - 1), mode="interp")
                    Qs = savgol_filter(Q, window_length=win, polyorder=min(poly, win - 1), mode="interp")
                else:
                    Es = smooth(E, window=7)
                    Qs = smooth(Q, window=7)
            else:
                Es = smooth(E, window=11)
                Qs = smooth(Q, window=11)
            # Compute derivative dQ/dE with guards against tiny dE and spikes
            with np.errstate(divide='ignore', invalid='ignore'):
                dQ_dE = np.gradient(Qs, Es)
            dE_local = np.gradient(Es)
            tiny = max(1e-6, np.nanmedian(np.abs(dE_local)) * 0.01)
            # Mask regions where voltage spacing is too small (unstable derivative)
            dQ_dE = np.where(np.abs(dE_local) < tiny, np.nan, dQ_dE)
            # Clip extreme outliers robustly
            finite_mask = np.isfinite(dQ_dE)
            if finite_mask.sum() > 20:
                lo, hi = np.nanpercentile(dQ_dE[finite_mask], [1.0, 99.0])
                dQ_dE = np.clip(dQ_dE, lo, hi)
            # Lightly smooth the derivative to suppress residual spikes
            dQ_dE = smooth(dQ_dE, window=min(101, max(7, (len(dQ_dE)//50)*2 + 1)))
            ax.plot(Es, dQ_dE, lw=1.4, color=colors[color_idx % 10], label=f"Cycle {cyc}")
            color_idx += 1

    ax.set_xlabel("E (V)")
    ax.set_ylabel("dQ/dE (mAh g^-1 V^-1)")
    ax.set_title(f"dQ/dE (cathodic) - {sample_id}")
    #if color_idx > 0:
    #ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


# New: simple voltage vs time diagnostic plot
def plot_voltage_over_time(sample_id: str, df: pd.DataFrame, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    t = df["time_s"].values
    e = df["voltage_v"].values
    i = df["current_a"].values

    color_e = "tab:blue"
    color_i = "tab:orange"

    line_e, = ax.plot(t, e, lw=1.2, color=color_e, label="Voltage (V)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Potential (V)", color=color_e)
    ax.tick_params(axis='y', labelcolor=color_e)

    ax2 = ax.twinx()
    line_i, = ax2.plot(t, i, lw=1.0, color=color_i, alpha=0.8, label="Current (A)")
    ax2.set_ylabel("Current (A)", color=color_i)
    ax2.tick_params(axis='y', labelcolor=color_i)

    ax.set_title(f"Voltage and Current vs Time - {sample_id}")

    lines = [line_e, line_i]
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="best", fontsize=8)

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


# -----------------------------
# Runner
# -----------------------------
def process_file(csv_path: Path, ce_records: List[Dict]) -> None:
    stem = csv_path.stem
    sample_id = extract_sample_id(stem)
    mass_mg = MASS_MG.get(sample_id, np.nan)
    if not np.isfinite(mass_mg) or mass_mg <= 0:
        print(f"[WARN] Mass for sample '{sample_id}' missing. Update MASS_MG in plot.py. Using 10 mg placeholder.")
        mass_mg = 10.0
    mass_g = mass_mg / 1000.0

    print(f"[INFO] Processing {csv_path.name} (sample={sample_id}, mass={mass_mg} mg)")

    df = read_coin_csv(csv_path)
    # Identify CC-like segments and integrate capacity
    df = identify_cc_segments(df)
    df = integrate_capacity_per_segment(df)
    df = estimate_cycles_from_segments(df)

    # CE
    ce_df = compute_ce_per_cycle(df, mass_g)
    if not ce_df.empty:
        ce_df.insert(0, "sample", sample_id)
        ce_records.extend(ce_df.to_dict(orient="records"))

    # Plots
    # 0) Voltage vs time (diagnostic)
    vt_path = FIG_DIR / f"T3.1_V_vs_t_{csv_path.stem}.png"
    plot_voltage_over_time(sample_id, df, vt_path)

    # 1) CV overlay (first three cycles if detectable)
    cv_cycles = try_extract_cv_cycles(df, max_cycles=3)
    cv_path = FIG_DIR / f"T3.1_CV_overlay_{csv_path.stem}.png"
    plot_cv_overlay(sample_id, cv_cycles, cv_path)

    # 2) Potential vs Capacity (overlay cycles and regimes)
    pvc_path = FIG_DIR / f"T3.1_Potential_vs_Capacity_{csv_path.stem}.png"
    plot_potential_vs_capacity(sample_id, df, mass_g, pvc_path)

    # 3) dQ/dE (cathodic only)
    dqdE_path = FIG_DIR / f"T3.1_dQdE_{csv_path.stem}.png"
    plot_dqdE_cathodic(sample_id, df, mass_g, dqdE_path)


def main(argv: Optional[List[str]] = None) -> int:
    csvs = sorted(HERE.glob("*.csv"))
    if not csvs:
        print(f"[ERROR] No CSV files found in {HERE}")
        return 1

    ce_records: List[Dict] = []
    for p in csvs:
        try:
            process_file(p, ce_records)
        except Exception as e:
            print(f"[ERROR] Failed {p.name}: {e}")

    if ce_records:
        ce_df = pd.DataFrame.from_records(ce_records)
        ce_df = ce_df[["sample", "cycle", "C_rate", "Qc_mAh_g", "Qd_mAh_g", "CE_percent"]]
        ce_out = RES_DIR / "T3.1_CE.csv"
        ce_df.to_csv(ce_out, index=False)
        print(f"[INFO] Wrote CE table -> {ce_out}")

    print("[INFO] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())