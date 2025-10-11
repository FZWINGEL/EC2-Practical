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
from matplotlib.ticker import AutoMinorLocator

try:
    from scipy.signal import savgol_filter
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# -----------------------------
# Configuration
# -----------------------------
# Theoretical capacity for graphite anode (adjust if needed)
# Updated to match expected C-rate currents: 0.1C=0.2376mA, 1C=2.3756mA, 2C=4.752mA
THEORETICAL_CAPACITY_MAH_G: float = 244.93

# Mass of active material per coin cell in mg, derived from Skript.md specs.
# Using graphite anode areal loading 6.3 mg/cm^2 and anode disc Ø14 mm:
# mass ≈ 6.3 mg/cm^2 × π × (0.7 cm)^2 ≈ 9.70 mg per cell.
# Keys are inferred from filenames: DS, FZ, HE
MASS_MG: Dict[str, float] = {
    "DS": 9.70,
    "FZ": 9.70,
    "HE": 9.70,
}

# dQdE smoothing configuration
class SmoothingConfig:
    """Configuration for dQdE smoothing methods."""
    
    # Available smoothing methods
    NONE = "none"           # No smoothing
    SAVITZKY_GOLAY = "savgol"  # Savitzky-Golay filter
    MOVING_AVERAGE = "moving_avg"  # Simple moving average
    
    def __init__(self, method: str = SAVITZKY_GOLAY, window_length: int = 21, polyorder: int = 3):
        """
        Initialize smoothing configuration.
        
        Args:
            method: Smoothing method ('none', 'savgol', 'moving_avg')
            window_length: Window length for smoothing (must be odd for Savitzky-Golay)
            polyorder: Polynomial order for Savitzky-Golay (ignored for other methods)
        """
        self.method = method
        self.window_length = window_length
        self.polyorder = polyorder
        
        # Ensure window_length is odd for Savitzky-Golay
        if method == self.SAVITZKY_GOLAY and window_length % 2 == 0:
            self.window_length = window_length + 1

# Default smoothing configuration - easily change this to modify behavior
DEFAULT_SMOOTHING = SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=21, polyorder=3)

# Easy configuration examples - uncomment one of these to change smoothing behavior:
# DEFAULT_SMOOTHING = SmoothingConfig(method=SmoothingConfig.NONE)  # No smoothing
# DEFAULT_SMOOTHING = SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=15, polyorder=2)  # Lighter Savitzky-Golay
# DEFAULT_SMOOTHING = SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=31, polyorder=4)  # Heavier Savitzky-Golay
# DEFAULT_SMOOTHING = SmoothingConfig(method=SmoothingConfig.MOVING_AVERAGE, window_length=15)  # Moving average

# Plot styling
plt.rcParams["figure.dpi"] = 300  # high DPI for PNG
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["text.usetex"] = False  # Use plain text, no LaTeX


# -----------------------------
# Paths
# -----------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"
PLOT_EXT = ".png"
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


def smooth(y: np.ndarray, config: SmoothingConfig = None) -> np.ndarray:
    """
    Apply smoothing to data based on configuration.
    
    Args:
        y: Input data array
        config: SmoothingConfig object specifying method and parameters
        
    Returns:
        Smoothed data array
    """
    if config is None:
        config = DEFAULT_SMOOTHING
        
    n = len(y)
    if n < 5:
        return y
    
    if config.method == SmoothingConfig.NONE:
        return y
    
    elif config.method == SmoothingConfig.SAVITZKY_GOLAY:
        if HAS_SCIPY and config.window_length < n and config.window_length % 2 == 1:
            try:
                return savgol_filter(
                    y, 
                    window_length=config.window_length, 
                    polyorder=min(config.polyorder, config.window_length - 2), 
                    mode="interp"
                )
            except Exception:
                # Fallback to moving average if Savitzky-Golay fails
                pass
        # Fallback to moving average if SciPy unavailable or window invalid
        window = min(config.window_length, max(3, n//10))
        if window % 2 == 0:
            window += 1
        return pd.Series(y, dtype=float).rolling(window=window, center=True, min_periods=1).mean().values
    
    elif config.method == SmoothingConfig.MOVING_AVERAGE:
        window = min(config.window_length, max(3, n//10))
        return pd.Series(y, dtype=float).rolling(window=window, center=True, min_periods=1).mean().values
    
    else:
        # Unknown method, return original data
        return y


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
    Compute per-segment capacity (mAh) for constant-sign current segments.

    Preference order for capacity source:
      1) Use instrument-measured charge columns if available
         - discharge: 'discharge_capacity_mAh_raw'
         - charge:    'charge_capacity_mAh_raw'
         - fallback:  'capacity_mAh_raw'
      2) Otherwise, numerically integrate current over time: Q = ∫ I dt

    The returned 'cap_mAh_seg' is zeroed at the start of each segment.
    """
    out = df.copy()
    out["cap_mAh_seg"] = np.nan

    has_qchg = "charge_capacity_mAh_raw" in out.columns and out["charge_capacity_mAh_raw"].notna().any()
    has_qdch = "discharge_capacity_mAh_raw" in out.columns and out["discharge_capacity_mAh_raw"].notna().any()
    has_qtot = "capacity_mAh_raw" in out.columns and out["capacity_mAh_raw"].notna().any()

    for sid, g in out.groupby("seg_id", sort=True):
        if len(g) == 0:
            continue
        is_dis = bool(g["is_discharge"].iloc[0])

        series = None
        if is_dis and has_qdch:
            series = g["discharge_capacity_mAh_raw"].astype(float).values
        elif (not is_dis) and has_qchg:
            series = g["charge_capacity_mAh_raw"].astype(float).values
        elif has_qtot:
            series = g["capacity_mAh_raw"].astype(float).values

        if series is not None and np.isfinite(series).any():
            # Use measured charge; reset to segment start
            dq = series - (series[0] if len(series) else 0.0)
        else:
            # Fallback: integrate current
            dq = np.cumsum(g["current_a"].values * g["dt_s"].values) / 3600.0 * 1000.0
            dq = dq - (dq[0] if len(dq) else 0.0)

        out.loc[g.index, "cap_mAh_seg"] = dq

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
        # Use nan-aware finite max to avoid RuntimeWarning when all values are NaN
        cap_end_mAh=(
            "cap_mAh_seg",
            lambda x: (
                float(np.nanmax(np.abs(x)))
                if np.isfinite(x.to_numpy(dtype=float)).any()
                else float("nan")
            ),
        ),
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
# Plotting
# -----------------------------


def plot_potential_vs_capacity(sample_id: str, df: pd.DataFrame, mass_g: float, outpath: Path) -> None:
    """
    Overlay charge/discharge segments across cycles in one plot: E vs Q (mAh/g).
    Colors represent C-rates (0.1C, 1C, 2C), with vertical lines for theoretical capacities.
    """
    fig, ax = plt.subplots(figsize=(8.5, 4.5))  # Wider to accommodate external legends
    # Professional ticks and grid
    ax.tick_params(axis="both", which="major", direction="out", length=5, width=0.8)
    ax.tick_params(axis="both", which="minor", direction="out", length=3, width=0.6)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(True, which="major", linestyle="-", linewidth=0.6, alpha=0.25)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.15)
    
    # Define colors for each C-rate
    c_rate_colors = {
        0.1: '#1f77b4',  # blue
        1.0: '#ff7f0e',  # orange  
        2.0: '#d62728'   # red
    }
    
    # Tolerance for C-rate matching (±20%)
    c_rate_tolerance = 0.2

    # Ensure normalized capacity per segment exists
    df = df.copy()
    # Convert integrated capacity (mAh) to specific capacity (mAh/g)
    df["cap_mAh_g_seg"] = df["cap_mAh_seg"] / mass_g

    # Track which C-rates are present for legend
    c_rates_present = set()
    labelled_crates = set()

    # Plot segments grouped by C-rate
    for cyc, gcyc in df.groupby("cycle_est", sort=True):
        if cyc == 0:
            continue
        # plot each segment in this cycle
        for sid, gseg in gcyc.groupby("seg_id", sort=True):
            if gseg["cap_mAh_g_seg"].abs().max() < 1e-6:
                continue
            
            # Calculate C-rate for this segment
            mean_current = gseg["current_a"].abs().mean()
            c_rate = estimate_c_rate(mean_current, mass_g, THEORETICAL_CAPACITY_MAH_G)
            
            # Find closest C-rate category
            closest_crate = None
            for target_crate in [0.1, 1.0, 2.0]:
                if abs(c_rate - target_crate) <= target_crate * c_rate_tolerance:
                    closest_crate = target_crate
                    break
            
            if closest_crate is None:
                continue  # Skip segments that don't match any C-rate category
                
            c_rates_present.add(closest_crate)
            color = c_rate_colors[closest_crate]
            
            is_discharge = bool(gseg["is_discharge"].iloc[0])
            
            # Label only once per C-rate (use discharge as representative)
            if is_discharge and closest_crate not in labelled_crates:
                lbl = f"{closest_crate}C"
                labelled_crates.add(closest_crate)
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

    ax.set_xlabel("Capacity [mAh/g]", fontsize=14)
    ax.set_ylabel("Potential [V]", fontsize=14)
    # Nice x-limits - ensure Graphite line at 372 mAh/g is visible
    try:
        qmax = float(np.nanmax(np.abs(df["cap_mAh_g_seg"].values)))
        ax.set_xlim(0, max(400.0, qmax * 1.05))  # Always show at least up to 400 mAh/g
    except Exception:
        ax.set_xlim(0, 400.0)  # Default range if data processing fails

    # Theoretical capacity markers (annotated lines, not part of legend)
    try:
        y0, y1 = ax.get_ylim()
        x_nmc = 200
        x_graphite = 372
        ax.axvline(x=x_nmc, color="black", linestyle=":", alpha=0.8, linewidth=1.5)
        ax.axvline(x=x_graphite, color="purple", linestyle=":", alpha=0.8, linewidth=1.5)
        ax.annotate(
            "NMC-811 (~200 mAh/g)",
            xy=(x_nmc, y1),
            xytext=(8, -8),
            textcoords="offset points",
            rotation=90,
            va="top",
            ha="left",
            color="black",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"),
        )
        ax.annotate(
            "Graphite (372 mAh/g)",
            xy=(x_graphite, y1),
            xytext=(5, -5),
            textcoords="offset points",
            rotation=90,
            va="top",
            ha="left",
            color="purple",
            fontsize=8,
        )
    except Exception:
        pass

    # Legends outside right: one for C-rates, one for style
    if c_rates_present:
        cr_list = sorted(list(c_rates_present))
        cr_handles = [
            Line2D([0], [0], color=c_rate_colors[c], lw=1.8, ls="-") for c in cr_list
        ]
        cr_labels = []
        for c in cr_list:
            if c == 0.1:
                cr_labels.append("C/10")
            else:
                cr_labels.append(f"{c:g}C")
        leg1 = ax.legend(
            cr_handles,
            cr_labels,
            title="C-Rates",
            fontsize=8,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            framealpha=0.95,
        )
        ax.add_artist(leg1)

    style_handles = [
        Line2D([0], [0], color="0.2", lw=1.6, ls="-", label="Discharge"),
        Line2D([0], [0], color="0.2", lw=1.2, ls="--", label="Charge"),
    ]
    ax.legend(
        style_handles,
        [h.get_label() for h in style_handles],
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.62),
        borderaxespad=0.0,
        framealpha=0.95,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_dqdE_cathodic(sample_id: str, df: pd.DataFrame, mass_g: float, outpath: Path, smoothing_config: SmoothingConfig = None) -> None:
    """
    Compute and plot dQ/dE vs E for cathodic (discharge) segments.
    
    Args:
        sample_id: Sample identifier
        df: DataFrame with processed data
        mass_g: Mass in grams
        outpath: Output path for the plot
        smoothing_config: SmoothingConfig object for dQdE smoothing
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
            # Use configurable smoothing for E and Q
            if smoothing_config is None:
                smoothing_config = DEFAULT_SMOOTHING
            Es = smooth(E, smoothing_config)
            Qs = smooth(Q, smoothing_config)
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
            # Use a lighter smoothing for the derivative (smaller window)
            derivative_config = SmoothingConfig(
                method=smoothing_config.method,
                window_length=min(101, max(7, (len(dQ_dE)//50)*2 + 1)),
                polyorder=smoothing_config.polyorder
            )
            dQ_dE = smooth(dQ_dE, derivative_config)
            ax.plot(Es, dQ_dE, lw=1.4, color=colors[color_idx % 10], label=f"Cycle {cyc}")
            color_idx += 1

    ax.set_xlabel("E [V]", fontsize=14)
    ax.set_ylabel("dQ/dE [mAh/(g·V)]", fontsize=14)
    # No legend for combined plot
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_dqdE_cathodic_by_crate(
    sample_id: str,
    df: pd.DataFrame,
    mass_g: float,
    target_c_rate: float,
    outpath: Path,
    smoothing_config: SmoothingConfig = None,
    target_currents_mA: Optional[List[float]] = None,
    current_tolerance: float = 0.15,
) -> None:
    """
    Compute and plot dQ/dE vs E for cathodic (discharge) segments filtered by
    either C-rate or explicit current levels.

    Args:
        sample_id: Sample identifier for plot title
        df: DataFrame with processed data
        mass_g: Mass in grams
        target_c_rate: Target C-rate to filter by (e.g., 0.1, 1.0, 2.0)
        outpath: Output path for the plot
        smoothing_config: SmoothingConfig object for dQdE smoothing
        target_currents_mA: If provided, match segments by mean |I| near any of these
            absolute current levels (in mA). Overrides C-rate filtering.
        current_tolerance: Fractional tolerance for matching current levels (default 0.15 → ±15%).
    """
    fig, ax = plt.subplots(figsize=(8, 4))  # Wider to accommodate external legend
    colors = plt.cm.tab10.colors
    color_idx = 0
    
    # Tolerance for C-rate matching (±20%)
    c_rate_tolerance = 0.2

    # Normalized capacity within segment
    df = df.copy()
    # Specific capacity in mAh/g for discharge segments
    df["cap_mAh_g_seg"] = df["cap_mAh_seg"] / mass_g

    # Calculate selection for each segment based on either C-rate or target currents
    segments_selected = []
    for cyc, gcyc in df.groupby("cycle_est", sort=True):
        if cyc == 0:
            continue
        for sid, gseg in gcyc.groupby("seg_id", sort=True):
            if not bool(gseg["is_discharge"].iloc[0]):
                continue
            
            # Mean current and C-rate for this segment
            mean_current_A = gseg["current_a"].abs().mean()
            mean_current_mA = mean_current_A * 1000.0
            c_rate = estimate_c_rate(mean_current_A, mass_g, THEORETICAL_CAPACITY_MAH_G)

            # Check selection criterion
            if target_currents_mA and len(target_currents_mA) > 0:
                for I_target in target_currents_mA:
                    if I_target is None or not np.isfinite(I_target):
                        continue
                    I_target_abs = abs(float(I_target))
                    if I_target_abs == 0:
                        continue
                    if abs(mean_current_mA - I_target_abs) <= I_target_abs * current_tolerance:
                        segments_selected.append((cyc, sid, gseg, c_rate))
                        break
            else:
                if abs(c_rate - target_c_rate) <= max(target_c_rate * c_rate_tolerance, 1e-9):
                    segments_selected.append((cyc, sid, gseg, c_rate))

    if not segments_selected:
        # No segments found for the requested selection
        if target_currents_mA and len(target_currents_mA) > 0:
            pretty_I = ", ".join(f"{i:.3g} mA" for i in target_currents_mA if np.isfinite(i))
            msg = f"No discharge segments found for I≈[{pretty_I}] (±{current_tolerance*100:.0f}%)"
        else:
            msg = f"No discharge segments found for {target_c_rate}C (±{c_rate_tolerance*100:.0f}%)"
        ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_xlabel("E [V]", fontsize=14)
        ax.set_ylabel("dQ/dE [mAh/(g·V)]", fontsize=14)
        fig.tight_layout()
        fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return

    # Plot segments that match the target C-rate
    for cyc, sid, gseg, c_rate in segments_selected:
        E = gseg["voltage_v"].values
        Q = gseg["cap_mAh_g_seg"].abs().values
        if len(E) < 5:
            continue
            
        # Use configurable smoothing for E and Q
        if smoothing_config is None:
            smoothing_config = DEFAULT_SMOOTHING
        Es = smooth(E, smoothing_config)
        Qs = smooth(Q, smoothing_config)
            
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
        # Use a lighter smoothing for the derivative (smaller window)
        derivative_config = SmoothingConfig(
            method=smoothing_config.method,
            window_length=min(101, max(7, (len(dQ_dE)//50)*2 + 1)),
            polyorder=smoothing_config.polyorder
        )
        dQ_dE = smooth(dQ_dE, derivative_config)
        ax.plot(Es, dQ_dE, lw=1.4, color=colors[color_idx % 10], 
                label=f"Cycle {cyc}")
        color_idx += 1

    ax.set_xlabel("E [V]", fontsize=14)
    ax.set_ylabel("dQ/dE [mAh/(g·V)]", fontsize=14)
    if color_idx > 0:
        ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# New: simple voltage vs time diagnostic plot
def plot_voltage_over_time(sample_id: str, df: pd.DataFrame, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    t_hours = df["time_s"].values / 3600.0
    e = df["voltage_v"].values

    # Grid and minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(True, which="major", linestyle="-", linewidth=0.6, alpha=0.25)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.15)

    # Single trace: Voltage vs Time
    ax.plot(t_hours, e, lw=1.2, color="tab:blue")
    ax.set_xlabel(r"Time [h]", fontsize=14, color='black')
    ax.set_ylabel(r"Potential [V]", fontsize=14, color='black')
    ax.tick_params(axis='both', colors='black', labelcolor='black')

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# -----------------------------
# Runner
# -----------------------------
def parse_dqde_cli(argv: Optional[List[str]]) -> Tuple[Optional[List[float]], float]:
    """
    Parse optional CLI args for dQ/dE current-based selection.
    --dqde-currents=2.38,0.24,4.75   (mA)
    --dqde-tol=0.15                  (fraction, ±15%)
    """
    currents: Optional[List[float]] = None
    tol: float = 0.15
    if not argv:
        return currents, tol
    for arg in argv:
        if arg.startswith("--dqde-currents="):
            try:
                val = arg.split("=", 1)[1]
                toks = re.split(r"[;,\s]+", val.strip())
                parsed: List[float] = []
                for t in toks:
                    if not t:
                        continue
                    try:
                        parsed.append(float(t))
                    except Exception:
                        pass
                if parsed:
                    currents = parsed
            except Exception:
                pass
        elif arg.startswith("--dqde-tol="):
            try:
                tol = float(arg.split("=", 1)[1])
            except Exception:
                pass
    return currents, tol


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

    # Voltage vs Time plot
    vt_path = FIG_DIR / f"T3.1_V_vs_t{PLOT_EXT}"
    plot_voltage_over_time(sample_id, df, vt_path)

    # Potential vs Capacity overlay
    pvc_path = FIG_DIR / f"T3.1_Potential_vs_Capacity_{PLOT_EXT}"
    try:
        plot_potential_vs_capacity(sample_id, df, mass_g, pvc_path)
    except Exception as e:
        print(f"[WARN] Potential vs Capacity plot failed: {e}")

    # dQ/dE overlay for all discharge segments
    try:
        dqdE_path = FIG_DIR / f"T3.1_dQdE_{PLOT_EXT}"
        plot_dqdE_cathodic(sample_id, df, mass_g, dqdE_path)
    except Exception as e:
        print(f"[WARN] dQ/dE overlay plot failed: {e}")

    # dQ/dE per target rate, optionally matched by explicit current levels via CLI
    try:
        currents_mA, tol = parse_dqde_cli(sys.argv)
    except Exception:
        currents_mA, tol = None, 0.15

    target_crates = [0.1, 1.0, 2.0]
    for idx, target_c_rate in enumerate(target_crates):
        outpath = FIG_DIR / f"T3.1_dQdE_{target_c_rate:.1f}C_{PLOT_EXT}"
        extra_kwargs = {}
        if currents_mA and idx < len(currents_mA) and np.isfinite(currents_mA[idx]):
            extra_kwargs["target_currents_mA"] = [float(currents_mA[idx])]
            extra_kwargs["current_tolerance"] = float(tol)
        try:
            plot_dqdE_cathodic_by_crate(
                sample_id,
                df,
                mass_g,
                target_c_rate,
                outpath,
                smoothing_config=None,
                **extra_kwargs,
            )
        except Exception as e:
            print(f"[WARN] dQ/dE {target_c_rate:.1f}C plot failed: {e}")


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