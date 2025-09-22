import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================
# Configuration and constants
# =============================================================
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"

# Glassy carbon working electrode (Tasks 3.2/3.3)
GC_DIAMETER_MM = 5.5
GC_RADIUS_CM = (GC_DIAMETER_MM / 10.0) / 2.0
GC_AREA_CM2 = math.pi * (GC_RADIUS_CM ** 2)

# Electrolyte / chemistry (Tasks 3.2 & 3.3)
N_ELECTRONS = 1  # Ferri/ferrocyanide is 1e-
T_K = 298.15
R_GAS = 8.314462618  # J/mol/K
F_CONST = 96485.33212  # C/mol

# 2.3 mM ferricyanide, used in Tasks 3.2 and 3.3
C_BULK_MOL_PER_L = 2.3e-3
C_BULK_MOL_PER_CM3 = C_BULK_MOL_PER_L / 1000.0  # mol/cm^3
C_BULK_mM = 2.3  # for Sand's eq. numeric form in the script (mM)

# For rotating-disk hydrodynamics (Task 3.3):
# Kinematic viscosity ν for 2.0 M KCl @ 25 °C, derived from NIST tables (Kestin et al., 1981).
# Steps:
#   1) Convert molarity (M, mol/L) → molality (m, mol/kg) using an estimated solution density ρ.
#        m = 1000*M / (1000*ρ - M*M_r), with ρ in g/mL and M_r in g/mol
#   2) Interpolate NIST ν between 2.0 and 2.5 mol/kg at 25 °C:
#        ν(2.0 m) = 0.8425 mm²/s,  ν(2.5 m) = 0.8324 mm²/s
#   3) Convert mm²/s → cm²/s (1 mm²/s = 0.01 cm²/s)
#
# Note: The exact ρ for 2.0 M KCl @ 25 °C varies slightly by source (~1.10–1.12 g/mL).
#       Using ρ = 1.10 g/mL gives m ≈ 2.103 mol/kg and ν ≈ 0.008404 cm²/s.

KCL_MOLAR_MASS_GPMOL = 74.5513      # g/mol
KCL_DENSITY_2M_25C_G_PER_ML = 1.10  # g/mL
M_KCL_MOLARITY = 2.0                 # mol/L

# 1) M -> m
m_KCL_molal = (1000.0 * M_KCL_MOLARITY) / (1000.0 * KCL_DENSITY_2M_25C_G_PER_ML - M_KCL_MOLARITY * KCL_MOLAR_MASS_GPMOL)

# 2) Interpolate ν(m) between the two nearest NIST tabulated points at 25 °C
m1, nu1_mm2s = 2.0, 0.8425  # 2.0 mol/kg
m2, nu2_mm2s = 2.5, 0.8324  # 2.5 mol/kg
nu_mm2s_interp = nu1_mm2s + (m_KCL_molal - m1) * (nu2_mm2s - nu1_mm2s) / (m2 - m1)

# 3) Convert to cm²/s
NU_CMS2 = nu_mm2s_interp * 0.01  # <-- kinematic viscosity to use downstream

# (Optional) One-line summary for your logs:
print(f"[ν] 2.0 M KCl @25°C: ρ={KCL_DENSITY_2M_25C_G_PER_ML:.3f} g/mL → m={m_KCL_molal:.3f} mol/kg → "
      f"ν≈{nu_mm2s_interp:.4f} mm²/s → NU_CMS2={NU_CMS2:.6f} cm²/s")


# Potentials vs Ag/AgCl unless otherwise specified

# =============================================================
# Helpers
# =============================================================
def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)


def setup_plot_style() -> None:
    """Set a consistent, clean plotting style across all figures."""
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.3,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "legend.frameon": True,
        "legend.fontsize": 9,
        "legend.framealpha": 0.85,
        "legend.edgecolor": "#00000055",
        "legend.title_fontsize": 9,
        "savefig.dpi": 300,
    })


def beautify_axes(ax: plt.Axes) -> None:
    """Apply minor ticks, outward ticks, and simplify spines on an Axes."""
    try:
        ax.minorticks_on()
    except Exception:
        pass
    ax.tick_params(which="both", direction="out", length=5, width=0.8)
    ax.tick_params(which="minor", length=3, width=0.6)
    # Hide top/right spines and standardize spine width
    if "top" in ax.spines:
        ax.spines["top"].set_visible(False)
    if "right" in ax.spines:
        ax.spines["right"].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def _detect_header_start(path: Path) -> int:
    """Return the 0-based line index where the tabular header starts.

    Heuristics for BioLogic .txt: header begins at a line that starts with
    one of: 'freq/Hz', 'time/s', or 'mode'.
    """
    header_markers = ("freq/Hz", "time/s", "mode")
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for idx, line in enumerate(f):
            first_token = line.strip().split("\t", 1)[0]
            if first_token in header_markers:
                return idx
    return 0


def read_biologic_table(path: Path) -> pd.DataFrame:
    header_idx = _detect_header_start(path)
    # Some BioLogic exports on Windows use cp1252/latin-1 (e.g., 'µ' in headers).
    # Try UTF-8 first, then fall back to latin-1 for robustness.
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            engine="python",
            skiprows=header_idx,
            comment="#",
            encoding="utf-8",
        )
    except UnicodeDecodeError:
        df = pd.read_csv(
            path,
            sep="\t",
            engine="python",
            skiprows=header_idx,
            comment="#",
            encoding="latin-1",
        )
    df.columns = [str(c).strip() for c in df.columns]
    return df


def get_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def label_from_filename(path: Path) -> str:
    return path.stem


def extract_rpm_from_name(path: Path) -> Optional[int]:
    m = re.search(r"(\d+)\s*rpm", path.name, flags=re.IGNORECASE)
    if not m:
        m = re.search(r"(\d+)rpm", path.name, flags=re.IGNORECASE)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def safe_save(fig: plt.Figure, filename: str) -> None:
    ensure_dirs()
    out = FIG_DIR / filename
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out.relative_to(ROOT)}")


def write_csv(df: pd.DataFrame, filename: str) -> None:
    ensure_dirs()
    out = RES_DIR / filename
    df.to_csv(out, index=False)
    print(f"Saved: {out.relative_to(ROOT)}")


# =============================================================
# Global diffusion-coefficient stash (in-memory for this run)
# =============================================================
GLOBAL_D_VALUES: List[float] = []

def _add_D(value: float, source: str = "") -> None:
    """Store a positive, finite diffusion coefficient [cm^2/s] computed in this run."""
    try:
        v = float(value)
    except Exception:
        return
    if np.isfinite(v) and v > 0:
        GLOBAL_D_VALUES.append(v)

def get_calculated_D() -> float:
    """Return the median D [cm^2/s] accumulated so far, or NaN if none."""
    return float(np.nanmedian(GLOBAL_D_VALUES)) if GLOBAL_D_VALUES else float("nan")


# =============================================================
# CV analysis (Task 3.2.1)
# =============================================================
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
    """Estimate scan rate |dE/dt| in V/s using robust median of derivative.
    Returns None if time is missing.
    """
    if t_s is None:
        return None
    e = np.asarray(e_v)
    t = np.asarray(t_s)
    if len(e) < 5 or len(t) < 5:
        return None
    # numerical derivative with guard against zero/neg dt
    dt = np.gradient(t)
    dt[dt == 0] = np.nan
    dEdt = np.gradient(e) / dt
    v = np.nanmedian(np.abs(dEdt))
    return float(v) if np.isfinite(v) else None


def _pick_last_cycle(df: pd.DataFrame) -> pd.DataFrame:
    cyc_col = get_column(df, ["cycle number", "cycle", "Cycle"]) 
    if cyc_col is None:
        return df
    cyc = pd.to_numeric(df[cyc_col], errors="coerce")
    if cyc.notna().any():
        last = int(cyc.max())
        return df[cyc == last]
    return df


def _find_cv_extrema(e_v: np.ndarray, i_A: np.ndarray) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Return (Epa, ipa) and (Epc, ipc) using simple extremum search on smoothed current.
    """
    # Light smoothing (moving average) to reduce noise without distorting peaks
    k = max(5, len(i_A) // 200)
    if k % 2 == 0:
        k += 1
    if k > 51:
        k = 51
    s = pd.Series(i_A).rolling(k, center=True, min_periods=1).mean().to_numpy()

    ipa = float(np.nanmax(s))
    Epa = float(e_v[np.nanargmax(s)])
    ipc = float(np.nanmin(s))
    Epc = float(e_v[np.nanargmin(s)])
    return (Epa, ipa), (Epc, ipc)


def analyze_task32_cvs() -> None:
    cv_dir = DATA_DIR / "Task 3.2 CV"
    files = sorted(cv_dir.glob("*_CV_C01.txt"))
    if not files:
        print("[CV] No CV files found; skipping.")
        return

    rows: List[CVPeak] = []

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for path in files:
        try:
            df = read_biologic_table(path)
        except Exception as e:
            print(f"[CV] Failed to read {path.name}: {e}")
            continue

        df_use = _pick_last_cycle(df)
        e_col = get_column(df_use, ["Ewe/V"]) or "Ewe/V"
        i_col = get_column(df_use, ["<I>/mA", "I/mA"]) or "<I>/mA"
        t_col = get_column(df_use, ["time/s", "Time/s"])
        if e_col not in df_use or i_col not in df_use:
            print(f"[CV] Missing columns in {path.name}; skipping.")
            continue

        e_v = pd.to_numeric(df_use[e_col], errors="coerce").to_numpy()
        i_mA = pd.to_numeric(df_use[i_col], errors="coerce").to_numpy()
        i_A = i_mA / 1000.0
        j_mA_cm2 = i_mA / GC_AREA_CM2

        # Scan rate estimate for legend (if time available)
        scan_rate = None
        if t_col:
            t_s = pd.to_numeric(df_use[t_col], errors="coerce").to_numpy()
            scan_rate = estimate_scan_rate_Vs(e_v, t_s)
        legend_label = (f"v≈{(scan_rate*1000):.0f} mV/s" if (scan_rate is not None and np.isfinite(scan_rate))
                        else label_from_filename(path))

        # Plot overlay (area-normalized current)
        ax.plot(e_v, j_mA_cm2, lw=1.1, label=legend_label)

        # Peak analysis
        (Epa, ipa_A), (Epc, ipc_A) = _find_cv_extrema(e_v, i_A)
        dEp = Epa - Epc
        ratio = abs(ipa_A) / abs(ipc_A) if ipc_A != 0 else np.nan

        rows.append(CVPeak(
            scan_rate_Vs=scan_rate if scan_rate is not None else float("nan"),
            ipa_A=float(ipa_A), Epa_V=float(Epa),
            ipc_A=float(ipc_A), Epc_V=float(Epc),
            delta_Ep_V=float(dEp), ipa_ipc_ratio=float(ratio)
        ))

    ax.set_xlabel("E (V vs Ag/AgCl)")
    ax.set_ylabel("j (mA cm$^{-2}$)")
    ax.set_title("Task 3.2 – CVs (last cycle per file)")
    ax.legend(title="Scan rate (approx.)", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    safe_save(fig, "T3.2_CV_overlay.png")

    # Summary table + Randles–Sevcik diffusion coefficient from |i_p| vs sqrt(v)
    if rows:
        df_sum = pd.DataFrame([r.__dict__ for r in rows])
        write_csv(df_sum, "T3.2_CV_peaks.csv")

        # Randles–Sevcik: i_p = 2.69e5 * n^(3/2) * A * C * D^(1/2) * v^(1/2)  (A, cm, mol/cm^3, V/s)
        df_rs = df_sum.dropna(subset=["scan_rate_Vs"])  # need v
        if len(df_rs) >= 2:
            v = np.sqrt(df_rs["scan_rate_Vs"].to_numpy())
            ip = np.abs(df_rs[["ipa_A", "ipc_A"]]).max(axis=1).to_numpy()  # use larger magnitude
            # Linear fit ip vs sqrt(v)
            A_mat = np.vstack([v, np.ones_like(v)]).T
            slope, intercept = np.linalg.lstsq(A_mat, ip, rcond=None)[0]
            const = 2.69e5 * (N_ELECTRONS ** 1.5) * GC_AREA_CM2 * C_BULK_MOL_PER_CM3
            D_half = slope / const
            D_RS = (D_half ** 2)  # cm^2/s
            rs_out = pd.DataFrame({
                "slope_A_per_Vs_sqrt": [slope],
                "intercept_A": [intercept],
                "D_RandlesSevcik_cm2_s": [D_RS],
            })
            write_csv(rs_out, "T3.2_CV_randles_sevcik.csv")
            _add_D(D_RS, "Randles-Sevcik (combined)")

            # ---------- NEW: Produce the required Randles–Sevcik plot (both peaks) ----------
            v_sqrt = v  # alias
            ipa_A = np.abs(df_rs["ipa_A"].to_numpy())
            ipc_A = np.abs(df_rs["ipc_A"].to_numpy())

            # Fits for anodic and cathodic peaks separately (in A)
            A_fit = np.vstack([v_sqrt, np.ones_like(v_sqrt)]).T
            slope_a, intercept_a = np.linalg.lstsq(A_fit, ipa_A, rcond=None)[0]
            slope_c, intercept_c = np.linalg.lstsq(A_fit, ipc_A, rcond=None)[0]

            # Diffusion coefficients from each slope
            D_half_a = slope_a / const
            D_half_c = slope_c / const
            D_RS_a = float(D_half_a ** 2) if np.isfinite(D_half_a) else np.nan
            D_RS_c = float(D_half_c ** 2) if np.isfinite(D_half_c) else np.nan
            _add_D(D_RS_a, "Randles-Sevcik (anodic)")
            _add_D(D_RS_c, "Randles-Sevcik (cathodic)")

            # Save a compact CSV with both lines (non-breaking with the original file)
            write_csv(pd.DataFrame({
                "slope_anodic_A_per_Vs_sqrt": [float(slope_a)],
                "intercept_anodic_A": [float(intercept_a)],
                "D_RandlesSevcik_anodic_cm2_s": [D_RS_a],
                "slope_cathodic_A_per_Vs_sqrt": [float(slope_c)],
                "intercept_cathodic_A": [float(intercept_c)],
                "D_RandlesSevcik_cathodic_cm2_s": [D_RS_c],
            }), "T3.2_CV_randles_sevcik_fits.csv")

            # Plot (use mA for readability but compute fits in A)
            fig_rs, ax_rs = plt.subplots(figsize=(6.0, 4.2))
            ax_rs.plot(v_sqrt, ipa_A * 1e3, "o", label="|i_pa| (anodic)")
            ax_rs.plot(v_sqrt, ipc_A * 1e3, "s", label="|i_pc| (cathodic)")

            xfit = np.linspace(np.nanmin(v_sqrt), np.nanmax(v_sqrt), 100)
            yfit_a_mA = (slope_a * xfit + intercept_a) * 1e3
            yfit_c_mA = (slope_c * xfit + intercept_c) * 1e3
            ax_rs.plot(xfit, yfit_a_mA, "-", label=f"fit anodic; D≈{D_RS_a:.2e} cm²/s")
            ax_rs.plot(xfit, yfit_c_mA, "--", label=f"fit cathodic; D≈{D_RS_c:.2e} cm²/s")

            ax_rs.set_xlabel(r"$\sqrt{v}$ (V$^{1/2}$ s$^{-1/2}$)")
            ax_rs.set_ylabel(r"$i_p$ (mA)")
            ax_rs.set_title("Task 3.2 – Randles–Sevcik (peak current vs $\\sqrt{v}$)")
            ax_rs.legend(title="Peaks and fits", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
            beautify_axes(ax_rs)
            safe_save(fig_rs, "T3.2_CV_RandlesSevcik.png")
        else:
            print("[CV] Need >=2 scan rates with time to estimate D via Randles–Sevcik.")


# =============================================================
# Chronopotentiometry (Task 3.2.3)
# =============================================================
@dataclass
class CPResult:
    label: str
    current_mA: float
    tau_s: float
    E_tau_over_4_V: float
    D_Sand_cm2_s: float


def _estimate_tau_s(t_s: np.ndarray, e_v: np.ndarray) -> float:
    """Estimate transition time τ as the point of maximum curvature of E(t).
    Fallback to 90% of range if needed.
    """
    # Smooth E slightly
    k = max(11, len(e_v) // 200)
    if k % 2 == 0:
        k += 1
    if k > 101:
        k = 101
    e_sm = pd.Series(e_v).rolling(k, center=True, min_periods=1).mean().to_numpy()
    de_dt = np.gradient(e_sm, t_s)
    d2e_dt2 = np.gradient(de_dt, t_s)

    # Focus on middle region to avoid startup artifacts
    lo, hi = int(0.05 * len(t_s)), int(0.95 * len(t_s))
    idx = np.argmax(np.abs(d2e_dt2[lo:hi])) + lo
    tau = float(t_s[idx])

    if not np.isfinite(tau) or tau <= 0:
        # Fallback: time when E reaches 90% of its dynamic range
        emin, emax = np.nanmin(e_sm), np.nanmax(e_sm)
        target = emin + 0.9 * (emax - emin)
        # find first crossing
        gt = np.where(e_sm >= target)[0]
        tau = float(t_s[gt[0]]) if len(gt) else float(t_s[-1])
    return tau


def analyze_task32_cp() -> None:
    cp_dir = DATA_DIR / "Task 3.2 CP"
    files = sorted(cp_dir.glob("*_CP_C01.txt"))
    if not files:
        print("[CP] No CP files found; skipping.")
        return

    results: List[CPResult] = []

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for path in files:
        try:
            df = read_biologic_table(path)
        except Exception as e:
            print(f"[CP] Failed to read {path.name}: {e}")
            continue
        t_col = get_column(df, ["time/s", "Time/s"]) or "time/s"
        e_col = get_column(df, ["Ewe/V"]) or "Ewe/V"
        i_col = get_column(df, ["<I>/mA", "I/mA"])  # optional but helps
        if t_col not in df or e_col not in df:
            print(f"[CP] Missing time/E columns in {path.name}; skipping.")
            continue
        t_s = pd.to_numeric(df[t_col], errors="coerce").to_numpy()
        e_v = pd.to_numeric(df[e_col], errors="coerce").to_numpy()
        tau = _estimate_tau_s(t_s, e_v)
        # Quarter-wave potential by interpolation at t = tau/4
        t_target = tau / 4.0
        E_tau4 = float(np.interp(t_target, t_s, e_v))

        # Current magnitude (mA): average if available
        if i_col and i_col in df:
            i_mA = float(np.nanmean(pd.to_numeric(df[i_col], errors="coerce")))
        else:
            i_mA = np.nan

        # Legend label using applied current, prefer µA when small
        file_label = label_from_filename(path)
        if np.isfinite(i_mA):
            if abs(i_mA) < 1.0:
                legend_label = f"I≈{i_mA*1000:.0f} µA"
            else:
                legend_label = f"I≈{i_mA:.2f} mA"
        else:
            legend_label = file_label

        ax.plot(t_s, e_v, lw=1.0, label=legend_label)

        # Sand's equation in the numeric form given in skript.md Eq. [1.1.5]
        # i[mA] * sqrt(tau[s]) / C*[mM] = 85.5 * n * A[cm^2] * sqrt(D[cm^2/s])
        # -> sqrt(D) = (i * sqrt(tau) / (85.5 * n * A * C*))
        if np.isfinite(i_mA):
            D_half = (abs(i_mA) * math.sqrt(tau)) / (85.5 * N_ELECTRONS * GC_AREA_CM2 * C_BULK_mM)
            D_sand = D_half ** 2
        else:
            D_sand = np.nan
        _add_D(D_sand, "Sand")

        results.append(CPResult(label=file_label, current_mA=i_mA, tau_s=tau, E_tau_over_4_V=E_tau4, D_Sand_cm2_s=D_sand))

    ax.set_xlabel("t (s)")
    ax.set_ylabel("E (V vs Ag/AgCl)")
    ax.set_title("Task 3.2 – Chronopotentiometry (E vs t)")
    ax.legend(title="Applied current", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    safe_save(fig, "T3.2_CP_transients.png")

    if results:
        df_res = pd.DataFrame([r.__dict__ for r in results])
        write_csv(df_res, "T3.2_CP_summary.csv")


# =============================================================
# EIS (Task 3.2.2) – Nyquist + Warburg + parameter extraction
# =============================================================
@dataclass
class EISResult:
    label: str
    Rs_area_Ohm_cm2: float
    Rct_area_Ohm_cm2: float
    Cdl_F_cm2: float
    f_peak_Hz: float
    k0_cm_s: float
    Aw_Ohm_cm2_sqrt_s: float
    D_Warburg_cm2_s: float


def _extract_eis_parameters(df: pd.DataFrame) -> Tuple[float, float, float, float]:
    """Return (Rs_area, Rct_area, Cdl_per_cm2, f_peak). Uses simple semicircle heuristics.
    Assumes single semicircle (Randles) measured near formal potential.
    """
    re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
    im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"  # already -Im
    f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
    z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
    z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
    f_hz = pd.to_numeric(df[f_col], errors="coerce").to_numpy()

    # Electrolyte resistance (high-f intercept) ≈ min Z'
    Rs_area = float(np.nanmin(z_re))
    # Low-f intercept ≈ max Z'
    R_low_area = float(np.nanmax(z_re))
    Rct_area = max(R_low_area - Rs_area, 0.0)

    # f_peak at top of semicircle: index of max -Z''
    idx_peak = int(np.nanargmax(z_im))
    f_peak = float(f_hz[idx_peak]) if np.isfinite(f_hz[idx_peak]) and f_hz[idx_peak] > 0 else np.nan
    Cdl = float(1.0 / (2.0 * math.pi * f_peak * (Rct_area / GC_AREA_CM2))) if np.isfinite(f_peak) and Rct_area > 0 else np.nan
    # Cdl per cm^2:
    Cdl_per_cm2 = Cdl / GC_AREA_CM2
    return Rs_area, Rct_area, Cdl_per_cm2, f_peak


def _warburg_fit(df: pd.DataFrame) -> Tuple[float, float]:
    """Return (Aw, D_Warburg) using Z' vs 1/sqrt(ω) slope (low-frequency tail).
    Uses Eq. [3.3.2.1]: Aw = 4RT / (n^2 F^2 C* D^{1/2}).
    """
    re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
    im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
    f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
    z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
    z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
    f_hz = pd.to_numeric(df[f_col], errors="coerce").to_numpy()
    omega = 2 * np.pi * f_hz

    # Select "Warburg" region where |Z' - Z''| is small (≈45°) AND low frequencies
    mask = np.isfinite(z_re) & np.isfinite(z_im) & np.isfinite(omega) & (omega > 0)
    z_re, z_im, omega = z_re[mask], z_im[mask], omega[mask]
    angle45 = np.abs(z_re - z_im) / np.maximum(1e-9, (np.abs(z_re) + np.abs(z_im))) < 0.15
    # take lowest 30% frequencies as candidate + angle criterion
    lf_cut = np.quantile(omega, 0.3) if len(omega) else np.inf
    sel = (omega <= lf_cut) & angle45

    if sel.sum() < 5:
        return float("nan"), float("nan")

    x = 1.0 / np.sqrt(omega[sel])  # s^{1/2}
    y = z_re[sel]
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    Aw = float(slope)  # Ω·cm^2·s^{1/2}

    # D from Eq. [3.3.2.1]: Aw = 4RT / (n^2 F^2 C* sqrt(D))  => sqrt(D) = 4RT / (n^2 F^2 C* Aw)
    D_half = (4.0 * R_GAS * T_K) / ( (N_ELECTRONS**2) * (F_CONST**2) * C_BULK_MOL_PER_CM3 * Aw )
    D_w = float(D_half ** 2) if np.isfinite(D_half) and D_half > 0 else float("nan")
    return Aw, D_w


def analyze_task32_eis() -> None:
    eis_dir = DATA_DIR / "Task 3.2 EIS"
    files = sorted(eis_dir.glob("*_C01.txt"))
    if not files:
        print("[EIS] No EIS files found; skipping.")
        return

    # Nyquist overlay
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    results: List[EISResult] = []

    for path in files:
        try:
            df = read_biologic_table(path)
        except Exception as e:
            print(f"[EIS] Failed to read {path.name}: {e}")
            continue

        re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
        im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
        if re_col not in df or im_col not in df:
            print(f"[EIS] Missing impedance columns in {path.name}; skipping.")
            continue

        z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
        z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
        ax.plot(z_re, z_im, marker="o", ms=2.5, lw=0.8, label=label_from_filename(path))

        Rs_area, Rct_area, Cdl_per_cm2, f_peak = _extract_eis_parameters(df)
        # k0 from area-normalized Rct: k0 = RT / (n^2 F^2 C* Rct_area)
        k0 = (R_GAS * T_K) / ((N_ELECTRONS**2) * (F_CONST**2) * C_BULK_MOL_PER_CM3 * Rct_area) if Rct_area > 0 else np.nan
        Aw, D_w = _warburg_fit(df)

        _add_D(D_w, "Warburg")

        results.append(EISResult(
            label=label_from_filename(path),
            Rs_area_Ohm_cm2=Rs_area,
            Rct_area_Ohm_cm2=Rct_area,
            Cdl_F_cm2=Cdl_per_cm2,
            f_peak_Hz=f_peak,
            k0_cm_s=float(k0),
            Aw_Ohm_cm2_sqrt_s=Aw,
            D_Warburg_cm2_s=D_w,
        ))

    ax.set_xlabel(r"Z' (Ω·cm$^2$)")
    ax.set_ylabel(r"-Z'' (Ω·cm$^2$)")
    ax.set_title("Task 3.2 – EIS Nyquist (area-normalized)")
    ax.legend(title="Dataset", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    ax.set_aspect("equal", adjustable="box")
    safe_save(fig, "T3.2_EIS_Nyquist.png")

    # Warburg tail diagnostic & linearization plots
    figw, axw = plt.subplots(figsize=(6.4, 4.2))
    for path in files:
        try:
            df = read_biologic_table(path)
        except Exception:
            continue
        re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
        im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
        f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
        if not ({re_col, im_col, f_col} <= set(df.columns)):
            continue
        z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
        z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
        omega = 2 * np.pi * pd.to_numeric(df[f_col], errors="coerce").to_numpy()
        mask = np.isfinite(z_re) & np.isfinite(z_im) & np.isfinite(omega) & (omega > 0)
        x = 1.0 / np.sqrt(omega[mask])
        axw.plot(x, z_re[mask], lw=1.0, label=f"Re: {label_from_filename(path)}")
        axw.plot(x, z_im[mask], lw=1.0, ls="--", label=f"-Im: {label_from_filename(path)}")

    axw.set_xlabel(r"$1/\sqrt{\omega}$ (s$^{1/2}$)")
    axw.set_ylabel(r"Z (Ω·cm$^2$)")
    axw.set_title("Task 3.2 – Warburg linearization")
    axw.legend(title="Components", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8, ncol=1)
    beautify_axes(axw)
    safe_save(figw, "T3.2_EIS_Warburg_linearization.png")

    if results:
        df_eis = pd.DataFrame([r.__dict__ for r in results])
        write_csv(df_eis, "T3.2_EIS_summary.csv")


# =============================================================
# Task 3.3 – Rotating electrode: LSV (Koutecky–Levich) + Tafel
# =============================================================
@dataclass
class KLPoint:
    E_V: float
    inv_sqrt_omega_s05: float
    inv_j_cm2_A_inv: float


def _interpolate_to_common_grid(curves: Dict[int, pd.DataFrame]) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    """Interpolate j(E) for each rpm onto the common E-range grid (overlap).
    Returns (E_grid, j_by_rpm[A/cm^2]).
    """
    # Determine overlap in potential
    Emins, Emaxs = [], []
    for rpm, df in curves.items():
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
        current_candidates = [c for c in ["<I>/mA", "I/mA"] if c in df.columns]
        if not current_candidates:
            raise KeyError("No current column '<I>/mA' or 'I/mA' found in LSV data.")
        i_mA = pd.to_numeric(df[current_candidates[0]], errors="coerce").to_numpy()
        j_A_cm2 = (i_mA / 1000.0) / GC_AREA_CM2
        # Ensure monotonic E for interpolation by sorting by potential
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

    # ---------------- LSV overlay ----------------
    lsv_files = sorted([p for p in lsv_dir.glob("*_LSV_C01.txt") if "PEIS" not in p.name])
    curves: Dict[int, pd.DataFrame] = {}

    if lsv_files:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        for path in lsv_files:
            try:
                df = read_biologic_table(path)
            except Exception as e:
                print(f"[LSV] Failed to read {path.name}: {e}")
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
                # store for KL analysis
                df_use = df[[e_col, i_col]].copy()
                df_use.columns = ["Ewe/V", i_col]
                # ensure current column name matches expected accessor in interpolation helper
                if i_col != "<I>/mA":
                    df_use["<I>/mA"] = df_use[i_col]
                curves[rpm] = df_use

        ax.set_xlabel("E (V vs Ag/AgCl)")
        ax.set_ylabel("j (mA cm$^{-2}$)")
        ax.set_title("Task 3.3 – LSVs at various rotation rates")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        safe_save(fig, "T3.3_LSV_overlay.png")
    else:
        print("[LSV] No LSV files found; skipping LSV plot.")

    # ---------------- Koutecky–Levich ----------------
    if len(curves) >= 3:
        try:
            E_grid, j_by_rpm = _interpolate_to_common_grid(curves)
        except Exception as e:
            print(f"[KL] Could not interpolate to common grid: {e}")
            E_grid, j_by_rpm = None, None

        if E_grid is not None:
            rpms = np.array(sorted(j_by_rpm.keys()))
            omega = 2 * np.pi * (rpms / 60.0)  # rad/s
            inv_sqrt_omega = 1.0 / np.sqrt(omega)

            # Choose a small set of evenly spaced potentials in the grid to report
            E_targets = np.linspace(E_grid.min() + 0.05*np.ptp(E_grid), E_grid.max() - 0.05*np.ptp(E_grid), 8)
            kl_rows = []
            jk_vs_E = []

            for E0 in E_targets:
                j_at_E = np.array([j_by_rpm[rpm][np.argmin(np.abs(E_grid - E0))] for rpm in rpms])  # A/cm^2
                # avoid zeros: restrict to magnitudes > small
                mask = np.isfinite(j_at_E) & (np.abs(j_at_E) > 1e-6)
                if mask.sum() < 3:
                    continue
                y = 1.0 / j_at_E[mask]
                x = inv_sqrt_omega[mask]
                A = np.vstack([x, np.ones_like(x)]).T
                slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                jk = 1.0 / intercept  # kinetic current density at E0 (A/cm^2)
                kl_rows.extend([
                    {"E_V": float(E0), "inv_sqrt_omega_s05": float(xi), "inv_j_cm2_A_inv": float(yi)}
                    for xi, yi in zip(x, y)
                ])
                jk_vs_E.append((E0, jk))

            if jk_vs_E:
                df_kl = pd.DataFrame(kl_rows)
                write_csv(df_kl, "T3.3_KouteckyLevich_points.csv")

                df_jk = pd.DataFrame({"E_V": [p[0] for p in jk_vs_E], "j_k_A_cm2": [p[1] for p in jk_vs_E]}).sort_values("E_V")
                write_csv(df_jk, "T3.3_jk_vs_E.csv")

                # Plot KL at a representative potential (middle of list)
                mid = len(E_targets) // 2
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
                    axkl.plot(x, slope*x + intercept, "-")
                    axkl.set_xlabel(r"$1/\sqrt{\omega}$ (s$^{1/2}$)")
                    axkl.set_ylabel(r"$1/j$ (cm$^2$ A$^{-1}$)")
                    axkl.set_title(f"Koutecky–Levich at E = {E_mid:.3f} V")
                    beautify_axes(axkl)
                    safe_save(figkl, "T3.3_KL_example.png")

                # --------- Tafel from j_k(E) ---------
                # Estimate E_eq as the potential where j crosses zero at the slowest rotation (first rpm)
                # Fallback to median zero-cross across rpms
                zero_Es = []
                for rpm, j in j_by_rpm.items():
                    sgn = np.sign(j)
                    zcross = np.where(np.diff(sgn))[0]
                    if len(zcross):
                        idx = zcross[0]
                        E0 = float(np.interp(0.0, [j[idx], j[idx+1]], [E_grid[idx], E_grid[idx+1]]))
                        zero_Es.append(E0)
                E_eq = float(np.median(zero_Es)) if zero_Es else float(df_jk["E_V"].iloc[np.argmin(np.abs(df_jk["j_k_A_cm2"]))])

                # Build Tafel arrays; enforce single branch (anodic or cathodic) and adapt window if needed
                jk = df_jk.copy()
                jk["eta_V"] = jk["E_V"] - E_eq
                jk = jk[np.isfinite(jk["eta_V"]) & np.isfinite(jk["j_k_A_cm2"])].copy()
                jk["sign"] = np.sign(jk["j_k_A_cm2"])
                sign_to_use = 1.0 if (jk["sign"] > 0).sum() >= (jk["sign"] < 0).sum() else -1.0
                branch = jk[(jk["sign"] == sign_to_use) & (np.abs(jk["j_k_A_cm2"]) > 1e-8)].copy()

                if len(branch) >= 3:
                    branch["log10_jk"] = np.log10(np.abs(branch["j_k_A_cm2"]))
                    candidates = [(0.03, 0.15), (0.02, 0.20), (0.02, 0.25), (0.015, 0.30), (0.01, 0.35)]
                    best = None
                    for eta_min, eta_max in candidates:
                        sel = (np.abs(branch["eta_V"]) >= eta_min) & (np.abs(branch["eta_V"]) <= eta_max)
                        if sel.sum() >= 3:
                            x = branch.loc[sel, "log10_jk"].to_numpy()
                            y = branch.loc[sel, "eta_V"].to_numpy()
                            A = np.vstack([x, np.ones_like(x)]).T
                            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                            y_pred = slope * x + intercept
                            ss_res = np.sum((y - y_pred) ** 2)
                            ss_tot = np.sum((y - np.mean(y)) ** 2)
                            r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                            if (best is None) or (r2 > best["r2"]):
                                best = {"eta_min": eta_min, "eta_max": eta_max, "slope": float(slope), "intercept": float(intercept), "r2": float(r2), "x": x, "y": y}

                    # Fallback: widen to |eta| >= 20 mV if needed
                    if best is None:
                        sel = np.abs(branch["eta_V"]) >= 0.02
                        if sel.sum() >= 3:
                            x = branch.loc[sel, "log10_jk"].to_numpy()
                            y = branch.loc[sel, "eta_V"].to_numpy()
                            A = np.vstack([x, np.ones_like(x)]).T
                            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
                            y_pred = slope * x + intercept
                            ss_res = np.sum((y - y_pred) ** 2)
                            ss_tot = np.sum((y - np.mean(y)) ** 2)
                            r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
                            best = {"eta_min": 0.02, "eta_max": float(np.abs(branch["eta_V"]).max()), "slope": float(slope), "intercept": float(intercept), "r2": float(r2), "x": x, "y": y}

                    if best is not None:
                        # η = a * log10(j_k) + b  => at η=0, log10(j0) = -b/a
                        log10_j0 = -best["intercept"] / best["slope"]
                        j0 = 10 ** log10_j0

                        # Plot Tafel
                        figtf, axtf = plt.subplots(figsize=(5.4, 4.0))
                        axtf.plot(best["x"], best["y"], "o", label="data")
                        xfit = np.linspace(np.min(best["x"]), np.max(best["x"]), 100)
                        axtf.plot(xfit, best["slope"]*xfit + best["intercept"], "-", label=f"fit; R²={best['r2']:.2f}; j0≈{j0:.2e} A/cm²")
                        axtf.set_xlabel(r"log$_{10}$(|j$_k$| / A cm$^{-2}$)")
                        axtf.set_ylabel(r"$\eta$ (V)")
                        axtf.set_title("Task 3.3 – Tafel from kinetic currents")
                        beautify_axes(axtf)
                        axtf.legend(title="Data and fit", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
                        safe_save(figtf, "T3.3_Tafel.png")

                        write_csv(pd.DataFrame({
                            "slope_V_per_decade": [best["slope"]],
                            "intercept_V": [best["intercept"]],
                            "E_eq_est_V": [E_eq],
                            "j0_A_cm2": [float(j0)],
                            "R2": [best["r2"]],
                            "branch_sign": [int(sign_to_use)],
                            "eta_min_V": [best["eta_min"]],
                            "eta_max_V": [best["eta_max"]],
                        }), "T3.3_Tafel_fit.csv")
                    else:
                        print("[Tafel] Not enough points in any reasonable region to fit.")
                else:
                    print("[Tafel] Not enough points on a single branch to fit.")
    else:
        print("[KL] Need LSV data at ≥3 rotation rates for Koutecky–Levich and Tafel.")

    # ---------------- PEIS overlay + qualitative parameters per rpm ----------------
    peis_files = sorted([p for p in lsv_dir.glob("*_PEIS_C01.txt")])
    if peis_files:
        # Use in-memory D that has been accumulated during this run (RS, Sand, Warburg).
        D_calc = get_calculated_D()
        if not np.isfinite(D_calc):
            print("[PEIS] No calculated D available yet (RS/Sand/Warburg). δ will be NaN. Run Task 3.2 first in the same run.")
        fig, ax = plt.subplots(figsize=(6.0, 4.8))
        peis_rows = []
        for path in peis_files:
            try:
                df = read_biologic_table(path)
            except Exception as e:
                print(f"[PEIS] Failed to read {path.name}: {e}")
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

            # Qualitative parameter extraction (as in Task 3.2)
            Rs_area, Rct_area, Cdl_per_cm2, f_peak = _extract_eis_parameters(df)

            # Diffusion layer thickness δ ~ 1.61 * D^{1/3} * ν^{1/6} * ω^{-1/2}
            # Always use the calculated D from this run (no assumed fallback).
            if rpm is not None and rpm > 0:
                omega = 2 * np.pi * (rpm / 60.0)
                delta_cm = (1.61 * (D_calc ** (1/3)) * (NU_CMS2 ** (1/6)) * (omega ** (-0.5))) if np.isfinite(D_calc) else float("nan")
            else:
                delta_cm = np.nan

            peis_rows.append({
                "label": label,
                "rpm": rpm,
                "Rs_area_Ohm_cm2": Rs_area,
                "Rct_area_Ohm_cm2": Rct_area,
                "Cdl_F_cm2": Cdl_per_cm2,
                "f_peak_Hz": f_peak,
                "delta_cm": delta_cm,
                "D_used_cm2_s": D_calc if np.isfinite(D_calc) else np.nan,
            })

        ax.set_xlabel(r"Z' (Ω·cm$^2$)")
        ax.set_ylabel(r"-Z'' (Ω·cm$^2$)")
        ax.set_title("Task 3.3 – PEIS Nyquist (area-normalized)")
        ax.legend(title="Rotation rate", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
        beautify_axes(ax)
        ax.set_aspect("equal", adjustable="box")
        safe_save(fig, "T3.3_EIS_Nyquist.png")

        if peis_rows:
            write_csv(pd.DataFrame(peis_rows), "T3.3_PEIS_summary.csv")
    else:
        print("[PEIS] No PEIS files found; skipping PEIS plot.")


# =============================================================
# Main
# =============================================================

def main() -> None:
    setup_plot_style()
    print(f"Working directory: {ROOT}")
    print(f"Electrode area (GC, d={GC_DIAMETER_MM} mm): {GC_AREA_CM2:.4f} cm^2")

    # NOTE: Task 3.1 (coin cell) data not yet collected -> not analyzed here.

    # Task 3.2 – Interface 1
    analyze_task32_cvs()             # CV overlay + peaks + Randles–Sevcik D
    analyze_task32_cp()              # CP E(t) + τ, E_{τ/4}, Sand D
    analyze_task32_eis()             # Nyquist + (Rs, Rct, Cdl, k0) + Warburg D

    # Task 3.3 – Interface 2
    analyze_task33_lsv_and_eis()     # LSV overlay + K–L + Tafel + PEIS params


if __name__ == "__main__":
    main()
