import math
from dataclasses import dataclass
from typing import List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from .config import (
    C_BULK_MOL_PER_CM3,
    DATA_DIR,
    GC_AREA_CM2,
    N_ELECTRONS,
    R_GAS,
    T_K,
    F_CONST,
)
from .diffusion import add_diffusion_coefficient
from .io_utils import read_biologic_table, write_csv
from .plotting import beautify_axes, safe_save
from .utils import get_column, label_from_filename

try:
    from scipy.optimize import least_squares  # type: ignore
except Exception:  # pragma: no cover - allow running without SciPy
    least_squares = None  # type: ignore

try:
    from impedance import preprocessing
    from impedance.models.circuits import CustomCircuit
    from impedance.models.circuits.fitting import calculateCircuitLength
    from impedance.validation import linKK
except ImportError:
    preprocessing = None
    CustomCircuit = None
    calculateCircuitLength = None
    linKK = None

@dataclass
class ECFFit:
    model: str
    Rs: float
    Rct: float
    Cdl: Optional[float]
    sigma: Optional[float]
    Q: Optional[float]
    alpha: Optional[float]
    rchi2: float
    n_points: int
    raw_parameters: Optional[np.ndarray] = None
    fitted_freq_hz: Optional[np.ndarray] = None
    fitted_impedance: Optional[np.ndarray] = None

def fit_with_impedance_py(df: pd.DataFrame, circuit_str: str = 'R0-p(R1,C1)-W1', initial_guess: Optional[List[float]] = None) -> ECFFit:
    if CustomCircuit is None:
        return ECFFit(
            model=circuit_str,
            Rs=float('nan'),
            Rct=float('nan'),
            Cdl=None,
            sigma=None,
            Q=None,
            alpha=None,
            rchi2=float('nan'),
            n_points=0,
            raw_parameters=None,
            fitted_freq_hz=None,
            fitted_impedance=None,
        )

    re_col = get_column(df, ['Re(Z)/Ohm']) or 'Re(Z)/Ohm'
    im_col = get_column(df, ['-Im(Z)/Ohm']) or '-Im(Z)/Ohm'
    f_col = get_column(df, ['freq/Hz']) or 'freq/Hz'

    f = pd.to_numeric(df[f_col], errors='coerce').to_numpy()
    Zre = pd.to_numeric(df[re_col], errors='coerce').to_numpy()
    Zim = -pd.to_numeric(df[im_col], errors='coerce').to_numpy()
    mask = np.isfinite(f) & np.isfinite(Zre) & np.isfinite(Zim) & (f > 0)
    if not mask.any():
        return ECFFit(
            model=circuit_str,
            Rs=float('nan'),
            Rct=float('nan'),
            Cdl=None,
            sigma=None,
            Q=None,
            alpha=None,
            rchi2=float('nan'),
            n_points=0,
            raw_parameters=None,
            fitted_freq_hz=None,
            fitted_impedance=None,
        )

    order = np.argsort(f[mask])[::-1]
    f_use = f[mask][order].astype(float)
    Zre_use = Zre[mask][order].astype(float)
    Zim_use = Zim[mask][order].astype(float)
    Z_use = Zre_use + 1j * Zim_use

    hf_window = max(3, len(Zre_use) // 6)
    rs_guess = float(np.nanmedian(Zre_use[:hf_window])) if len(Zre_use) else float('nan')
    lf_guess = float(np.nanmedian(Zre_use[-hf_window:])) if len(Zre_use) else float('nan')
    if np.isfinite(rs_guess) and np.isfinite(lf_guess):
        rct_guess = float(max(lf_guess - rs_guess, 1.0))
    else:
        rct_guess = 100.0

    im_mag = np.abs(Zim_use)
    idx_peak = int(np.nanargmax(im_mag)) if len(im_mag) else 0
    f_peak = float(f_use[idx_peak]) if 0 <= idx_peak < len(f_use) else float('nan')
    if np.isfinite(f_peak) and f_peak > 0 and rct_guess > 0:
        c_guess = float(1.0 / (2.0 * math.pi * f_peak * rct_guess))
    else:
        c_guess = 1e-5

    omega_low = 2.0 * math.pi * max(f_use[-1], 1e-3)
    sigma_base = abs(lf_guess - rs_guess) if np.isfinite(lf_guess) and np.isfinite(rs_guess) else rct_guess
    sigma_guess = float(max(sigma_base, 1.0) / max(np.sqrt(omega_low), 1e-6))
    sigma_guess = max(sigma_guess, 1.0)

    def _fallback_result() -> ECFFit:
        rs_area = rs_guess * GC_AREA_CM2 if np.isfinite(rs_guess) else float('nan')
        rct_area = rct_guess * GC_AREA_CM2 if np.isfinite(rct_guess) else float('nan')
        cdl_area = c_guess / GC_AREA_CM2 if np.isfinite(c_guess) and c_guess > 0 else None
        sigma_area = sigma_guess * GC_AREA_CM2 if np.isfinite(sigma_guess) else None
        return ECFFit(
            model=circuit_str,
            Rs=rs_area,
            Rct=rct_area,
            Cdl=cdl_area,
            sigma=sigma_area,
            Q=None,
            alpha=None,
            rchi2=float('nan'),
            n_points=len(f_use),
            raw_parameters=None,
            fitted_freq_hz=None,
            fitted_impedance=None,
        )

    guess = list(initial_guess) if initial_guess else []
    if not guess:
        if circuit_str == 'R0-p(R1,C1)-W1' or ('W' in circuit_str and 'p(R1,C1)' in circuit_str):
            guess = [
                float(max(rs_guess, 1e-3)),
                float(max(rct_guess, 1.0)),
                float(max(c_guess, 1e-8)),
                float(max(sigma_guess, 1.0)),
            ]
        elif circuit_str == 'R0-p(R1,CPE1)-W1':
            # Randles with CPE: Rs - p(Rct, CPE) - Warburg
            alpha_guess = 0.85
            # Use peak frequency if available; otherwise a HF-mid value
            f_ref = f_peak if ('f_peak' in locals() and np.isfinite(f_peak) and f_peak > 0) else (f_use[len(f_use)//4] if len(f_use) else 1.0)
            w_ref = 2.0 * math.pi * max(float(f_ref), 1e-12)
            q_guess = float(max(min(1.0 / (max(rct_guess, 1.0) * (w_ref ** alpha_guess)), 1e-1), 1e-10))
            guess = [
                float(max(rs_guess, 1e-3)),
                float(max(rct_guess, 1.0)),
                q_guess,
                float(min(max(alpha_guess, 0.5), 0.99)),
                float(max(sigma_guess, 1.0)),
            ]
        elif circuit_str == 'R0-p(R1,C1)-p(R2,C2)':
            # Two-semicircle heuristic from Z''(f) peaks
            logf = np.log10(np.maximum(f_use, 1e-12))
            # Work with -Im(Z) which is positive for capacitive arcs
            zim_pos = -Zim_use
            zim_smooth = pd.Series(zim_pos).rolling(window=max(5, len(zim_pos) // 50), center=True, min_periods=1).median().to_numpy()

            # Find local maxima indices in smoothed Zim
            candidate_idx = []
            for i in range(1, len(zim_smooth) - 1):
                if zim_smooth[i] >= zim_smooth[i - 1] and zim_smooth[i] >= zim_smooth[i + 1]:
                    candidate_idx.append(i)

            if not candidate_idx and len(zim_smooth):
                candidate_idx = [int(np.nanargmax(zim_smooth))]

            # Filter by amplitude threshold and enforce separation in decades
            peaks: List[int] = []
            if candidate_idx:
                amp_thr = 0.10 * float(np.nanmax(zim_smooth))
                cand_sorted = sorted(candidate_idx, key=lambda k: float(zim_smooth[k]), reverse=True)
                min_sep_dec = 0.20
                for idx in cand_sorted:
                    if not np.isfinite(zim_smooth[idx]) or zim_smooth[idx] < amp_thr:
                        continue
                    if all(abs(logf[idx] - logf[p]) >= min_sep_dec for p in peaks):
                        peaks.append(idx)
                    if len(peaks) >= 2:
                        break

            # Build guesses
            R1_guess, R2_guess = float(max(rct_guess * 0.6, 1.0)), float(max(rct_guess * 0.4, 1.0))
            C1_guess, C2_guess = float(max(c_guess, 1e-8)), float(max(c_guess, 1e-8))

            if len(peaks) >= 1:
                f_pk1 = float(f_use[peaks[0]])
                tau1 = 1.0 / (2.0 * math.pi * max(f_pk1, 1e-12))
                R1_pre = float(max(2.0 * zim_smooth[peaks[0]], 1e-6))
            else:
                f_pk1 = float(f_peak) if ("f_peak" in locals() and np.isfinite(f_peak)) else float('nan')
                tau1 = 1.0 / (2.0 * math.pi * max(f_pk1, 1e-12)) if np.isfinite(f_pk1) else 1.0
                R1_pre = float(max(rct_guess * 0.7, 1.0))

            if len(peaks) >= 2:
                f_pk2 = float(f_use[peaks[1]])
                tau2 = 1.0 / (2.0 * math.pi * max(f_pk2, 1e-12))
                R2_pre = float(max(2.0 * zim_smooth[peaks[1]], 1e-6))
            else:
                # Place a second arc roughly one decade lower in frequency
                tau2 = tau1 * 10.0
                R2_pre = float(max(rct_guess * 0.3, 1.0))

            # Scale R1,R2 so that R1+R2 ~= total diameter (lf-hf)
            r_sum_pre = R1_pre + R2_pre
            if np.isfinite(r_sum_pre) and r_sum_pre > 0 and np.isfinite(rct_guess) and rct_guess > 0:
                scale_r = float(rct_guess / r_sum_pre)
                R1_guess = float(max(R1_pre * scale_r, 1.0))
                R2_guess = float(max(R2_pre * scale_r, 1.0))

            # Capacitances from tau = R*C
            C1_guess = float(max(min(tau1 / max(R1_guess, 1e-6), 1e-1), 1e-9))
            C2_guess = float(max(min(tau2 / max(R2_guess, 1e-6), 1e-1), 1e-9))

            # Final guardrails
            R1_guess = float(max(min(R1_guess, 1e7), 1e-3))
            R2_guess = float(max(min(R2_guess, 1e7), 1e-3))

            guess = [
                float(max(rs_guess, 1e-3)),
                R1_guess,
                C1_guess,
                R2_guess,
                C2_guess,
            ]
        elif circuit_str == 'R0-p(R1,CPE1)-p(R2,CPE2)':
            # Two CPE arcs: estimate R1,R2 from peaks; alpha~0.85; Q from f_peak
            logf = np.log10(np.maximum(f_use, 1e-12))
            zim_pos = -Zim_use
            zim_smooth = pd.Series(zim_pos).rolling(window=max(5, len(zim_pos) // 50), center=True, min_periods=1).median().to_numpy()

            candidate_idx = []
            for i in range(1, len(zim_smooth) - 1):
                if zim_smooth[i] >= zim_smooth[i - 1] and zim_smooth[i] >= zim_smooth[i + 1]:
                    candidate_idx.append(i)
            if not candidate_idx and len(zim_smooth):
                candidate_idx = [int(np.nanargmax(zim_smooth))]

            peaks: List[int] = []
            if candidate_idx:
                amp_thr = 0.10 * float(np.nanmax(zim_smooth))
                cand_sorted = sorted(candidate_idx, key=lambda k: float(zim_smooth[k]), reverse=True)
                min_sep_dec = 0.20
                for idx in cand_sorted:
                    if not np.isfinite(zim_smooth[idx]) or zim_smooth[idx] < amp_thr:
                        continue
                    if all(abs(logf[idx] - logf[p]) >= min_sep_dec for p in peaks):
                        peaks.append(idx)
                    if len(peaks) >= 2:
                        break

            # Rough R estimates from peak heights and total diameter scaling
            if len(peaks) >= 1:
                f_pk1 = float(f_use[peaks[0]])
                R1_pre = float(max(2.0 * zim_smooth[peaks[0]], 1e-6))
            else:
                f_pk1 = float(f_use[int(np.argmax(zim_smooth))]) if len(zim_smooth) else float('nan')
                R1_pre = float(max(rct_guess * 0.7, 1.0))
            if len(peaks) >= 2:
                f_pk2 = float(f_use[peaks[1]])
                R2_pre = float(max(2.0 * zim_smooth[peaks[1]], 1e-6))
            else:
                f_pk2 = f_pk1 / 10.0 if np.isfinite(f_pk1) and f_pk1 > 0 else 1.0
                R2_pre = float(max(rct_guess * 0.3, 1.0))

            r_sum_pre = R1_pre + R2_pre
            if np.isfinite(r_sum_pre) and r_sum_pre > 0 and np.isfinite(rct_guess) and rct_guess > 0:
                scale_r = float(rct_guess / r_sum_pre)
                R1_guess = float(max(R1_pre * scale_r, 1.0))
                R2_guess = float(max(R2_pre * scale_r, 1.0))
            else:
                R1_guess = float(max(rct_guess * 0.6, 1.0))
                R2_guess = float(max(rct_guess * 0.4, 1.0))

            alpha1_guess = 0.85
            alpha2_guess = 0.85
            # Q from (omega_pk)^alpha = 1/(R*Q)  => Q = 1/(R * omega_pk^alpha)
            w1 = 2.0 * math.pi * max(f_pk1, 1e-12)
            w2 = 2.0 * math.pi * max(f_pk2, 1e-12)
            Q1_guess = float(max(min(1.0 / (max(R1_guess, 1e-6) * (w1 ** alpha1_guess)), 1e-1), 1e-10))
            Q2_guess = float(max(min(1.0 / (max(R2_guess, 1e-6) * (w2 ** alpha2_guess)), 1e-1), 1e-10))

            # Guardrails
            R1_guess = float(max(min(R1_guess, 1e7), 1e-3))
            R2_guess = float(max(min(R2_guess, 1e7), 1e-3))

            guess = [
                float(max(rs_guess, 1e-3)),
                R1_guess,
                Q1_guess,
                float(min(max(alpha1_guess, 0.5), 0.99)),
                R2_guess,
                Q2_guess,
                float(min(max(alpha2_guess, 0.5), 0.99)),
            ]
        else:
            guess = [float(max(rs_guess, 1e-3))]

    param_len = int(calculateCircuitLength(circuit_str)) if calculateCircuitLength is not None else len(guess)
    if len(guess) < param_len:
        guess.extend([float(max(c_guess, 1e-8))] * (param_len - len(guess)))
    elif len(guess) > param_len and param_len > 0:
        guess = guess[:param_len]

    # Construct simple physical bounds to stabilize the fit
    bounds = None
    try:
        if circuit_str == 'R0-p(R1,C1)-W1':
            rs_lo = max(1e-6, rs_guess * 0.01 if np.isfinite(rs_guess) else 1e-6)
            rct_lo = 1e-6
            c_lo = 1e-10
            sigma_lo = 1e-6
            rs_hi = max(rs_lo * 100.0, 1e3)
            rct_hi = max(rct_guess * 100.0, 1e3)
            c_hi = 1e-1
            sigma_hi = 1e3
            bounds = ([rs_lo, rct_lo, c_lo, sigma_lo], [rs_hi, rct_hi, c_hi, sigma_hi])
        elif circuit_str == 'R0-p(R1,CPE1)-W1':
            rs_lo = max(1e-6, rs_guess * 0.01 if np.isfinite(rs_guess) else 1e-6)
            r_lo = 1e-6
            q_lo = 1e-10
            a_lo = 0.5
            sigma_lo = 1e-6
            rs_hi = max(rs_lo * 100.0, 1e6)
            r_hi = max(rct_guess * 100.0, 1e8)
            q_hi = 1e-1
            a_hi = 0.99
            sigma_hi = 1e3
            bounds = ([rs_lo, r_lo, q_lo, a_lo, sigma_lo], [rs_hi, r_hi, q_hi, a_hi, sigma_hi])
        elif circuit_str == 'R0-p(R1,C1)-p(R2,C2)':
            rs_lo = max(1e-6, rs_guess * 0.01 if np.isfinite(rs_guess) else 1e-6)
            r_lo = 1e-6
            c_lo = 1e-10
            rs_hi = max(rs_lo * 100.0, 1e6)
            r_hi = max(rct_guess * 100.0, 1e8)
            c_hi = 1e-1
            bounds = ([rs_lo, r_lo, c_lo, r_lo, c_lo], [rs_hi, r_hi, c_hi, r_hi, c_hi])
        elif circuit_str == 'R0-p(R1,CPE1)-p(R2,CPE2)':
            rs_lo = max(1e-6, rs_guess * 0.01 if np.isfinite(rs_guess) else 1e-6)
            r_lo = 1e-6
            q_lo = 1e-10
            a_lo = 0.5
            rs_hi = max(rs_lo * 100.0, 1e6)
            r_hi = max(rct_guess * 100.0, 1e8)
            q_hi = 1e-1
            a_hi = 0.99
            bounds = ([rs_lo, r_lo, q_lo, a_lo, r_lo, q_lo, a_lo], [rs_hi, r_hi, q_hi, a_hi, r_hi, q_hi, a_hi])
    except Exception:
        bounds = None

    try:
        circuit = CustomCircuit(initial_guess=guess, circuit=circuit_str)
    except Exception:
        return _fallback_result()

    try:
        # Frequency weighting: down-weight low-frequency large-|Z| points to help resolve HF arc
        weights = None
        try:
            mag = np.abs(Z_use)
            weights = 1.0 / np.maximum(mag, np.median(mag) * 1e-3)
        except Exception:
            weights = None

        if bounds is not None:
            try:
                if weights is not None:
                    fitted = circuit.fit(f_use, Z_use, bounds=bounds, weights=weights)
                else:
                    fitted = circuit.fit(f_use, Z_use, bounds=bounds)
            except TypeError:
                fitted = circuit.fit(f_use, Z_use)
        else:
            try:
                if weights is not None:
                    fitted = circuit.fit(f_use, Z_use, weights=weights)
                else:
                    fitted = circuit.fit(f_use, Z_use)
            except TypeError:
                fitted = circuit.fit(f_use, Z_use)
    except Exception:
        return _fallback_result()

    params_raw = np.array(fitted.parameters_, dtype=float)
    Z_model = fitted.predict(f_use)

    res_real = Z_model.real - Z_use.real
    res_imag = Z_model.imag - Z_use.imag
    ss_res = float(np.sum(res_real**2 + res_imag**2))
    dof = max(1, 2 * len(f_use) - params_raw.size)
    rchi2 = ss_res / dof

    Rs_raw = float(params_raw[0]) if params_raw.size >= 1 else float('nan')
    # Map parameters by circuit
    if circuit_str == 'R0-p(R1,C1)-W1':
        Rct_raw = float(params_raw[1]) if params_raw.size >= 2 else float('nan')
        Cdl_raw = float(params_raw[2]) if params_raw.size >= 3 else None
        sigma_raw = float(params_raw[3]) if params_raw.size >= 4 else None
        Q_raw, alpha_raw = None, None
    elif circuit_str == 'R0-p(R1,CPE1)-W1':
        Rct_raw = float(params_raw[1]) if params_raw.size >= 2 else float('nan')
        Q_raw = float(params_raw[2]) if params_raw.size >= 3 else None
        alpha_raw = float(params_raw[3]) if params_raw.size >= 4 else None
        sigma_raw = float(params_raw[4]) if params_raw.size >= 5 else None
        Cdl_raw = None
    elif circuit_str == 'R0-p(R1,C1)-p(R2,C2)':
        r1 = float(params_raw[1]) if params_raw.size >= 2 else float('nan')
        c1 = float(params_raw[2]) if params_raw.size >= 3 else float('nan')
        r2 = float(params_raw[3]) if params_raw.size >= 4 else float('nan')
        c2 = float(params_raw[4]) if params_raw.size >= 5 else float('nan')
        Rct_raw = (r1 + r2) if np.isfinite(r1) and np.isfinite(r2) else float('nan')
        Cdl_raw = c1  # store first branch C as representative Cdl
        sigma_raw = None
        Q_raw, alpha_raw = None, None
    elif circuit_str == 'R0-p(R1,CPE1)-p(R2,CPE2)':
        r1 = float(params_raw[1]) if params_raw.size >= 2 else float('nan')
        q1 = float(params_raw[2]) if params_raw.size >= 3 else float('nan')
        a1 = float(params_raw[3]) if params_raw.size >= 4 else float('nan')
        r2 = float(params_raw[4]) if params_raw.size >= 5 else float('nan')
        q2 = float(params_raw[5]) if params_raw.size >= 6 else float('nan')
        a2 = float(params_raw[6]) if params_raw.size >= 7 else float('nan')
        Rct_raw = (r1 + r2) if np.isfinite(r1) and np.isfinite(r2) else float('nan')
        Cdl_raw = None
        sigma_raw = None
        # store first branch Q, alpha as representatives
        Q_raw = q1
        alpha_raw = a1
    else:
        Rct_raw = float(params_raw[1]) if params_raw.size >= 2 else float('nan')
        Cdl_raw = float(params_raw[2]) if params_raw.size >= 3 else None
        sigma_raw = float(params_raw[3]) if ('W' in circuit_str and params_raw.size >= 4) else None
        Q_raw, alpha_raw = None, None

    scale = GC_AREA_CM2
    Rs_area = Rs_raw * scale if np.isfinite(Rs_raw) else float('nan')
    Rct_area = Rct_raw * scale if np.isfinite(Rct_raw) else float('nan')
    Cdl_area = Cdl_raw / scale if (Cdl_raw is not None and np.isfinite(Cdl_raw) and Cdl_raw > 0) else None
    sigma_area = sigma_raw * scale if (sigma_raw is not None and np.isfinite(sigma_raw)) else None
    Q_area = Q_raw / scale if (Q_raw is not None and np.isfinite(Q_raw) and Q_raw > 0) else None
    alpha_val = float(alpha_raw) if (alpha_raw is not None and np.isfinite(alpha_raw)) else None

    fit_impedance_area = Z_model * scale

    result = ECFFit(
        model=circuit_str,
        Rs=Rs_area,
        Rct=Rct_area,
        Cdl=Cdl_area,
        sigma=sigma_area,
        Q=Q_area,
        alpha=alpha_val,
        rchi2=float(rchi2),
        n_points=len(f_use),
        raw_parameters=params_raw,
        fitted_freq_hz=f_use,
        fitted_impedance=fit_impedance_area,
    )
    result.parameters_ = params_raw
    return result

@dataclass
class EISResult:
    label: str
    Rs_area_fit_Ohm_cm2: float
    Rct_area_fit_Ohm_cm2: float
    Rs_area_heur_Ohm_cm2: float
    Rct_area_heur_Ohm_cm2: float
    Cdl_F_cm2: float
    f_peak_Hz: float
    k0_cm_s: float
    Aw_Ohm_cm2_sqrt_s: float
    D_Warburg_cm2_s: float
    Rs_area_nls_Ohm_cm2: float = float("nan")
    Rct_area_nls_Ohm_cm2: float = float("nan")
    Cdl_nls_F_cm2: float = float("nan")
    sigma_Ohm_cm2_s05: float = float("nan")
    rchi2: float = float("nan")
    model_label: str = ""


def _identify_semicircle_region(
    z_re: np.ndarray,
    z_im: np.ndarray,
    f_hz: np.ndarray,
) -> Tuple[int, int, int, float, float]:
    """Identify indices bounding the primary semicircle and estimate Rs and f_peak."""
    smoothed = pd.Series(z_im).rolling(window=max(5, len(z_im) // 100), center=True, min_periods=1).median().to_numpy()
    d_im = np.diff(smoothed)
    peak_candidates = np.where((d_im[:-1] > 0) & (d_im[1:] <= 0))[0] + 1
    idx_peak1 = int(peak_candidates[0]) if len(peak_candidates) else int(np.nanargmax(smoothed))

    hf_cut = int(max(3, round(0.15 * len(z_re))))
    Rs_area = float(np.nanmedian(z_re[:hf_cut])) if hf_cut > 0 else float(np.nanmin(z_re))

    im_peak = float(smoothed[idx_peak1])
    thresh = 0.05 * im_peak
    hf_below_idx = np.where(smoothed[:idx_peak1] <= thresh)[0]
    start_idx = int(hf_below_idx[-1]) if len(hf_below_idx) else 0
    lf_below_idx = np.where(smoothed[idx_peak1 + 1 :] <= thresh)[0]
    if len(lf_below_idx):
        end_idx = int(idx_peak1 + 1 + lf_below_idx[0])
    else:
        after = d_im[idx_peak1:]
        pos_turn = np.where(after > 0)[0]
        end_idx = int(idx_peak1 + (pos_turn[0] if len(pos_turn) else max(3, len(z_re) // 3)))
        end_idx = min(end_idx, len(z_re) - 1)

    f_peak = float(f_hz[idx_peak1]) if np.isfinite(f_hz[idx_peak1]) and (f_hz[idx_peak1] > 0) else float("nan")
    return idx_peak1, start_idx, end_idx, Rs_area, f_peak


def _extract_eis_parameters_heuristic(
    df: pd.DataFrame,
) -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, Tuple[int, int, int]]:
    re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
    im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
    f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
    z_re_raw = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
    z_im_raw = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
    f_hz_raw = pd.to_numeric(df[f_col], errors="coerce").to_numpy()

    mask = np.isfinite(z_re_raw) & np.isfinite(z_im_raw) & np.isfinite(f_hz_raw) & (f_hz_raw > 0)
    if not mask.any():
        return float("nan"), float("nan"), float("nan"), float("nan"), z_re_raw, z_im_raw, f_hz_raw, (0, 0, 0)

    order = np.argsort(f_hz_raw[mask])[::-1]
    z_re = z_re_raw[mask][order]
    z_im = z_im_raw[mask][order]
    f_hz = f_hz_raw[mask][order]

    idx_peak1, start_idx, end_idx, Rs_area, f_peak = _identify_semicircle_region(z_re, z_im, f_hz)

    if end_idx <= start_idx:
        Rct_area_est = 2.0 * float(z_im[idx_peak1])
        Cdl_per_cm2 = float(1.0 / (2.0 * math.pi * f_peak * Rct_area_est)) if (np.isfinite(f_peak) and Rct_area_est > 0) else float("nan")
        return Rs_area, max(Rct_area_est, 0.0), Cdl_per_cm2, f_peak, z_re, z_im, f_hz, (idx_peak1, start_idx, end_idx)

    seg_len = end_idx - start_idx + 1
    left_span = max(1, int(0.15 * seg_len))
    right_span = max(1, int(0.15 * seg_len))
    Rs_est = float(np.nanmedian(z_re[start_idx : start_idx + left_span]))
    Rlf_candidates = z_re[end_idx - right_span + 1 : end_idx + 1]
    R_low_area = float(np.nanmedian(Rlf_candidates)) if len(Rlf_candidates) else float(z_re[end_idx])
    Rct_area = max(R_low_area - Rs_est, 0.0)
    Cdl_per_cm2 = float(1.0 / (2.0 * math.pi * f_peak * Rct_area)) if (np.isfinite(f_peak) and Rct_area > 0) else float("nan")
    return float(Rs_est), float(Rct_area), float(Cdl_per_cm2), f_peak, z_re, z_im, f_hz, (idx_peak1, start_idx, end_idx)


def _circle_fit_semicircle(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    try:
        A = np.column_stack([x, y, np.ones_like(x)])
        b = -(x ** 2 + y ** 2)
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        a, bcoef, c = coef
        x0 = -a / 2.0
        y0 = -bcoef / 2.0
        r_sq = x0 ** 2 + y0 ** 2 - c
        radius = float(np.sqrt(r_sq)) if r_sq > 0 else float("nan")
        return float(x0), float(y0), radius
    except Exception:
        return float("nan"), float("nan"), float("nan")


def _extract_eis_parameters_both(df: pd.DataFrame) -> Tuple[float, float, float, float, float, float]:
    Rs_h, Rct_h, Cdl_h, f_peak, z_re, z_im, f_hz, (idx_peak1, start_idx, end_idx) = _extract_eis_parameters_heuristic(df)

    Rs_fit, Rct_fit = float("nan"), float("nan")
    if (end_idx > start_idx) and (end_idx - start_idx + 1) >= 5 and np.isfinite(z_re[start_idx : end_idx + 1]).all() and np.isfinite(z_im[start_idx : end_idx + 1]).all():
        x_arc = z_re[start_idx : end_idx + 1]
        y_arc = z_im[start_idx : end_idx + 1]
        x0, y0, radius = _circle_fit_semicircle(x_arc, y_arc)
        if np.isfinite(x0) and np.isfinite(y0) and np.isfinite(radius) and radius > 0:
            delta_sq = radius ** 2 - y0 ** 2
            if delta_sq >= 0:
                delta = float(np.sqrt(delta_sq))
                x_left, x_right = x0 - delta, x0 + delta
                Rs_fit = float(min(x_left, x_right))
                Rct_fit = float(max(x_left, x_right) - min(x_left, x_right))

    Rs_best, Rct_best = Rs_h, Rct_h
    zim_max = float(np.nanmax(z_im)) if len(z_im) else float("nan")
    if np.isfinite(Rct_fit) and Rct_fit > 0:
        if not np.isfinite(zim_max) or (0.3 <= (Rct_fit / max(zim_max * 2.0, 1e-9)) <= 3.0):
            Rs_best, Rct_best = Rs_fit, Rct_fit

    Cdl_best = float(1.0 / (2.0 * math.pi * f_peak * Rct_best)) if (np.isfinite(f_peak) and Rct_best and Rct_best > 0) else Cdl_h
    return float(Rs_best), float(Rct_best), float(Cdl_best), float(f_peak), float(Rs_fit), float(Rct_fit)


def _warburg_fit(df: pd.DataFrame) -> Tuple[float, float]:
    re_col = get_column(df, ["Re(Z)/Ohm"]) or "Re(Z)/Ohm"
    im_col = get_column(df, ["-Im(Z)/Ohm"]) or "-Im(Z)/Ohm"
    f_col = get_column(df, ["freq/Hz"]) or "freq/Hz"
    z_re = pd.to_numeric(df[re_col], errors="coerce").to_numpy() * GC_AREA_CM2
    z_im = pd.to_numeric(df[im_col], errors="coerce").to_numpy() * GC_AREA_CM2
    f_hz = pd.to_numeric(df[f_col], errors="coerce").to_numpy()
    omega = 2 * np.pi * f_hz

    mask = np.isfinite(z_re) & np.isfinite(z_im) & np.isfinite(omega) & (omega > 0)
    z_re, z_im, omega = z_re[mask], z_im[mask], omega[mask]
    if len(omega) == 0:
        return float("nan"), float("nan")

    angle45 = np.abs(z_re - z_im) / np.maximum(1e-9, (np.abs(z_re) + np.abs(z_im))) < 0.20
    lf_cut = np.quantile(omega, 0.35)
    sel = (omega <= lf_cut) & angle45

    def _fit_slope(x_vals: np.ndarray, y_vals: np.ndarray) -> float:
        A_loc = np.vstack([x_vals, np.ones_like(x_vals)]).T
        slope_loc, _ = np.linalg.lstsq(A_loc, y_vals, rcond=None)[0]
        return float(slope_loc)

    if sel.sum() < 5:
        order = np.argsort(omega)
        num = max(6, min(20, int(len(omega) * 0.30)))
        idx = order[:num]
    else:
        idx = np.where(sel)[0]

    x = 1.0 / np.sqrt(omega[idx])
    y_re = z_re[idx]
    y_im = z_im[idx]

    slope_re = _fit_slope(x, y_re) if len(x) >= 2 else float("nan")
    slope_im = _fit_slope(x, y_im) if len(x) >= 2 else float("nan")

    slopes = [s for s in [slope_re, slope_im] if np.isfinite(s) and s > 0]
    if len(slopes) == 0:
        return float("nan"), float("nan")

    if len(slopes) == 2 and (max(slopes) / min(slopes) > 2.0):
        Aw = float(min(slopes))
    else:
        Aw = float(np.mean(slopes))

    D_half = (4.0 * R_GAS * T_K) / (((N_ELECTRONS ** 2) * (F_CONST ** 2) * C_BULK_MOL_PER_CM3) * Aw)
    D_w = float(D_half ** 2) if np.isfinite(D_half) and D_half > 0 else float("nan")
    return Aw, D_w


def analyze_task32_eis() -> None:
    eis_dir = DATA_DIR / "Task 3.2 EIS"
    files = sorted(eis_dir.glob("*_C01.txt"))
    if not files:
        print("[EIS] No EIS files found; skipping.")
        return

    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    results: List[EISResult] = []

    from pathlib import Path
    debug_dir = Path("results/debugplots")
    debug_dir.mkdir(exist_ok=True, parents=True)

    for path in files:
        try:
            df = read_biologic_table(path)
        except Exception as exc:
            print(f"[EIS] Failed to read {path.name}: {exc}")
            continue

        re_col = get_column(df, ['Re(Z)/Ohm']) or 'Re(Z)/Ohm'
        im_col = get_column(df, ['-Im(Z)/Ohm']) or '-Im(Z)/Ohm'
        f_col = get_column(df, ['freq/Hz']) or 'freq/Hz'
        if re_col not in df or im_col not in df:
            print(f"[EIS] Missing impedance columns in {path.name}; skipping.")
            continue

        label = label_from_filename(path)
        z_re = pd.to_numeric(df[re_col], errors='coerce').to_numpy() * GC_AREA_CM2
        z_im = pd.to_numeric(df[im_col], errors='coerce').to_numpy() * GC_AREA_CM2
        ax.plot(z_re, z_im, marker='o', ms=2.5, lw=0.8, label=label)

        chosen = fit_with_impedance_py(df, circuit_str='R0-p(R1,CPE1)-W1')
        z_fit_plot = None
        if chosen.fitted_impedance is not None:
            z_fit = np.asarray(chosen.fitted_impedance)
            finite_mask = np.isfinite(z_fit.real) & np.isfinite(z_fit.imag)
            if finite_mask.any():
                z_fit_plot = z_fit[finite_mask]
                ax.plot(z_fit_plot.real, -z_fit_plot.imag, lw=1.2, alpha=0.9, label=f"{label} (fit)")

        fig_debug, ax_debug = plt.subplots(figsize=(5.0, 4.0))
        ax_debug.plot(z_re, z_im, 'o', ms=3, label='data')
        if z_fit_plot is not None:
            ax_debug.plot(z_fit_plot.real, -z_fit_plot.imag, '-', lw=1.5, label='fit')

        ax_debug.set_xlabel("Z' (ohm cm^2)")
        ax_debug.set_ylabel("-Z'' (ohm cm^2)")
        ax_debug.set_title(f"Debug EIS Fit: {label}")
        ax_debug.text(
            0.05,
            0.90,
            f'rchi2 = {chosen.rchi2:.2e}\nN = {chosen.n_points}',
            transform=ax_debug.transAxes,
            fontsize=10,
        )
        ax_debug.legend()
        beautify_axes(ax_debug)
        ax_debug.set_aspect('equal', adjustable='box')
        fig_debug.savefig(debug_dir / f'debug_EIS_fit_{label}.png', dpi=150, bbox_inches='tight')
        plt.close(fig_debug)

        Rs_heur, Rct_heur, Cdl_h, f_peak, z_re_h, z_im_h, f_hz_h, (idx_peak1, start_idx, end_idx) = _extract_eis_parameters_heuristic(df)

        Rs_fit, Rct_fit = float('nan'), float('nan')
        if (end_idx > start_idx) and (end_idx - start_idx + 1) >= 5 and np.isfinite(z_re_h[start_idx : end_idx + 1]).all() and np.isfinite(z_im_h[start_idx : end_idx + 1]).all():
            x_arc = z_re_h[start_idx : end_idx + 1]
            y_arc = z_im_h[start_idx : end_idx + 1]
            x0, y0, radius = _circle_fit_semicircle(x_arc, y_arc)
            if np.isfinite(x0) and np.isfinite(y0) and np.isfinite(radius) and radius > 0:
                delta_sq = radius ** 2 - y0 ** 2
                if delta_sq >= 0:
                    delta = float(np.sqrt(delta_sq))
                    x_left, x_right = x0 - delta, x0 + delta
                    Rs_fit = float(min(x_left, x_right))
                    Rct_fit = float(max(x_left, x_right) - min(x_left, x_right))

        Rs_best, Rct_best = Rs_heur, Rct_heur
        zim_max = float(np.nanmax(z_im_h)) if len(z_im_h) else float('nan')
        if np.isfinite(Rct_fit) and Rct_fit > 0:
            if not np.isfinite(zim_max) or (0.3 <= (Rct_fit / max(zim_max * 2.0, 1e-9)) <= 3.0):
                Rs_best, Rct_best = Rs_fit, Rct_fit
        Cdl_per_cm2 = float(1.0 / (2.0 * math.pi * f_peak * Rct_best)) if (np.isfinite(f_peak) and Rct_best and Rct_best > 0) else Cdl_h

        # Effective Cdl from CPE if present
        if (chosen.Cdl is not None) and np.isfinite(chosen.Cdl) and (chosen.Cdl > 0):
            cdl_nls_cm2 = float(chosen.Cdl)
        elif (
            (chosen.Q is not None)
            and (chosen.alpha is not None)
            and np.isfinite(chosen.Q)
            and np.isfinite(chosen.alpha)
            and (0.5 <= float(chosen.alpha) < 1.0)
            and np.isfinite(chosen.Rct)
            and (chosen.Rct > 0)
        ):
            a = float(chosen.alpha)
            q = float(chosen.Q)
            r = float(chosen.Rct)
            cdl_nls_cm2 = float((q ** (1.0 / a)) / (r ** (1.0 - 1.0 / a)))
        else:
            cdl_nls_cm2 = float('nan')

        Rct_for_k0 = chosen.Rct if np.isfinite(chosen.Rct) and (chosen.Rct > 0) else Rct_best
        k0 = (R_GAS * T_K) / ((N_ELECTRONS ** 2) * (F_CONST ** 2) * C_BULK_MOL_PER_CM3 * Rct_for_k0) if (Rct_for_k0 and Rct_for_k0 > 0) else np.nan

        Aw, D_w = _warburg_fit(df)
        add_diffusion_coefficient(D_w, 'Warburg')

        results.append(
            EISResult(
                label=label,
                Rs_area_fit_Ohm_cm2=Rs_fit,
                Rct_area_fit_Ohm_cm2=Rct_fit,
                Rs_area_heur_Ohm_cm2=Rs_heur,
                Rct_area_heur_Ohm_cm2=Rct_heur,
                Cdl_F_cm2=Cdl_per_cm2,
                f_peak_Hz=f_peak,
                k0_cm_s=float(k0),
                Aw_Ohm_cm2_sqrt_s=Aw,
                D_Warburg_cm2_s=D_w,
                Rs_area_nls_Ohm_cm2=chosen.Rs,
                Rct_area_nls_Ohm_cm2=chosen.Rct,
                Cdl_nls_F_cm2=cdl_nls_cm2,
                sigma_Ohm_cm2_s05=(chosen.sigma if chosen.sigma is not None else float('nan')),
                rchi2=chosen.rchi2,
                model_label=chosen.model,
            )
        )
    ax.set_xlabel("Z' (ohm cm^2)")
    ax.set_ylabel("-Z'' (ohm cm^2)")
    ax.set_title("Task 3.2 - EIS Nyquist (area-normalized)")
    ax.legend(title="Dataset", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8)
    beautify_axes(ax)
    ax.set_aspect("equal", adjustable="box")
    safe_save(fig, "T3.2_EIS_Nyquist.png")

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

    axw.set_xlabel("1/sqrt(omega) (s^0.5)")
    axw.set_ylabel("Z (ohm cm^2)")
    axw.set_title("Task 3.2 - Warburg linearization")
    axw.legend(title="Components", loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, fontsize=8, ncol=1)
    beautify_axes(axw)
    safe_save(figw, "T3.2_EIS_Warburg_linearization.png")

    if results:
        df_eis = pd.DataFrame(
            [
                {
                    "label": r.label,
                    "Rs_area_nls_Ohm_cm2": r.Rs_area_nls_Ohm_cm2,
                    "Rct_area_nls_Ohm_cm2": r.Rct_area_nls_Ohm_cm2,
                    "Cdl_nls_F_cm2": r.Cdl_nls_F_cm2,
                    "sigma_Ohm_cm2_s05": r.sigma_Ohm_cm2_s05,
                    "k0_cm_s": r.k0_cm_s,
                    "Aw_Ohm_cm2_sqrt_s": r.Aw_Ohm_cm2_sqrt_s,
                    "D_Warburg_cm2_s": r.D_Warburg_cm2_s,
                    "rchi2": r.rchi2,
                    "model_label": r.model_label,
                }
                for r in results
            ]
        )
        write_csv(df_eis, "T3.2_EIS_summary.csv")






