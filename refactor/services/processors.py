from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, List, Optional, Tuple

import numpy as np

from ..config import AnalysisConfig
from ..models import (
    ElectrochemicalMeasurement,
    CVData,
    CPData,
    EISData,
    PEISData,
    FitResult,
)

import pandas as pd


class DataProcessor(ABC):
    """Base class for all measurement processors."""

    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config

    @abstractmethod
    def process(self, measurement: ElectrochemicalMeasurement) -> FitResult | dict:
        ...

    def process_batch(self, measurements: Iterable[ElectrochemicalMeasurement]) -> List[FitResult | dict]:
        return [self.process(measurement) for measurement in measurements]


class CVProcessor(DataProcessor):
    def process(self, measurement: ElectrochemicalMeasurement) -> dict:
        assert isinstance(measurement, CVData)
        e_v = measurement.potential_V
        i_A = measurement.current_A

        ipa = float(np.nanmax(i_A))
        ipc = float(np.nanmin(i_A))
        Epa = float(e_v[np.nanargmax(i_A)])
        Epc = float(e_v[np.nanargmin(i_A)])
        delta_ep = Epa - Epc
        ratio = abs(ipa) / abs(ipc) if ipc != 0 else np.nan

        return {
            "ipa_A": ipa,
            "ipc_A": ipc,
            "Epa_V": Epa,
            "Epc_V": Epc,
            "delta_Ep_V": float(delta_ep),
            "ipa_ipc_ratio": float(ratio),
            "scan_rate_Vs": measurement.scan_rate_Vs if measurement.scan_rate_Vs is not None else float("nan"),
        }


class CPProcessor(DataProcessor):
    def process(self, measurement: ElectrochemicalMeasurement) -> dict:
        assert isinstance(measurement, CPData)
        t_s = measurement.time_s if measurement.time_s is not None else np.arange(len(measurement.potential_V), dtype=float)
        e_v = measurement.potential_V
        i_series = measurement.current_A * 1000.0  # convert back to mA

        segments = self._segment_by_current_plateaus(t_s, i_series)
        results = []
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

            tau = self._estimate_tau_s(t_seg, e_seg)
            if not np.isfinite(tau) or tau <= 0:
                tau = self._fallback_tau(t_seg, e_seg)

            if np.isfinite(tau) and tau > 0:
                e_tau4 = float(np.interp(tau / 4.0, t_seg, e_seg))
            else:
                e_tau4 = float("nan")

            if np.isfinite(i_applied_mA) and np.isfinite(tau) and tau > 0:
                D_half = (abs(i_applied_mA) * np.sqrt(tau)) / (85.5 * self.config.parameters.get("n_electrons", 1) * self.config.parameters.get("electrode_area_cm2", 1.0) * self.config.parameters.get("c_bulk_mM", 2.3))
                D_sand = float(D_half ** 2)
            else:
                D_sand = float("nan")

            results.append(
                {
                    "label": f"seg {seg_idx}",
                    "current_mA": i_applied_mA,
                    "tau_s": tau,
                    "E_tau_over_4_V": e_tau4,
                    "D_Sand_cm2_s": D_sand,
                }
            )

        return {"segments": results}

    def _segment_by_current_plateaus(self, t_s: np.ndarray, i_mA: Optional[np.ndarray]) -> List[Tuple[int, int]]:
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

    def _estimate_tau_s(self, t_s: np.ndarray, e_v: np.ndarray) -> float:
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

    def _fallback_tau(self, t_s: np.ndarray, e_v: np.ndarray) -> float:
        try:
            de_dt = np.gradient(pd.Series(e_v).rolling(11, center=True, min_periods=1).mean().to_numpy(), t_s)
            early = de_dt[: max(5, len(de_dt) // 10)]
            threshold = 3.0 * float(np.nanmedian(np.abs(early))) if np.isfinite(np.nanmedian(np.abs(early))) else float("nan")
            idx = int(np.argmax(de_dt > threshold)) if np.any(de_dt > threshold) else (len(t_s) - 1)
            return float(t_s[idx])
        except Exception:
            return float("nan")


class EISProcessor(DataProcessor):
    def process(self, measurement: ElectrochemicalMeasurement) -> dict:
        assert isinstance(measurement, EISData)
        max_im_idx = np.nanargmax(-measurement.impedance_imag_Ohm)
        f_peak = float(measurement.frequency_Hz[max_im_idx])
        return {
            "type": "eis",
            "f_peak_Hz": f_peak,
        }


class PEISProcessor(EISProcessor):
    def process(self, measurement: ElectrochemicalMeasurement) -> dict:
        result = super().process(measurement)
        assert isinstance(measurement, PEISData)
        result["rotation_speed_rpm"] = measurement.rotation_speed_rpm
        return result


