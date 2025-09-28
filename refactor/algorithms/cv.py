from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

import numpy as np

from ..config import AnalysisConfig
from ..models import CVData, FitResult


@dataclass(slots=True)
class RandlesSevcikAnalyzer:
    config: AnalysisConfig

    def analyze(self, cv_results: Iterable[dict]) -> FitResult:
        scan_rates = []
        peak_currents = []
        for result in cv_results:
            v = result.get("scan_rate_Vs")
            ip = max(abs(result.get("ipa_A", 0.0)), abs(result.get("ipc_A", 0.0)))
            if v is None or v <= 0 or not np.isfinite(ip):
                continue
            scan_rates.append(np.sqrt(v))
            peak_currents.append(ip)

        if len(scan_rates) < 2:
            return FitResult(model_name="Randles-Sevcik", parameters={}, goodness_of_fit=None)

        A = np.vstack([scan_rates, np.ones_like(scan_rates)]).T
        slope, intercept = np.linalg.lstsq(A, peak_currents, rcond=None)[0]
        params = {
            "slope_A_per_Vs_sqrt": float(slope),
            "intercept_A": float(intercept),
        }
        return FitResult(model_name="Randles-Sevcik", parameters=params, goodness_of_fit=float(intercept))


