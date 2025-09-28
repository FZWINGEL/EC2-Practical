from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from ..config import AnalysisConfig
from ..models import FitResult


@dataclass(slots=True)
class TafelAnalyzer:
    config: AnalysisConfig

    def analyze(self, eta_V: np.ndarray, jk_A_cm2: np.ndarray) -> FitResult:
        mask = np.isfinite(eta_V) & np.isfinite(jk_A_cm2) & (np.abs(jk_A_cm2) > 1e-12)
        if mask.sum() < 3:
            return FitResult(model_name="Tafel", parameters={})
        x = np.log10(np.abs(jk_A_cm2[mask]))
        y = eta_V[mask]
        A = np.vstack([x, np.ones_like(x)]).T
        slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
        params = {
            "slope_V_per_decade": float(slope),
            "intercept_V": float(intercept),
            "j0_A_cm2": float(10 ** (-intercept / slope)) if slope != 0 else float("nan"),
        }
        return FitResult(model_name="Tafel", parameters=params)


