from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List

import numpy as np

from ..config import AnalysisConfig
from ..models import FitResult


@dataclass(slots=True)
class KouteckyLevichAnalyzer:
    config: AnalysisConfig

    def analyze(self, j_by_rpm: Dict[int, np.ndarray], potentials: np.ndarray) -> FitResult:
        rpms = np.array(sorted(j_by_rpm.keys()))
        omega = 2 * np.pi * (rpms / 60.0)
        inv_sqrt_omega = 1.0 / np.sqrt(omega)
        rows = []
        for idx, E in enumerate(potentials):
            j_values = np.array([j_by_rpm[rpm][idx] for rpm in rpms])
            mask = np.isfinite(j_values) & (np.abs(j_values) > 1e-9)
            if mask.sum() < 3:
                continue
            x = inv_sqrt_omega[mask]
            y = 1.0 / j_values[mask]
            A = np.vstack([x, np.ones_like(x)]).T
            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
            rows.append((E, slope, intercept))
        params = {
            "num_points": len(rows),
        }
        return FitResult(model_name="Koutecky-Levich", parameters=params)


