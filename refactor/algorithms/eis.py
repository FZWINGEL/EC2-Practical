from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

import numpy as np

from ..config import AnalysisConfig
from ..models import EISData, FitResult


@dataclass(slots=True)
class CircuitModel:
    name: str
    parameter_names: List[str]


@dataclass(slots=True)
class EISCircuitFitter:
    config: AnalysisConfig

    MODELS = {
        "randles": CircuitModel(name="R0-p(R1,C1)-W1", parameter_names=["Rs", "Rct", "Cdl", "sigma"]),
    }

    def fit(self, data: EISData, model: str = "randles") -> FitResult:
        circuit = self.MODELS.get(model)
        if circuit is None:
            raise ValueError(f"Unknown circuit model: {model}")

        Rs = float(np.nanmedian(data.impedance_real_Ohm[:5]))
        Rct = float(np.nanmax(data.impedance_real_Ohm) - Rs)
        Cdl = float(1.0 / (2.0 * np.pi * data.frequency_Hz[np.nanargmax(-data.impedance_imag_Ohm)] * max(Rct, 1e-9)))
        params = {"Rs": Rs, "Rct": Rct, "Cdl": Cdl, "sigma": float("nan")}
        return FitResult(model_name=circuit.name, parameters=params, goodness_of_fit=None)


