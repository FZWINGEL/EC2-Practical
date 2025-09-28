from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np


@dataclass(slots=True)
class ElectrochemicalMeasurement:
    potential_V: np.ndarray
    current_A: np.ndarray
    metadata: Dict[str, object] = field(default_factory=dict)
    time_s: Optional[np.ndarray] = None


@dataclass(slots=True)
class CVData(ElectrochemicalMeasurement):
    scan_rate_Vs: Optional[float] = None
    cycle_number: Optional[int] = None


@dataclass(slots=True)
class CPData(ElectrochemicalMeasurement):
    applied_current_A: Optional[float] = None


@dataclass(slots=True)
class EISData:
    frequency_Hz: np.ndarray
    impedance_real_Ohm: np.ndarray
    impedance_imag_Ohm: np.ndarray
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass(slots=True)
class PEISData(EISData):
    rotation_speed_rpm: Optional[int] = None


@dataclass(slots=True)
class FitResult:
    model_name: str
    parameters: Dict[str, float]
    parameter_errors: Optional[Dict[str, float]] = None
    goodness_of_fit: Optional[float] = None


