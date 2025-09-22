from typing import List

import numpy as np


_GLOBAL_D_VALUES: List[float] = []


def add_diffusion_coefficient(value: float, source: str = "") -> None:
    """Store a positive, finite diffusion coefficient computed in this run."""
    try:
        numeric = float(value)
    except Exception:
        return
    if np.isfinite(numeric) and numeric > 0:
        _GLOBAL_D_VALUES.append(numeric)


def get_calculated_diffusion() -> float:
    """Return the median diffusion coefficient accumulated so far, or NaN."""
    if not _GLOBAL_D_VALUES:
        return float("nan")
    return float(np.nanmedian(_GLOBAL_D_VALUES))

