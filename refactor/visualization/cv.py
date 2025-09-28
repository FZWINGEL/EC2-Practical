from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from ..models import CVData
from .base import Plotter


class CVPlotter(Plotter):
    def create(self, curves: Iterable[CVData]) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.0, 4.2))
        for curve in curves:
            ax.plot(curve.potential_V, curve.current_A, lw=1.0, label=curve.metadata.get("source", "CV"))
        ax.set_xlabel("E (V vs reference)")
        ax.set_ylabel("i (A)")
        ax.legend(loc="best")
        ax.grid(True, ls=":", lw=0.5)
        return fig


