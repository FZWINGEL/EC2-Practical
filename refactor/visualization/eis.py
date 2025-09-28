from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt

from ..models import EISData
from .base import Plotter


class EISPlotter(Plotter):
    def create(self, spectra: Iterable[EISData]) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.0, 4.2))
        for spectrum in spectra:
            ax.plot(spectrum.impedance_real_Ohm, -spectrum.impedance_imag_Ohm, "o-", lw=1.0, label=spectrum.metadata.get("source", "EIS"))
        ax.set_xlabel("Z' (Ω)")
        ax.set_ylabel("-Z'' (Ω)")
        ax.legend(loc="best")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, ls=":", lw=0.5)
        return fig


