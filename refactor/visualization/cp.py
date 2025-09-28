from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt

from ..models import CPData
from .base import Plotter


class CPPlotter(Plotter):
    def create(self, traces: Iterable[CPData]) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        for trace in traces:
            if trace.time_s is None:
                continue
            ax.plot(trace.time_s, trace.potential_V, lw=1.0, label=trace.metadata.get("source", "CP"))
        ax.set_xlabel("t (s)")
        ax.set_ylabel("E (V vs reference)")
        ax.legend(loc="best")
        ax.grid(True, ls=":", lw=0.5)
        return fig


