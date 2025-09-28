from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .base import Plotter


class TafelPlotter(Plotter):
    def create(self, eta_V: np.ndarray, jk_A_cm2: np.ndarray) -> plt.Figure:
        mask = np.isfinite(eta_V) & np.isfinite(jk_A_cm2) & (np.abs(jk_A_cm2) > 1e-12)
        fig, ax = plt.subplots(figsize=(5.4, 4.0))
        ax.plot(np.log10(np.abs(jk_A_cm2[mask])), eta_V[mask], "o")
        ax.set_xlabel("log10(|j_k| / A cm^-2)")
        ax.set_ylabel("η (V)")
        ax.grid(True, ls=":", lw=0.5)
        return fig


