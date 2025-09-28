from __future__ import annotations

from typing import Dict

import matplotlib.pyplot as plt
import numpy as np

from .base import Plotter


class KouteckyLevichPlotter(Plotter):
    def create(self, j_by_rpm: Dict[int, np.ndarray], potential_index: int, potentials: np.ndarray) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(5.4, 4.0))
        rpms = np.array(sorted(j_by_rpm.keys()))
        omega = 2 * np.pi * (rpms / 60.0)
        inv_sqrt_omega = 1.0 / np.sqrt(omega)
        j_values = np.array([j_by_rpm[rpm][potential_index] for rpm in rpms])
        mask = np.isfinite(j_values) & (np.abs(j_values) > 1e-9)
        ax.plot(inv_sqrt_omega[mask], 1.0 / j_values[mask], "o", label=f"E={potentials[potential_index]:.3f} V")
        ax.set_xlabel("1/√ω (s^0.5)")
        ax.set_ylabel("1/j (cm^2 A^-1)")
        ax.grid(True, ls=":", lw=0.5)
        ax.legend(loc="best")
        return fig


