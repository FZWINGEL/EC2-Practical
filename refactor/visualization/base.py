from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from ..config import AnalysisConfig


class Plotter(ABC):
    """Base class for all figure generators."""

    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config

    @abstractmethod
    def create(self, *args: Any, **kwargs: Any) -> plt.Figure:
        raise NotImplementedError

    def save(self, figure: plt.Figure, filename: str) -> Path:
        path = self.config.resolve_figures(filename)
        # Ensure filename has .png extension
        if not filename.endswith('.png'):
            filename = filename.rsplit('.', 1)[0] + '.png'
            path = self.config.resolve_figures(filename)
        figure.savefig(path, dpi=300, bbox_inches="tight", format='png', facecolor="white")
        plt.close(figure)
        return path


