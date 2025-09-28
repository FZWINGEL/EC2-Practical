from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable

from ..models import ElectrochemicalMeasurement, FitResult


class DataRepository(ABC):
    """Abstract base class for loading and storing electrochemical data."""

    def __init__(self, root: Path) -> None:
        self.root = root.resolve()

    @abstractmethod
    def load(self, path: Path) -> ElectrochemicalMeasurement:
        """Load a measurement from disk."""

    @abstractmethod
    def save_results(self, results: Iterable[FitResult], destination: Path) -> None:
        """Persist analysis results to disk."""


