from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from ..config import AnalysisConfig
from .core import AnalysisPipeline
from .reporting import TaskReport


@dataclass
class AnalysisApp:
    config: AnalysisConfig
    pipeline: AnalysisPipeline

    def run(self) -> Iterable[TaskReport]:
        self.config.ensure_output_dirs()
        return self.pipeline.run()


