from __future__ import annotations

import pandas as pd

from ..services.task32_cv import run_cv_analysis
from ..config import AnalysisConfig
from ..models import CVData
from ..pipeline import TaskReport
from ..repositories import BiologicRepository
from ..services import CVProcessor
from ..visualization import CVPlotter


def build_cv_task(config: AnalysisConfig):
    def _run() -> TaskReport:
        return run_cv_analysis(config)

    return _run


