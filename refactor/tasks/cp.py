from __future__ import annotations

import pandas as pd

from ..config import AnalysisConfig
from ..models import CPData
from ..pipeline import TaskReport
from ..repositories import BiologicRepository
from ..services.task32_cp import run_cp_analysis


def build_cp_task(config: AnalysisConfig):
    def _run() -> TaskReport:
        return run_cp_analysis(config)

    return _run


