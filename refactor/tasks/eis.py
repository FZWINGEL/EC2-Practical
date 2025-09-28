from __future__ import annotations

import pandas as pd

from ..services.task32_eis import run_eis_analysis


def build_eis_task(config: AnalysisConfig):
    def _run() -> TaskReport:
        return run_eis_analysis(config)

    return _run


