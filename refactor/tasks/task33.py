from __future__ import annotations

from ..services.task33 import run_task33_lsv_eis


def build_task33(config: AnalysisConfig):
    def _run() -> TaskReport:
        return run_task33_lsv_eis(config)

    return _run


