from __future__ import annotations

from pathlib import Path
from .config import AnalysisConfig
from .pipeline import AnalysisApp, AnalysisPipeline
from .tasks import build_cp_task, build_cv_task, build_eis_task, build_task33


def build_analysis_app(root: Path) -> AnalysisApp:
    config = AnalysisConfig.from_root(root)
    pipeline = AnalysisPipeline()
    pipeline.add_task("Task 3.2 - CV", build_cv_task(config))
    pipeline.add_task("Task 3.2 - CP", build_cp_task(config))
    pipeline.add_task("Task 3.2 - EIS", build_eis_task(config))
    pipeline.add_task("Task 3.3 - LSV/EIS", build_task33(config))
    return AnalysisApp(config=config, pipeline=pipeline)


