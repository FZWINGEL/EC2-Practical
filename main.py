from modules.config import (
    GC_AREA_CM2,
    GC_DIAMETER_MM,
    ROOT,
    ensure_dirs,
    hydrodynamics_summary,
)
from modules.pipeline import AnalysisPipeline
from modules.plotting import setup_plot_style
from modules.task32_cp import analyze_task32_cp
from modules.task32_cv import analyze_task32_cvs
from modules.task32_eis import analyze_task32_eis
from modules.task33 import analyze_task33_lsv_and_eis


def main() -> None:
    setup_plot_style()
    ensure_dirs()
    print(f"Working directory: {ROOT}")
    print(f"Electrode area (GC, d={GC_DIAMETER_MM} mm): {GC_AREA_CM2:.4f} cm^2")
    print(hydrodynamics_summary())

    pipeline = AnalysisPipeline()
    pipeline.add_task("Task 3.2 - CV", analyze_task32_cvs)
    pipeline.add_task("Task 3.2 - CP", analyze_task32_cp)
    pipeline.add_task("Task 3.2 - EIS", analyze_task32_eis)
    pipeline.add_task("Task 3.3 - LSV/EIS", analyze_task33_lsv_and_eis)

    reports = pipeline.run()
    AnalysisPipeline.summarize(reports)


if __name__ == "__main__":
    main()
