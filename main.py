from modules.config import (
    GC_AREA_CM2,
    GC_DIAMETER_MM,
    ROOT,
    hydrodynamics_summary,
)
from modules.plotting import setup_plot_style
from modules.task32_cp import analyze_task32_cp
from modules.task32_cv import analyze_task32_cvs
from modules.task32_eis import analyze_task32_eis
from modules.task33 import analyze_task33_lsv_and_eis


def main() -> None:
    setup_plot_style()
    print(f"Working directory: {ROOT}")
    print(f"Electrode area (GC, d={GC_DIAMETER_MM} mm): {GC_AREA_CM2:.4f} cm^2")
    print(hydrodynamics_summary())

    analyze_task32_cvs()
    analyze_task32_cp()
    analyze_task32_eis()
    analyze_task33_lsv_and_eis()


if __name__ == "__main__":
    main()

