import matplotlib.pyplot as plt

from .config import FIG_DIR, ROOT, ensure_dirs


def setup_plot_style() -> None:
    """Set a consistent, clean plotting style across all figures."""
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.grid": True,
            "grid.linestyle": ":",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.3,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.2,
            "legend.frameon": True,
            "legend.fontsize": 9,
            "legend.framealpha": 0.85,
            "legend.edgecolor": "#00000055",
            "legend.title_fontsize": 9,
            "savefig.dpi": 300,
        }
    )


def beautify_axes(ax: plt.Axes) -> None:
    """Apply minor ticks, outward ticks, and simplify spines on an Axes."""
    try:
        ax.minorticks_on()
    except Exception:
        pass
    ax.tick_params(which="both", direction="out", length=5, width=0.8)
    ax.tick_params(which="minor", length=3, width=0.6)
    if "top" in ax.spines:
        ax.spines["top"].set_visible(False)
    if "right" in ax.spines:
        ax.spines["right"].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def safe_save(fig: plt.Figure, filename: str) -> None:
    """Save a figure into the figures directory with tight layout handling."""
    ensure_dirs()
    out_path = FIG_DIR / filename
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path.relative_to(ROOT)}")

