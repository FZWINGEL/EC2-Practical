from pathlib import Path
from typing import Optional

import pandas as pd

from .config import RES_DIR, ROOT, ensure_dirs


HEADER_MARKERS = ("freq/Hz", "time/s", "mode")


def _detect_header_start(path: Path) -> int:
    """Return the zero-based line index where the tabular header starts."""
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle):
            first_token = line.strip().split("\t", 1)[0]
            if first_token in HEADER_MARKERS:
                return idx
    return 0


def read_biologic_table(path: Path) -> pd.DataFrame:
    """Read a BioLogic export (tab-delimited) into a DataFrame."""
    header_idx = _detect_header_start(path)
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            engine="python",
            skiprows=header_idx,
            comment="#",
            encoding="utf-8",
        )
    except UnicodeDecodeError:
        df = pd.read_csv(
            path,
            sep="\t",
            engine="python",
            skiprows=header_idx,
            comment="#",
            encoding="latin-1",
        )
    df.columns = [str(col).strip() for col in df.columns]
    return df


def write_csv(df: pd.DataFrame, filename: str) -> None:
    """Write a DataFrame into the results directory."""
    ensure_dirs()
    out_path = RES_DIR / filename
    df.to_csv(out_path, index=False)
    print(f"Saved: {(out_path).relative_to(ROOT)}")

