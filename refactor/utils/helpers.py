from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

import pandas as pd


def get_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first matching column name from the candidate list."""
    for column in candidates:
        if column in df.columns:
            return column
    return None


def label_from_filename(path: Path) -> str:
    """Return a human-readable label derived from a file stem."""
    return path.stem


def extract_rpm_from_name(path: Path | str) -> Optional[int]:
    """Extract the rotation rate (rpm) from a filename, if present."""
    name = path if isinstance(path, str) else path.name
    match = re.search(r"(\d+)\s*rpm", name, flags=re.IGNORECASE)
    if not match:
        match = re.search(r"(\d+)rpm", name, flags=re.IGNORECASE)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return None
    return None


