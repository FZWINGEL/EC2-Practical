from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from ..models import ElectrochemicalMeasurement, CVData, CPData, EISData, PEISData, FitResult
from .base import DataRepository


class BiologicRepository(DataRepository):
    """Repository for handling BioLogic instrument exports."""

    def load(self, path: Path) -> ElectrochemicalMeasurement:
        df = self._read_table(path)
        if "freq/Hz" in df.columns:
            return self._load_eis(df, path)
        if "mode" in df.columns and "chrono" in "".join(df["mode"].astype(str)).lower():
            return self._load_cp(df, path)
        if "cycle" in "".join(df.columns).lower():
            return self._load_cv(df, path)
        raise ValueError(f"Unrecognized BioLogic file format: {path.name}")

    def save_results(self, results: Iterable[FitResult], destination: Path) -> None:
        rows = []
        for result in results:
            row = {"model": result.model_name}
            row.update(result.parameters)
            if result.goodness_of_fit is not None:
                row["chi2"] = result.goodness_of_fit
            rows.append(row)
        df = pd.DataFrame(rows)
        destination.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(destination, index=False)

    def _read_table(self, path: Path) -> pd.DataFrame:
        try:
            df = pd.read_csv(path, sep="\t", comment="#", encoding="utf-8")
        except UnicodeDecodeError:
            df = pd.read_csv(path, sep="\t", comment="#", encoding="latin-1")
        df.columns = [str(col).strip() for col in df.columns]
        return df

    def _load_cv(self, df: pd.DataFrame, path: Path) -> CVData:
        e_v = pd.to_numeric(df.get("Ewe/V"), errors="coerce").to_numpy()
        i_mA = pd.to_numeric(df.get("<I>/mA"), errors="coerce").to_numpy()
        current_A = i_mA / 1000.0
        t_s = pd.to_numeric(df.get("time/s"), errors="coerce").to_numpy() if "time/s" in df else None
        return CVData(
            potential_V=e_v,
            current_A=current_A,
            time_s=t_s,
            metadata={"source": path.name},
        )

    def _load_cp(self, df: pd.DataFrame, path: Path) -> CPData:
        e_v = pd.to_numeric(df.get("Ewe/V"), errors="coerce").to_numpy()
        i_mA = pd.to_numeric(df.get("<I>/mA"), errors="coerce").to_numpy()
        current_A = i_mA / 1000.0
        t_s = pd.to_numeric(df.get("time/s"), errors="coerce").to_numpy() if "time/s" in df else None
        mode = df.get("mode")
        applied_current = float(pd.to_numeric(df.get("I Range"), errors="coerce").median()) if mode is not None else None
        return CPData(
            potential_V=e_v,
            current_A=current_A,
            time_s=t_s,
            applied_current_A=applied_current,
            metadata={"source": path.name},
        )

    def _load_eis(self, df: pd.DataFrame, path: Path) -> EISData:
        freq = pd.to_numeric(df.get("freq/Hz"), errors="coerce").to_numpy()
        re = pd.to_numeric(df.get("Re(Z)/Ohm"), errors="coerce").to_numpy()
        im_raw = pd.to_numeric(df.get("-Im(Z)/Ohm"), errors="coerce").to_numpy()
        im = -im_raw
        metadata = {"source": path.name}
        rpm = self._extract_rpm(path.name)
        if rpm is not None:
            return PEISData(
                frequency_Hz=freq,
                impedance_real_Ohm=re,
                impedance_imag_Ohm=im,
                rotation_speed_rpm=rpm,
                metadata=metadata,
            )
        return EISData(frequency_Hz=freq, impedance_real_Ohm=re, impedance_imag_Ohm=im, metadata=metadata)

    def _extract_rpm(self, name: str) -> int | None:
        import re

        match = re.search(r"(\d+)\s*rpm", name, flags=re.IGNORECASE)
        if not match:
            match = re.search(r"(\d+)rpm", name, flags=re.IGNORECASE)
        if match:
            try:
                return int(match.group(1))
            except ValueError:
                return None
        return None


