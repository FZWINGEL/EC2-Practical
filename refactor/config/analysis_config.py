from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional


@dataclass(slots=True)
class AnalysisConfig:
    """Central configuration object for the refactored analysis pipeline."""

    root: Path
    data_dir: Path = field(init=False)
    figures_dir: Path = field(init=False)
    results_dir: Path = field(init=False)
    parameters: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.root = self.root.resolve()
        self.data_dir = self.root / "Data"
        self.figures_dir = self.root / "figures"
        self.results_dir = self.root / "results"

    @classmethod
    def from_root(cls, root: Path, **overrides: Any) -> "AnalysisConfig":
        params = dict(overrides)
        return cls(root=root, parameters=params)

    def with_overrides(self, **overrides: Any) -> "AnalysisConfig":
        merged = dict(self.parameters)
        merged.update(overrides)
        return AnalysisConfig(root=self.root, parameters=merged)

    def get(self, key: str, default: Optional[Any] = None) -> Any:
        return self.parameters.get(key, default)

    def ensure_output_dirs(self) -> None:
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)

    def resolve_data(self, relative: Path | str) -> Path:
        return (self.data_dir / relative).resolve()

    def resolve_results(self, relative: Path | str) -> Path:
        return (self.results_dir / relative).resolve()

    def resolve_figures(self, relative: Path | str) -> Path:
        return (self.figures_dir / relative).resolve()

    def iter_task_parameters(self, prefix: str) -> Iterable[tuple[str, Any]]:
        for key, value in self.parameters.items():
            if key.startswith(prefix):
                yield key, value

    def update_from_mapping(self, mapping: Mapping[str, Any]) -> None:
        self.parameters.update(mapping)


