from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class TaskReport:
    name: str
    figures: List[Path] = field(default_factory=list)
    tables: List[Path] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def record_figure(self, path: Optional[Path]) -> None:
        if path is not None:
            self.figures.append(Path(path))

    def record_table(self, path: Optional[Path]) -> None:
        if path is not None:
            self.tables.append(Path(path))

    def add_message(self, message: str) -> None:
        if message:
            self.messages.append(message)

    def add_warning(self, message: str) -> None:
        if message:
            self.warnings.append(message)

    def add_error(self, message: str) -> None:
        if message:
            self.errors.append(message)


