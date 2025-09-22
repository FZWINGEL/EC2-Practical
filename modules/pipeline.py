from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Sequence


TaskCallable = Callable[[], Optional["TaskReport"]]


@dataclass
class TaskReport:
    """Structured output for a single analysis task."""

    name: str
    figures: List[Path] = field(default_factory=list)
    tables: List[Path] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def record_figure(self, path: Optional[Path]) -> Optional[Path]:
        """Add a figure artifact to the report and return the path."""
        if path is None:
            return None
        resolved = Path(path)
        if resolved not in self.figures:
            self.figures.append(resolved)
        return resolved

    def record_table(self, path: Optional[Path]) -> Optional[Path]:
        """Add a tabular artifact to the report and return the path."""
        if path is None:
            return None
        resolved = Path(path)
        if resolved not in self.tables:
            self.tables.append(resolved)
        return resolved

    def add_message(self, message: str) -> None:
        if message:
            self.messages.append(message)

    def add_warning(self, message: str) -> None:
        if message:
            self.warnings.append(message)

    def add_error(self, message: str) -> None:
        if message:
            self.errors.append(message)


@dataclass
class AnalysisTask:
    """Container for an analysis callable and its display name."""

    name: str
    run: TaskCallable

    def execute(self) -> TaskReport:
        result = self.run()
        if result is None:
            return TaskReport(name=self.name)
        if not result.name:
            result.name = self.name
        return result


class AnalysisPipeline:
    """Coordinate execution of analysis tasks and collect their reports."""

    def __init__(self) -> None:
        self._tasks: List[AnalysisTask] = []

    def add_task(self, name: str, func: TaskCallable) -> None:
        self._tasks.append(AnalysisTask(name=name, run=func))

    def extend(self, tasks: Sequence[AnalysisTask]) -> None:
        self._tasks.extend(tasks)

    def run(self) -> List[TaskReport]:
        reports: List[TaskReport] = []
        for task in self._tasks:
            print(f"\n=== {task.name} ===")
            try:
                report = task.execute()
            except KeyboardInterrupt:  # pragma: no cover - allow graceful stop
                raise
            except Exception as exc:  # pragma: no cover - runtime guard
                print(f"[Pipeline] {task.name} failed: {exc}")
                report = TaskReport(name=task.name)
                report.add_error(str(exc))
                report.add_warning(traceback.format_exc())
            else:
                for message in report.messages:
                    print(f"[{task.name}] {message}")
                for warning in report.warnings:
                    print(f"[{task.name}] Warning: {warning}")
                for error in report.errors:
                    print(f"[{task.name}] Error: {error}")
            reports.append(report)
        return reports

    @staticmethod
    def summarize(reports: Sequence[TaskReport]) -> None:
        """Print a short summary report to stdout."""
        print("\n=== Pipeline summary ===")
        if not reports:
            print("No tasks executed.")
            return
        for report in reports:
            print(f"- {report.name}:")
            if report.errors:
                print(f"  errors: {len(report.errors)}")
            if report.warnings:
                print(f"  warnings: {len(report.warnings)}")
            if report.tables:
                print("  tables:")
                for path in report.tables:
                    print(f"    - {path}")
            if report.figures:
                print("  figures:")
                for path in report.figures:
                    print(f"    - {path}")
            if not any([report.errors, report.warnings, report.tables, report.figures]):
                print("  no recorded artifacts.")
