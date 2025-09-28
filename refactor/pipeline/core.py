from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, List, Optional, Sequence

from .reporting import TaskReport


TaskCallable = Callable[[], Optional[TaskReport]]


@dataclass(slots=True)
class PipelineTask:
    name: str
    callable: TaskCallable

    def run(self) -> TaskReport:
        report = self.callable()
        if report is None:
            return TaskReport(name=self.name)
        if not report.name:
            report.name = self.name
        return report


@dataclass
class AnalysisPipeline:
    tasks: List[PipelineTask] = field(default_factory=list)

    def add_task(self, name: str, callable: TaskCallable) -> None:
        self.tasks.append(PipelineTask(name=name, callable=callable))

    def extend(self, tasks: Sequence[PipelineTask]) -> None:
        self.tasks.extend(tasks)

    def run(self) -> List[TaskReport]:
        reports: List[TaskReport] = []
        for task in self.tasks:
            reports.append(task.run())
        return reports


