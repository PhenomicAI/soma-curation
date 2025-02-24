import abc
import concurrent.futures
from typing import Any, List, Callable, Tuple, TypeVar, Generic

from ..sc_logging import init_worker_logging

T = TypeVar("T")


class ExecutionResult(Generic[T]):
    """
    Holds the outcome of executing tasks.
    """

    def __init__(self) -> None:
        self.successes: List[T] = []
        self.failures: List[Tuple[Any, Exception]] = []

    @property
    def all_successful(self) -> bool:
        return len(self.failures) == 0

    @property
    def num_failures(self) -> int:
        return len(self.failures)

    @property
    def num_successes(self) -> int:
        return len(self.successes)

    def __repr__(self) -> str:
        return f"ExecutionResult(successes={len(self.successes)}, failures={len(self.failures)})"


class ExecutorBase(abc.ABC):
    """Abstract base class for concurrency executors."""

    @abc.abstractmethod
    def run(self, tasks: List[Any], func: Callable[..., T]) -> ExecutionResult[T]:
        pass


class MultiprocessingExecutor(ExecutorBase):
    """
    Run tasks in parallel using ProcessPoolExecutor, capturing errors gracefully.
    Optionally takes an `initializer` function (plus `initargs`) that is invoked
    *once* in each worker process before it starts running tasks.
    """

    def __init__(
        self,
        processes: int = 4,
        init_worker_logging: Callable[[str, str], None] = init_worker_logging,
        init_args: Tuple[str] = (),
    ):
        self.processes = processes
        self.init_worker_logging = init_worker_logging
        self.init_args = init_args

    def run(self, tasks: Tuple[Any], func: Callable[..., T]) -> ExecutionResult[T]:
        result = ExecutionResult[T]()
        if isinstance(tasks[0], str):
            raise ValueError("Tasks need to be tuples!")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.processes, initializer=self.init_worker_logging, initargs=self.init_args
        ) as executor:
            future_to_task = {executor.submit(func, *task): task for task in tasks}
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    output = future.result()
                    result.successes.append(output)
                except Exception as exc:
                    result.failures.append((task, exc))

        return result
