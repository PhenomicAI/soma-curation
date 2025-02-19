# executor.py
import abc
import concurrent.futures
from typing import Any, List, Callable, Tuple, TypeVar, Generic

# A generic type variable for our task results
T = TypeVar("T")


class ExecutionResult(Generic[T]):
    """
    Holds the outcome of executing tasks.

    Attributes:
        successes (List[T]): The list of task outputs that completed successfully.
        failures (List[Tuple[Any, Exception]]): A list of tuples, each containing the task input and the exception raised.
    """

    def __init__(self) -> None:
        self.successes: List[T] = []
        self.failures: List[Tuple[Any, Exception]] = []

    @property
    def all_successful(self) -> bool:
        """Return True if no task has failed."""
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
    def run(self, tasks: List[Any], func: Callable[[Any], T]) -> ExecutionResult[T]:
        """
        Execute `func` on each element of `tasks` and return an ExecutionResult that
        includes both successes and failures.
        """
        pass


class MultiprocessingExecutor(ExecutorBase):
    """Run tasks in parallel using ProcessPoolExecutor and capture errors gracefully."""

    def __init__(self, processes: int = 4):
        self.processes = processes

    def run(self, tasks: List[Any], func: Callable[[Any], T]) -> ExecutionResult[T]:
        result = ExecutionResult[T]()
        # Use ProcessPoolExecutor to submit tasks in parallel.
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.processes) as executor:
            future_to_task = {executor.submit(func, task): task for task in tasks}
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    output = future.result()
                    result.successes.append(output)
                except Exception as exc:
                    result.failures.append((task, exc))
        return result
