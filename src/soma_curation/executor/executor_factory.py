# concurrency_config.py
from pydantic import BaseModel, Field
from typing import Literal

from .executors import ExecutorBase, SerialExecutor, MultiprocessingExecutor


class ExecutorFactory(BaseModel):
    """Factory method for selecting an executor"""

    mode: Literal["serial", "multiprocessing", "tiledb-cloud", "beam"] = "serial"
    processes: int = Field(default=2, description="Number of processes (if multiprocess).")

    def create_executor(self) -> ExecutorBase:
        """Factory method to build an ExecutorBase instance from the config."""
        if self.mode == "serial":
            return SerialExecutor()
        elif self.mode == "multiprocessing":
            return MultiprocessingExecutor(processes=self.processes)
        elif self.mode == "tiledb-cloud":
            raise NotImplementedError("We will work to integrate tiledb-cloud in the future.")
        elif self.mode == "beam":
            raise NotImplementedError("We will work to integrate Apache Beam in the future.")
        else:
            raise ValueError(f"Unsupported mode: {self.mode}")
