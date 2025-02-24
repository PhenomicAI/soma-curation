import logging

from pathlib import Path

# Create a logger
logger = logging.getLogger("soma_curation")


def configure_logging(
    level: int = logging.WARNING, add_file_handler: bool = False, log_dir: str = None, log_file: str = None
) -> None:
    """
    Configure logging with the specified level.

    This function sets the logging level and ensures that the logger has a handler.

    Args:
    - level (int): Logging level (e.g., logging.DEBUG, logging.INFO)
    """
    _set_level(level)


def warning() -> None:
    """
    Set logging to a WARNING level.

    This function sets the logging level to WARNING.
    """
    _set_level(logging.WARNING)


def info() -> None:
    """
    Set logging to an INFO level.

    This function sets the logging level to INFO.
    """
    _set_level(logging.INFO)


def debug() -> None:
    """
    Set logging to a DEBUG level.

    This function sets the logging level to DEBUG.
    """
    _set_level(logging.DEBUG)


def error() -> None:
    """
    Set logging to an ERROR level.

    This function sets the logging level to ERROR.
    """
    _set_level(logging.ERROR)


def critical() -> None:
    """
    Set logging to a CRITICAL level.

    This function sets the logging level to CRITICAL.
    """
    _set_level(logging.CRITICAL)


def _set_level(level: int, add_file_handler: bool = False, log_dir: Path = None, log_file: str = None) -> None:
    """
    Set the logging level for the logger.

    This private function sets the logging level and adds a StreamHandler with a specific formatter if no handlers are present.

    Args:
    - level (int): Logging level (e.g., logging.DEBUG, logging.INFO)
    """
    logger.setLevel(level)
    if add_file_handler:
        if log_dir is None or log_file is None:
            raise ValueError("Specify a log directory and filename if you want to store logs in a file!")
        else:
            log_dir.mkdir(parents=True, exist_ok=True)
            logs_file = log_dir / log_file
            handler = logging.FileHandler(logs_file)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
            logger.addHandler(handler)

    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)


def init_worker_logging(log_level: int, log_dir: str, log_file: str):
    """
    Called once per worker process. Configures file logging with the same
    log level, directory, and filename as the main process.
    """
    logger.setLevel(log_level)

    # If you want a separate handler per worker, that's possible.
    # They will all write to the same file, leading to interleaved logs.
    log_path = Path(log_dir) / log_file
    fh = logging.FileHandler(log_path, mode="a")
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - [pid %(process)d] %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)


# Initialize logging with a default level (optional)
configure_logging(logging.WARNING)
