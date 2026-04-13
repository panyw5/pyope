"""Compute backend configuration for pyope."""

import os
from contextlib import contextmanager
from typing import Iterator

from .wolfram_dependency import get_missing_wolframscript_message, require_wolframscript


SUPPORTED_BACKENDS = {"sympy", "wolfram"}

_current_backend = "sympy"


def get_compute_backend() -> str:
    """Return the current compute backend name."""
    return _current_backend


def set_compute_backend(name: str, max_worker_number: int = 1) -> None:
    """Set the process-wide compute backend."""
    normalized = name.strip().lower()
    if normalized not in SUPPORTED_BACKENDS:
        supported = ", ".join(sorted(SUPPORTED_BACKENDS))
        raise ValueError(
            f"Unsupported compute backend '{name}'. Expected one of: {supported}"
        )
    if not isinstance(max_worker_number, int) or max_worker_number <= 0:
        raise ValueError("max_worker_number must be a positive integer")
    if normalized == "wolfram":
        try:
            require_wolframscript('`set_compute_backend("wolfram")`')
        except FileNotFoundError as exc:
            raise ValueError(
                get_missing_wolframscript_message('`set_compute_backend("wolfram")`')
            ) from exc

    global _current_backend
    _current_backend = normalized
    if normalized == "wolfram":
        os.environ["PYOPE_WL_MAX_WORKERS"] = str(max_worker_number)


@contextmanager
def compute_backend(name: str, max_worker_number: int = 1) -> Iterator[None]:
    """Temporarily switch the compute backend within a context."""
    previous = get_compute_backend()
    previous_max_workers = os.environ.get("PYOPE_WL_MAX_WORKERS")
    set_compute_backend(name, max_worker_number=max_worker_number)
    try:
        yield
    finally:
        if previous == "wolfram":
            restore_workers = previous_max_workers or "1"
            set_compute_backend(previous, max_worker_number=int(restore_workers))
        else:
            set_compute_backend(previous)
            if previous_max_workers is None:
                os.environ.pop("PYOPE_WL_MAX_WORKERS", None)
            else:
                os.environ["PYOPE_WL_MAX_WORKERS"] = previous_max_workers
