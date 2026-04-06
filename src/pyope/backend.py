"""Compute backend configuration for pyope."""

from contextlib import contextmanager
from typing import Iterator


SUPPORTED_BACKENDS = {"sympy", "wolfram"}

_current_backend = "sympy"


def get_compute_backend() -> str:
    """Return the current compute backend name."""
    return _current_backend


def set_compute_backend(name: str) -> None:
    """Set the process-wide compute backend."""
    normalized = name.strip().lower()
    if normalized not in SUPPORTED_BACKENDS:
        supported = ", ".join(sorted(SUPPORTED_BACKENDS))
        raise ValueError(
            f"Unsupported compute backend '{name}'. Expected one of: {supported}"
        )

    global _current_backend
    _current_backend = normalized


@contextmanager
def compute_backend(name: str) -> Iterator[None]:
    """Temporarily switch the compute backend within a context."""
    previous = get_compute_backend()
    set_compute_backend(name)
    try:
        yield
    finally:
        set_compute_backend(previous)
