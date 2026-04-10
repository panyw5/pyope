"""Compute backend configuration for pyope.

The Wolfram backend is currently incomplete and intentionally blocked from
public selection.
"""

from contextlib import contextmanager
from typing import Iterator


SUPPORTED_BACKENDS = {"sympy", "wolfram"}
DISABLED_BACKENDS = {"wolfram"}

_current_backend = "sympy"


def get_compute_backend() -> str:
    """Return the current compute backend name."""
    return _current_backend


def set_compute_backend(name: str) -> None:
    """Set the process-wide compute backend.

    Note:
        The ``wolfram`` backend is not ready for public use yet and is rejected
        here even though it remains listed as a known backend name for internal
        compatibility.
    """
    normalized = name.strip().lower()
    if normalized not in SUPPORTED_BACKENDS:
        supported = ", ".join(sorted(SUPPORTED_BACKENDS))
        raise ValueError(
            f"Unsupported compute backend '{name}'. Expected one of: {supported}"
        )
    if normalized in DISABLED_BACKENDS:
        raise ValueError(
            f"Compute backend '{name}' is not ready yet and must not be used. "
            "Please use 'sympy' instead."
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
