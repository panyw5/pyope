"""Helpers for locating the Wolfram command-line runtime."""

from __future__ import annotations

import shutil


def get_missing_wolframscript_message(context: str) -> str:
    """Return a user-facing installation hint for missing wolframscript."""
    return (
        f"{context} requires `wolframscript`, but it was not found in PATH. "
        "Install Wolfram Engine or Mathematica so the `wolframscript` command is available. "
        "See https://www.wolfram.com/wolframscript/ for installation instructions. "
        "Then verify the installation with `wolframscript -version`. "
        "On macOS, if you installed a Wolfram app bundle, you may need to add its `Contents/MacOS` directory to PATH. "
        'If you want to keep working without it, switch back to `set_compute_backend("sympy")`.'
    )


def require_wolframscript(context: str) -> str:
    """Return the wolframscript executable path or raise a clear error."""
    executable = shutil.which("wolframscript")
    if executable is None:
        raise FileNotFoundError(get_missing_wolframscript_message(context))
    return executable
