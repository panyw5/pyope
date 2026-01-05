"""
PyOPE - Python Operator Product Expansion Library

A high-performance Python library for computing Operator Product Expansions (OPE)
in Vertex Operator Algebras (VOA).

Based on the Mathematica package OPEdefs by Kris Thielemans.
"""

__version__ = "0.1.0"
__author__ = "PyOPE Contributors"

# Local operator base classes
from .local_operator import (
    LocalOperator,
    OperatorSum,
    OperatorProduct,
    is_local_operator,
    extract_scalar_operator,
)

# Core operator classes
from .operators import (
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
)

from .ope_data import OPEData
from .registry import OPERegistry, ope_registry

# Public API functions
from .api import OPE, NO, bracket
from .operators import d, dn

# Special operators
from .constants import One, Delta

__all__ = [
    # Version
    "__version__",

    # Local operator base
    "LocalOperator",
    "OperatorSum",
    "OperatorProduct",
    "is_local_operator",
    "extract_scalar_operator",

    # Operator classes
    "Operator",
    "BasisOperator",
    "DerivativeOperator",
    "NormalOrderedOperator",

    # OPE data structures
    "OPEData",
    "OPERegistry",
    "ope_registry",

    # API functions
    "OPE",
    "NO",
    "d",
    "dn",
    "bracket",

    # Constants
    "One",
    "Delta",
]
