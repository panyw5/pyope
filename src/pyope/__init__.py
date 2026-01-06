"""
PyOPE - Python Operator Product Expansion Library

A Python library for computing Operator Product Expansions (OPE)
in Vertex Operator Algebras (VOA).

基于 Mathematica 包 OPEdefs 的 Python 实现。
"""

__version__ = "0.1.0"
__author__ = "PyOPE Contributors"

# 已实现的模块
from .operators import (
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    d,
    dn,
)
from .local_operator import (
    LocalOperator,
    OperatorSum,
    OperatorProduct,
    is_local_operator,
    extract_scalar_operator,
    get_operator_parity,
    simplify_operator_expr,
    collect_operator_terms,
)
from .constants import (
    ConstantOperator,
    One,
    Zero,
    Delta,
)
from .ope_data import OPEData

# Registry 和 API 模块
from .registry import OPERegistry, ope_registry, Bosonic, Fermionic
from .api import OPE, NO, bracket, MakeOPE

__all__ = [
    # Version info
    "__version__",
    "__author__",
    # Operators
    "Operator",
    "BasisOperator",
    "DerivativeOperator",
    "NormalOrderedOperator",
    "d",
    "dn",
    # Local operators
    "LocalOperator",
    "OperatorSum",
    "OperatorProduct",
    "is_local_operator",
    "extract_scalar_operator",
    "get_operator_parity",
    "simplify_operator_expr",
    "collect_operator_terms",
    # Constants
    "ConstantOperator",
    "One",
    "Zero",
    "Delta",
    # OPE Data
    "OPEData",
    # Registry
    "OPERegistry",
    "ope_registry",
    "Bosonic",
    "Fermionic",
    # API
    "OPE",
    "NO",
    "bracket",
    "MakeOPE",
]
