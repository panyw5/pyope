"""
PyOPE - Python Operator Product Expansion Library

A Python library for computing Operator Product Expansions (OPE)
in Vertex Operator Algebras (VOA).

基于 Mathematica 包 OPEdefs 的 Python 实现。
"""

__version__ = "0.1.0"
__author__ = "PyOPE Contributors"

# 已实现的模块
from .api import NO, OPE, MakeOPE, bracket, normal_product

# 缓存模块
from .cache import get_ope_cache
from .constants import (
    ConstantOperator,
    Delta,
    One,
    Zero,
)

# Jacobi 恒等式模块
from .jacobi import check_jacobi_identity, verify_jacobi_identity
from .local_operator import (
    LocalOperator,
    OperatorProduct,
    OperatorSum,
    collect_operator_terms,
    extract_scalar_operator,
    get_operator_parity,
    is_local_operator,
    simplify_operator_expr,
)

# Null states 计算模块
from .null_states import (
    CoefficientExtractor,
    CoefficientMatrixBuilder,
    FockSpaceBasis,
    GroupedNullStatesCalculator,
    NullStatesCalculator,
    OperatorEnumerator,
    OperatorExpander,
    QuantumNumberCalculator,
    QuantumNumberGrouper,
    calculate_null_states,
    enumerate_fock_basis,
    extract_coefficients,
)
from .ope_data import OPEData
from .operators import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    Operator,
    d,
    dn,
)

# Registry 和 API 模块
from .registry import Bosonic, Fermionic, OPERegistry, ope_registry

# Simplification 模块
from .simplify import collect_normal_ordered_terms, expand_nested_no, simplify

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
    "normal_product",
    # Simplification
    "simplify",
    "collect_normal_ordered_terms",
    "expand_nested_no",
    # Jacobi identity
    "check_jacobi_identity",
    "verify_jacobi_identity",
    # Cache
    "get_ope_cache",
    # Null states
    "CoefficientExtractor",
    "FockSpaceBasis",
    "OperatorExpander",
    "CoefficientMatrixBuilder",
    "NullStatesCalculator",
    "QuantumNumberCalculator",
    "QuantumNumberGrouper",
    "GroupedNullStatesCalculator",
    "OperatorEnumerator",
    "extract_coefficients",
    "enumerate_fock_basis",
    "calculate_null_states",
]
