"""
PyOPE - Python Operator Product Expansion Library

A Python library for computing Operator Product Expansions (OPE)
in Vertex Operator Algebras (VOA).

基于 Mathematica 包 OPEdefs 的 Python 实现。
"""

__version__ = "0.1.0"
__author__ = "PyOPE Contributors"

# 已实现的模块
from .api import NO, NO_product, OPE, MakeOPE, bracket, normal_product

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

from .ope_data import OPEData
from .operator_spaces import (
    C2NullSearcher,
    C2Space,
    DescendantSpace,
    LocalOperatorBasis,
    make_realized,
    RealizedGenerator,
    SingularVectorAnalyzer,
    independent_under_realization,
    realize,
    realize_and_simplify,
    realized_coordinates,
)
from .operators import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    Operator,
    d,
    dn,
)

# Registry 和 API 模块
from .registry import Bosonic, Fermionic, OPERegistry, clear_registry, ope_registry

# Simplification 模块
from .simplify import (
    collect_normal_ordered_terms,
    combine_normal_ordered_terms,
    expand_nested_no,
    simplify,
)

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
    "C2NullSearcher",
    "C2Space",
    "DescendantSpace",
    "LocalOperatorBasis",
    "make_realized",
    "RealizedGenerator",
    "SingularVectorAnalyzer",
    "realize",
    "realize_and_simplify",
    "realized_coordinates",
    "independent_under_realization",
    # Registry
    "OPERegistry",
    "ope_registry",
    "clear_registry",
    "Bosonic",
    "Fermionic",
    # API
    "OPE",
    "NO",
    "NO_product",
    "bracket",
    "MakeOPE",
    "normal_product",
    # Simplification
    "simplify",
    "collect_normal_ordered_terms",
    "combine_normal_ordered_terms",
    "expand_nested_no",
    # Jacobi identity
    "check_jacobi_identity",
    "verify_jacobi_identity",
    # Cache
    "get_ope_cache",
]
