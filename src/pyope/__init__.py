"""
PyOPE - Python Operator Product Expansion Library

A Python library for computing Operator Product Expansions (OPE)
in Vertex Operator Algebras (VOA).

基于 Mathematica 包 OPEdefs 的 Python 实现。
"""

__version__ = "0.1.0.post1"
__author__ = "PyOPE Contributors"

# 已实现的模块
from .api import NO, NO_product, OPE, MakeOPE, bracket, normal_product
from .backend import compute_backend, get_compute_backend, set_compute_backend
from .wolfram_backend import canonicalize_exprs as wolfram_canonicalize_exprs
from .wolfram_backend import evaluate_expr as evaluate_with_wolfram
from .wolfram_backend import evaluate_exprs as evaluate_many_with_wolfram
from .wolfram_backend import simplify_expr as simplify_with_wolfram
from .c2 import AbstractC2Reducer, C2ReductionWitness, GenericC2Reducer

# 缓存模块
from .cache import get_ope_cache
from .compact_ope import compact_family_poles
from .constants import (
    ConstantOperator,
    Delta,
    One,
    Zero,
)
from .descendants import DescendantSpace
from .finite_algebra import (
    AlgebraElement,
    AlgebraValidationIssue,
    AssociativityValidationResult,
    FiniteDimensionalAlgebra,
    IdentityValidationResult,
    IrreducibleRepresentationClassification,
    MultiplicationTableEntry,
    OneDimensionalRepresentation,
    build_finite_dimensional_algebra,
)
from .free_field_c2 import DerivativeKillingFreeFieldC2Reducer, FreeFieldC2Reducer

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
    simplify_with_sympy,
)
from .null_search import C2NullSearcher as QuotientC2NullSearcher
from .null_search import NullSearchResult
from .ope_data import OPEData
from .operator_spaces import (
    C2NullSearcher,
    C2Space,
    LocalOperatorBasis,
    LocalOperatorCanonicalizer,
    RealizedGenerator,
    SparseLinearContext,
    clear_realize_cache,
    independent_under_realization,
    list_independent_ops,
    list_zero_relations,
    make_realized,
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
from .quasiprimary import qp, quasiprimary_product
from .realizations import (
    DerivativeKillingRealizationBackend,
    IdentityRealizationBackend,
    RealizationBackend,
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
from .singularity import SingularVectorAnalyzer
from .zhu import (
    AbstractZhuReducer,
    GenericZhuReducer,
    ZhuReductionWitness,
    ZhuSpace,
    zhu_circle_product,
    zhu_star_product,
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
    "simplify_with_sympy",
    "collect_operator_terms",
    # Constants
    "ConstantOperator",
    "One",
    "Zero",
    "Delta",
    # OPE Data
    "OPEData",
    "compact_family_poles",
    "quasiprimary_product",
    "qp",
    "AbstractC2Reducer",
    "C2NullSearcher",
    "C2ReductionWitness",
    "C2Space",
    "DescendantSpace",
    "GenericC2Reducer",
    "AlgebraElement",
    "AlgebraValidationIssue",
    "AssociativityValidationResult",
    "build_finite_dimensional_algebra",
    "DerivativeKillingFreeFieldC2Reducer",
    "DerivativeKillingRealizationBackend",
    "FiniteDimensionalAlgebra",
    "FreeFieldC2Reducer",
    "IdentityRealizationBackend",
    "IdentityValidationResult",
    "IrreducibleRepresentationClassification",
    "LocalOperatorBasis",
    "LocalOperatorCanonicalizer",
    "MultiplicationTableEntry",
    "NullSearchResult",
    "OneDimensionalRepresentation",
    "QuotientC2NullSearcher",
    "RealizationBackend",
    "AbstractZhuReducer",
    "GenericZhuReducer",
    "SparseLinearContext",
    "ZhuReductionWitness",
    "ZhuSpace",
    "list_independent_ops",
    "list_zero_relations",
    "make_realized",
    "RealizedGenerator",
    "SingularVectorAnalyzer",
    "realize",
    "realize_and_simplify",
    "realized_coordinates",
    "clear_realize_cache",
    "independent_under_realization",
    "zhu_circle_product",
    "zhu_star_product",
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
    "get_compute_backend",
    "set_compute_backend",
    "compute_backend",
    "evaluate_many_with_wolfram",
    "evaluate_with_wolfram",
    "wolfram_canonicalize_exprs",
    "simplify_with_wolfram",
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
