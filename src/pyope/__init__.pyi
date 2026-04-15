from .api import MakeOPE, NO, NO_product, OPE, bracket
from .cache import get_ope_cache
from .constants import ConstantOperator, Delta, One, Zero
from .jacobi import check_jacobi_identity, verify_jacobi_identity
from .local_operator import (
    LocalOperator,
    OperatorProduct,
    OperatorSum,
    collect_operators_coefficients,
    extract_scalar_operator,
    get_operator_parity,
    is_local_operator,
    simplify_with_sympy,
)
from .ope_data import OPEData
from .compact_ope import compact_family_poles
from .quasiprimary import qp, quasiprimary_product
from .c2 import AbstractC2Reducer, C2ReductionWitness, GenericC2Reducer
from .descendants import DescendantSpace
from .free_field_c2 import DerivativeKillingFreeFieldC2Reducer, FreeFieldC2Reducer
from .finite_algebra import (
    AlgebraElement,
    AlgebraValidationIssue,
    AssociativityValidationResult,
    build_finite_dimensional_algebra,
    FiniteDimensionalAlgebra,
    IdentityValidationResult,
    IrreducibleRepresentationClassification,
    MultiplicationTableEntry,
    OneDimensionalRepresentation,
)
from .null_search import C2NullSearcher as QuotientC2NullSearcher, NullSearchResult
from .realizations import (
    DerivativeKillingRealizationBackend,
    IdentityRealizationBackend,
    RealizationBackend,
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
from .operator_spaces import (
    C2NullSearcher,
    C2Space,
    LocalOperatorBasis,
    LocalOperatorCanonicalizer,
    SparseLinearContext,
    list_independent_op_indices,
    list_independent_ops,
    list_zero_relations,
    make_realized,
    RealizedGenerator,
    realize,
    realized_coordinates,
    clear_realize_cache,
)
from .operators import (
    BasicOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    Operator,
    d,
    dn,
)
from .registry import (
    Bosonic,
    Fermionic,
    OPERegistry,
    clear_registry,
    ope_registry,
)
from .simplify import (
    collect_normal_ordered_terms,
    combine_normal_ordered_terms,
    expand_nested_no,
    simplify,
)
from .wolfram_backend import (
    op_to_wolfram_string,
    simplify_with_wolfram,
)

__version__: str
__author__: str
