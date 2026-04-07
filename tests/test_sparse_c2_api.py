import sympy as sp

from pyope import (
    AbstractC2Reducer,
    BasisOperator,
    Bosonic,
    DescendantSpace,
    Fermionic,
    GenericC2Reducer,
    LocalOperatorBasis,
    NullSearchResult,
    QuotientC2NullSearcher,
    SparseLinearContext,
    NO,
    d,
)


def test_local_operator_basis_exposes_sparse_terms_without_basis_enumeration():
    J = BasisOperator("J_sparse_terms", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])

    def fail_basis(*args, **kwargs):
        raise AssertionError("basis() should not be called for sparse_terms()")

    basis.basis = fail_basis  # type: ignore[method-assign]

    terms = basis.sparse_terms(3 * NO(J, J) + d(J))

    assert terms[NO(J, J)] == 3
    assert terms[d(J)] == 1


def test_sparse_terms_distributes_symbolic_scalars_over_operator_sums():
    b = BasisOperator("b_sparse_symbolic", fermionic=True, conformal_weight=2)
    c = BasisOperator("c_sparse_symbolic", fermionic=True, conformal_weight=-1)
    Fermionic(b, c)

    basis = LocalOperatorBasis([b, c], max_occurence=4)
    x = NO(b, c)
    y = NO(b, d(c))
    a0, a1 = sp.symbols("a0:2")

    terms = basis.sparse_terms(5 * x + y - a0 * (x + y) - a1 * (x - y))

    assert terms == {
        x: -a0 - a1 + 5,
        y: -a0 + a1 + 1,
    }


def test_sparse_linear_context_finds_independent_subset_and_relations():
    J = BasisOperator("J_sparse_context", conformal_weight=1)
    Bosonic(J)

    basis = LocalOperatorBasis([J])
    context = SparseLinearContext(basis)
    jj = NO(J, J)
    dj = d(J)

    independent = context.independent_subset([jj, dj, jj + dj])
    relations = context.zero_relations([jj, dj, jj + dj])

    assert len(independent) == 2
    assert len(relations) == 1
    assert basis.canonicalize(relations[0]["relation"]) == 0


def test_generic_c2_reducer_detects_c2_generator_and_obstruction():
    T = BasisOperator("T_generic_c2", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T])
    reducer = GenericC2Reducer(basis)

    assert isinstance(reducer, AbstractC2Reducer)
    assert reducer.is_zero_mod_c2(NO(d(T), T))
    assert reducer.quotient_normal_form(NO(T, T)) == NO(T, T)


def test_quotient_null_searcher_precheck_reports_status_and_remainder():
    T = BasisOperator("T_precheck", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], stress_tensor=T)
    reducer = GenericC2Reducer(basis)
    searcher = QuotientC2NullSearcher(canonicalizer=basis, c2_reducer=reducer)

    result_zero = searcher.quotient_precheck(NO(d(T), T))
    result_obstruction = searcher.search_stress_tensor_nilpotency(2)

    assert isinstance(result_zero, NullSearchResult)
    assert result_zero.status == "needs_lift"
    assert result_zero.quotient_remainder == 0
    assert result_obstruction.status == "obstructed"
    assert result_obstruction.quotient_remainder == NO(T, T)


def test_quotient_null_searcher_solves_descendant_lift_from_sources():
    T = BasisOperator("T_lift", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], stress_tensor=T)
    reducer = GenericC2Reducer(basis)
    searcher = QuotientC2NullSearcher(
        canonicalizer=basis,
        c2_reducer=reducer,
        descendants=DescendantSpace(basis),
    )

    result = searcher.search_from_sources(4, [T], NO(T, T))

    assert isinstance(result, NullSearchResult)
    assert result is not None
    assert result.status == "solved"
    assert result.null_operator == NO(T, T)
    assert result.c2_remainder == 0


def test_quotient_null_searcher_records_stress_tensor_search_metadata():
    T = BasisOperator("T_nilpotency_meta", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], stress_tensor=T)
    reducer = GenericC2Reducer(basis)
    searcher = QuotientC2NullSearcher(canonicalizer=basis, c2_reducer=reducer)

    result = searcher.search_stress_tensor_nilpotency(2)

    assert isinstance(result, NullSearchResult)
    assert result.n == 2
    assert result.target == NO(T, T)


def test_quotient_null_searcher_marks_failed_singularity_when_hook_rejects():
    class RejectAllConstraints:
        def positive_mode_constraints(self, expr):
            return {"T": {1: expr}}

    T = BasisOperator("T_singularity_reject", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], stress_tensor=T)
    reducer = GenericC2Reducer(basis)
    searcher = QuotientC2NullSearcher(
        canonicalizer=basis,
        c2_reducer=reducer,
        descendants=DescendantSpace(basis),
        singular_constraints=RejectAllConstraints(),
    )

    result = searcher.search_from_sources(4, [T], NO(T, T))

    assert result is not None
    assert result.status == "failed_singularity"
    assert result.is_singular is False


def test_null_search_result_legacy_payload_preserves_structured_fields():
    T = BasisOperator("T_legacy_payload", conformal_weight=2)
    Bosonic(T)

    basis = LocalOperatorBasis([T], stress_tensor=T)
    reducer = GenericC2Reducer(basis)
    searcher = QuotientC2NullSearcher(
        canonicalizer=basis,
        c2_reducer=reducer,
        descendants=DescendantSpace(basis),
    )

    result = searcher.search_stress_tensor_nilpotency(2, [T])
    payload = result.legacy_payload(basis)

    assert payload["status"] == "solved"
    assert payload["n"] == 2
    assert payload["target"] == NO(T, T)
    assert payload["null_operator"] == NO(T, T)
    assert payload["quotient_remainder"] == NO(T, T)
    assert payload["obstruction"] == NO(T, T)
