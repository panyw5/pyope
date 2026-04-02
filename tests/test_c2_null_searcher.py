from pyope import (
    BasisOperator,
    Bosonic,
    C2NullSearcher,
    LocalOperatorBasis,
    NO,
    NullSearchResult,
    d,
)


def test_c2_null_searcher_finds_exact_descendant_target():
    T = BasisOperator("T_c2_search_exact", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=4)
    searcher = C2NullSearcher(basis_builder, stress_tensor=T)

    result = searcher.search_from_sources(4, [T], NO(T, T))

    assert result is not None
    assert result["status"] == "solved"
    assert result["null_operator"] == NO(T, T)
    assert result["c2_remainder"] == 0
    assert result["quotient_remainder"] == NO(T, T)
    assert result["obstruction"] == NO(T, T)


def test_c2_null_searcher_search_stress_tensor_nilpotency_uses_no_product():
    T = BasisOperator("T_c2_search_n", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=4)
    searcher = C2NullSearcher(basis_builder, stress_tensor=T)

    result = searcher.search_stress_tensor_nilpotency(2, [T])

    assert result is not None
    assert result["n"] == 2
    assert result["target"] == NO(T, T)
    assert result["null_operator"] == NO(T, T)
    assert result["status"] == "solved"


def test_c2_null_searcher_returns_none_when_unsolved():
    T = BasisOperator("T_c2_search_none", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=4)
    searcher = C2NullSearcher(basis_builder, stress_tensor=T)

    result = searcher.search_from_sources(4, [], NO(T, T))

    assert result is None


def test_c2_null_searcher_exposes_quotient_precheck_bridge():
    T = BasisOperator("T_c2_precheck_bridge", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=5)
    searcher = C2NullSearcher(basis_builder, stress_tensor=T)

    result = searcher.quotient_precheck(NO(d(T), T))

    assert isinstance(result, NullSearchResult)
    assert result.status == "needs_lift"
    assert result.quotient_remainder == 0


def test_c2_null_searcher_search_from_sources_uses_new_lift_core():
    T = BasisOperator("T_c2_lift_bridge", conformal_weight=2)
    Bosonic(T)

    basis_builder = LocalOperatorBasis([T], max_weight=5)
    searcher = C2NullSearcher(basis_builder, stress_tensor=T)

    result = searcher.search_from_sources(4, [T], NO(T, T))

    assert result is not None
    assert result["null_operator"] == NO(T, T)
    assert result["c2_remainder"] == 0
    assert len(result["descendant_basis"]) == len(result["descendant_coeffs"])
    assert NO(T, T) in result["descendant_basis"]
    assert result["descendant_coeffs"][result["descendant_basis"].index(NO(T, T))] == 1
