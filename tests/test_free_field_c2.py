from pyope import (
    AbstractC2Reducer,
    BasisOperator,
    Bosonic,
    DerivativeKillingFreeFieldC2Reducer,
    DerivativeKillingRealizationBackend,
    FreeFieldC2Reducer,
    GenericC2Reducer,
    IdentityRealizationBackend,
    LocalOperatorCanonicalizer,
    RealizationBackend,
    NO,
    d,
)


def test_identity_realization_backend_matches_canonicalizer_interface():
    T = BasisOperator("T_identity_backend", conformal_weight=2)
    Bosonic(T)

    canonicalizer = LocalOperatorCanonicalizer([T], stress_tensor=T, max_weight=5)
    backend = IdentityRealizationBackend(canonicalizer)

    assert isinstance(backend, RealizationBackend)
    assert backend.realize(NO(T, T)) == NO(T, T)
    assert backend.sparse_terms(NO(T, T))[NO(T, T)] == 1


def test_derivative_killing_backend_removes_derivative_terms_in_quotient():
    T = BasisOperator("T_derivative_backend", conformal_weight=2)
    Bosonic(T)

    canonicalizer = LocalOperatorCanonicalizer([T], stress_tensor=T, max_weight=5)
    backend = DerivativeKillingRealizationBackend(canonicalizer)

    assert backend.quotient_normal_form(d(T)) == 0
    assert backend.quotient_normal_form(NO(d(T), T)) == 0
    assert backend.quotient_normal_form(NO(T, d(T))) == 0
    assert backend.quotient_normal_form(NO(T, T)) == NO(T, T)


def test_free_field_c2_reducer_falls_back_to_generic_reducer_by_default():
    T = BasisOperator("T_free_field_fallback", conformal_weight=2)
    Bosonic(T)

    canonicalizer = LocalOperatorCanonicalizer([T], stress_tensor=T, max_weight=5)
    generic = GenericC2Reducer(canonicalizer)
    reducer = FreeFieldC2Reducer(canonicalizer)

    assert isinstance(reducer, AbstractC2Reducer)
    assert reducer.quotient_normal_form(NO(d(T), T)) == generic.quotient_normal_form(
        NO(d(T), T)
    )
    assert reducer.solve_c2_witness(NO(T, T)).remainder == NO(T, T)


def test_derivative_killing_free_field_reducer_uses_backend_rule_before_fallback():
    T = BasisOperator("T_derivative_ff", conformal_weight=2)
    Bosonic(T)

    canonicalizer = LocalOperatorCanonicalizer([T], stress_tensor=T, max_weight=5)
    reducer = DerivativeKillingFreeFieldC2Reducer(canonicalizer)

    assert isinstance(reducer, AbstractC2Reducer)
    assert reducer.quotient_normal_form(d(T)) == 0
    witness = reducer.solve_c2_witness(NO(d(T), T))
    assert witness.remainder == 0
    assert witness.c2_part == NO(d(T), T)
