import shutil

import pytest
import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    MakeOPE,
    NO,
    OPE,
    One,
    dn,
    expand_with_wolfram,
    simplify,
    wolfram_canonicalize_exprs,
    set_compute_backend,
    get_compute_backend,
    compute_backend,
)
from pyope.backend import SUPPORTED_BACKENDS
from pyope.wolfram_backend import (
    WolframBackendError,
    _decode_expr,
    chunk_exprs_for_wolfram,
)


def test_default_backend_is_sympy():
    assert get_compute_backend() == "sympy"


def test_set_compute_backend_rejects_unknown_backend():
    with pytest.raises(ValueError):
        set_compute_backend("unknown")


def test_supported_backends_include_wolfram():
    assert SUPPORTED_BACKENDS == {"sympy", "wolfram"}


def test_compute_backend_context_manager_restores_previous_backend():
    set_compute_backend("sympy")
    with compute_backend("wolfram"):
        assert get_compute_backend() == "wolfram"
    assert get_compute_backend() == "sympy"


def test_sympy_backend_still_computes_normal_ordered_product():
    T = BasisOperator("T")
    J = BasisOperator("J")
    set_compute_backend("sympy")

    result = NO(T, J)

    assert repr(result) == "NO(T, J)"


def test_sympy_backend_still_computes_zero_ope():
    T = BasisOperator("T")
    J = BasisOperator("J")
    set_compute_backend("sympy")

    result = OPE(T, J)

    assert result.max_pole == 0


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_backend_computes_registered_binary_ope():
    T = BasisOperator("T")
    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c / 2 * One, 0, 2 * T, 0])

    with compute_backend("wolfram"):
        result = OPE(T, T)

    assert result.pole(4) == c * One / 2
    assert result.pole(2) == 2 * T


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_backend_computes_binary_no_node():
    T = BasisOperator("T")
    J = BasisOperator("J")

    with compute_backend("wolfram"):
        result = NO(T, J)

    assert repr(result) == "NO(T, J)"


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_backend_preserves_linear_no_structure():
    T = BasisOperator("T")
    J = BasisOperator("J")

    with compute_backend("wolfram"):
        result = NO(T + J, J)

    expected = NO(T, J) + NO(J, J)
    assert sorted(map(repr, sp.Add.make_args(result))) == sorted(
        map(repr, sp.Add.make_args(expected))
    )


def test_wolfram_decoder_rejects_illegal_operator_multiplication():
    T = BasisOperator("T")
    J = BasisOperator("J")

    with pytest.raises(WolframBackendError, match="Illegal operator multiplication"):
        _decode_expr("(T * J)", {"T", "J"})


def test_wolfram_decoder_preserves_derivative_protocol_shape():
    T = BasisOperator("T")

    result = _decode_expr("dn(2, T)", {"T"})

    assert repr(result) == "d^2(T)"


def test_wolfram_decoder_supports_nested_container_payloads():
    T = BasisOperator("T")
    J = BasisOperator("J")

    result = _decode_expr(
        'PyDict([("ops", [NO(T, J), PyTuple([dn(1, T), "tag"])])])',
        {"T", "J"},
    )

    assert result == {"ops": [NO(T, J), (dn(1, T), "tag")]}


def test_chunk_exprs_for_wolfram_respects_item_limit():
    exprs = [sp.Integer(i) for i in range(5)]

    chunks = chunk_exprs_for_wolfram(exprs, max_items=2, max_chars=100)

    assert chunks == [[0, 1], [2, 3], [4]]


def test_chunk_exprs_for_wolfram_raises_for_oversized_expression():
    x = sp.Symbol("x" * 40)

    with pytest.raises(WolframBackendError, match="transport size limit"):
        chunk_exprs_for_wolfram([x], max_items=4, max_chars=20)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_expand_with_wolfram_round_trips_expression_lists():
    T = BasisOperator("Tbatcheval")
    J = BasisOperator("Jbatcheval")

    exprs = [NO(T, J), dn(1, T) + J]

    result = expand_with_wolfram(exprs)

    assert isinstance(result, list)
    assert len(result) == 2
    assert repr(result[0]) == "NO(Tbatcheval, Jbatcheval)"
    assert result[1] == J + dn(1, T)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_expand_with_wolfram_round_trips_nested_containers():
    T = BasisOperator("Tnested")
    J = BasisOperator("Jnested")

    payload = {
        "ops": [NO(T, J), (dn(1, T), J + dn(1, T))],
        "meta": {"label": "trial"},
    }

    result = expand_with_wolfram(payload)

    assert result["ops"][0] == NO(T, J)
    assert result["ops"][1] == (dn(1, T), J + dn(1, T))
    assert result["meta"] == {"label": "trial"}


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_expand_with_wolfram_plus_simplify_round_trips_nested_containers():
    T = BasisOperator("Tsimplify_nested")
    J = BasisOperator("Jsimplify_nested")
    Bosonic(T, J)

    payload = {
        "ops": [NO(T, J) + NO(T, J), (NO(T, J) + NO(T, J), NO(T, J))],
        "meta": {"label": "trial"},
    }

    expanded = expand_with_wolfram(payload)
    result = {
        "ops": [simplify(item) for item in expanded["ops"][:1]]
        + [tuple(simplify(part) for part in expanded["ops"][1])],
        "meta": expanded["meta"],
    }

    assert result == {
        "ops": [2 * NO(T, J), (2 * NO(T, J), NO(T, J))],
        "meta": {"label": "trial"},
    }


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_batch_canonicalize_reorders_normal_products():
    A = BasisOperator("Abatchcanon", conformal_weight=1)
    B = BasisOperator("Bbatchcanon", conformal_weight=1)
    Bosonic(A, B)

    result = wolfram_canonicalize_exprs([NO(B, A), NO(A, B)])

    assert len(result) == 2
    assert result[0] == NO(A, B)
    assert result[1] == NO(A, B)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_backend_round_trips_derivative_payloads():
    T = BasisOperator("T")

    with compute_backend("wolfram"):
        result = NO(sp.Rational(3, 2) * T, dn(2, T))

    assert result == sp.Rational(3, 2) * NO(T, dn(2, T))


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_backend_round_trips_nested_no_protocol_shape():
    T = BasisOperator("T")
    J = BasisOperator("J")

    expr = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(J, J))))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, dn(1, J))))
        + sp.Rational(2, 3) * NO(T, NO(T, NO(dn(1, T), J)))
    )

    with compute_backend("wolfram"):
        result = NO(T, expr)

    expected = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(T, NO(J, J)))))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, NO(T, dn(1, J)))))
        + sp.Rational(2, 3) * NO(T, NO(T, NO(T, NO(dn(1, T), J))))
    )
    assert result == expected


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_whole_expression_eval_preserves_nested_no_shape():
    T = BasisOperator("T")
    J = BasisOperator("J")

    expr = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(J, J))))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, dn(1, J))))
        + sp.Rational(2, 3) * NO(dn(1, T), NO(T, NO(T, J)))
    )

    result = expand_with_wolfram(expr)

    assert result != 0
    assert "NO(" in repr(result)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_expand_with_wolfram_plus_simplify_handles_null_state_style_expr():
    b = BasisOperator("b")
    c = BasisOperator("c")
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")

    from pyope import Bosonic, Fermionic, clear_registry

    clear_registry()
    Fermionic(b, c)
    Bosonic(beta, gamma)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = 3 * NO(beta, gamma) + 2 * NO(b, c)
    G = NO(b, gamma)
    Gbar = 3 * NO(beta, dn(1, c)) + 2 * NO(dn(1, beta), c)
    T = (
        -sp.Rational(3, 2) * NO(beta, dn(1, gamma))
        - sp.Rational(1, 2) * NO(dn(1, beta), gamma)
        - 2 * NO(b, dn(1, c))
        - NO(dn(1, b), c)
    )
    W = beta
    Wbar = (
        sp.Rational(10, 9) * dn(2, gamma)
        - sp.Rational(4, 3) * NO(gamma, b, dn(1, c))
        + sp.Rational(2, 3) * NO(gamma, dn(1, b), c)
        - sp.Rational(8, 3) * NO(dn(1, gamma), b, c)
        - 4 * NO(beta, gamma, dn(1, gamma))
        + sp.Rational(2, 3) * NO(dn(1, beta), gamma, gamma)
        + 2 * NO(beta, gamma, gamma, b, c)
        + NO(beta, beta, gamma, gamma, gamma)
    )

    expr = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(J, J))))
        + sp.Rational(4, 3) * NO(T, NO(T, NO(T, T)))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, dn(1, J))))
        - sp.Rational(8, 9) * NO(T, NO(T, NO(J, NO(G, Gbar))))
        + sp.Rational(2, 3) * NO(T, NO(T, NO(G, dn(1, Gbar))))
        - 3 * NO(T, NO(T, NO(W, dn(1, Wbar))))
        - sp.Rational(1, 9) * NO(T, NO(T, NO(dn(1, J), NO(J, J))))
        - sp.Rational(4, 9) * NO(T, NO(T, NO(dn(1, J), dn(1, J))))
        - sp.Rational(2, 9) * NO(T, NO(T, NO(dn(1, G), Gbar)))
        + 5 * NO(T, NO(T, NO(dn(1, W), Wbar)))
        + NO(T, NO(T, NO(dn(2, J), J)))
    )

    result = simplify(expand_with_wolfram(expr))

    assert result == 0
