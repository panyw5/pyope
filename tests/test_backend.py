import shutil

import pytest
import sympy as sp

from pyope import (
    BasisOperator,
    MakeOPE,
    NO,
    OPE,
    One,
    dn,
    simplify_with_wolfram,
    wolfram_evaluate_expr,
    set_compute_backend,
    get_compute_backend,
    compute_backend,
)
from pyope.backend import SUPPORTED_BACKENDS
from pyope.wolfram_backend import WolframBackendError, _decode_expr


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

    result = wolfram_evaluate_expr(expr)

    assert result != 0
    assert "NO(" in repr(result)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_whole_expression_simplify_handles_null_state_style_expr():
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

    result = simplify_with_wolfram(expr)

    assert result == 0
