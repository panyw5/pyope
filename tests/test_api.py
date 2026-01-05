"""
Unit tests for API functions (OPE, NO, bracket).
"""

import pytest
import sympy as sp
from pyope.api import OPE, NO, bracket
from pyope.operators import BasisOperator, d
from pyope.ope_data import OPEData
from pyope.constants import One
from pyope.registry import ope_registry


@pytest.fixture(autouse=True)
def clear_ope_registry():
    """Clear OPE registry before each test to ensure test isolation."""
    ope_registry.clear()
    yield
    ope_registry.clear()


class TestOPEFunction:
    """Tests for OPE function."""

    def test_ope_linearity_right(self):
        """Test OPE linearity: OPE(A, B+C) = OPE(A,B) + OPE(A,C)."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        K = BasisOperator("K", bosonic=True)

        # Define OPEs
        OPE[T, J] = OPEData({1: J})
        OPE[T, K] = OPEData({1: K})

        # Test linearity
        result = OPE(T, J + K)
        expected = OPE(T, J) + OPE(T, K)

        assert result == expected

    def test_ope_linearity_left(self):
        """Test OPE linearity: OPE(A+B, C) = OPE(A,C) + OPE(B,C)."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        K = BasisOperator("K", bosonic=True)

        # Define OPEs
        OPE[T, K] = OPEData({1: T})
        OPE[J, K] = OPEData({1: J})

        # Test linearity
        result = OPE(T + J, K)
        expected = OPE(T, K) + OPE(J, K)

        assert result == expected

    def test_ope_scalar_multiplication(self):
        """Test OPE with scalar: OPE(c*A, B) = c*OPE(A,B)."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        c = sp.Symbol("c")

        # Define OPE
        OPE[T, J] = OPEData({1: J})

        # Test scalar multiplication
        result = OPE(c * T, J)
        expected = c * OPE(T, J)

        assert result == expected

    def test_ope_with_zero(self):
        """Test OPE with zero."""
        T = BasisOperator("T", bosonic=True)

        result = OPE(T, 0)
        assert result.is_zero()

        result = OPE(0, T)
        assert result.is_zero()

    def test_ope_undefined(self):
        """Test OPE returns zero for undefined pairs."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        # No OPE defined
        result = OPE(T, J)
        assert result.is_zero()

    def test_ope_make(self):
        """Test OPE.make() function."""
        T = BasisOperator("T", bosonic=True)
        dT = d(T)
        c = sp.Symbol("c")

        # Create OPE from list
        ope = OPE.make([c / 2 * One, 0, 2 * T, dT])

        assert ope.pole(4) == c / 2 * One
        assert ope.pole(3) == 0
        assert ope.pole(2) == 2 * T
        assert ope.pole(1) == dT


class TestNOFunction:
    """Tests for NO (normal ordering) function."""

    def test_no_creation(self):
        """Test creating normal ordered products."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        result = NO(T, J)

        assert result.left == T
        assert result.right == J

    def test_no_linearity(self):
        """Test NO linearity: NO(A, B+C) = NO(A,B) + NO(A,C)."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        K = BasisOperator("K", bosonic=True)

        result = NO(T, J + K)
        # Should expand to NO(T,J) + NO(T,K)
        assert isinstance(result, sp.Add)

    def test_no_scalar_multiplication(self):
        """Test NO with scalar: NO(c*A, B) = c*NO(A,B)."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)
        c = sp.Symbol("c")

        result = NO(c * T, J)
        # Should be c * NO(T, J)
        assert isinstance(result, sp.Mul)

    def test_no_with_zero(self):
        """Test NO with zero."""
        T = BasisOperator("T", bosonic=True)

        result = NO(T, 0)
        assert result == 0

        result = NO(0, T)
        assert result == 0


class TestBracket:
    """Tests for bracket (commutator/anticommutator) function."""

    def test_commutator(self):
        """Test commutator [A,B] = AB - BA."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        result = bracket(T, J, anticommutator=False)
        # Should be NO(T,J) - NO(J,T)
        assert isinstance(result, sp.Add)

    def test_anticommutator(self):
        """Test anticommutator {A,B} = AB + BA."""
        psi = BasisOperator("ψ", bosonic=False)
        chi = BasisOperator("χ", bosonic=False)

        result = bracket(psi, chi, anticommutator=True)
        # Should be NO(psi,chi) + NO(chi,psi)
        assert isinstance(result, sp.Add)
