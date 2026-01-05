"""
Unit tests for Operator classes.
"""

import pytest
import sympy as sp
from pyope.operators import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    d,
    dn,
)


class TestBasisOperator:
    """Tests for BasisOperator class."""

    def test_create_bosonic_operator(self):
        """Test creating a bosonic operator."""
        T = BasisOperator("T", bosonic=True)
        assert T.name == "T"
        assert T.is_bosonic
        assert not T.is_fermionic
        assert T.parity == 0

    def test_create_fermionic_operator(self):
        """Test creating a fermionic operator."""
        psi = BasisOperator("ψ", bosonic=False)
        assert psi.name == "ψ"
        assert psi.is_fermionic
        assert not psi.is_bosonic
        assert psi.parity == 1

    def test_operator_equality(self):
        """Test operator equality."""
        T1 = BasisOperator("T", bosonic=True)
        T2 = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        assert T1 == T2
        assert T1 != J

    def test_operator_hash(self):
        """Test operator hashing."""
        T1 = BasisOperator("T", bosonic=True)
        T2 = BasisOperator("T", bosonic=True)

        # Same operators should have same hash
        assert hash(T1) == hash(T2)

        # Can be used in dictionaries
        d = {T1: "value"}
        assert d[T2] == "value"

    def test_indexed_operator(self):
        """Test indexed operators."""
        J = BasisOperator("J", bosonic=True, indexed=True)
        i = sp.Symbol("i")
        j = sp.Symbol("j")

        J_i = J[i]
        J_j = J[j]

        assert J_i.base_name == "J"  # Use base_name instead of name
        assert J_i.indices == (i,)
        assert J_j.indices == (j,)
        assert J_i != J_j

    def test_arithmetic_operations(self):
        """Test arithmetic operations on operators."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        # Addition
        expr = T + J
        assert isinstance(expr, sp.Add)

        # Multiplication
        expr = 2 * T
        assert isinstance(expr, sp.Mul)

        # Subtraction
        expr = T - J
        assert isinstance(expr, sp.Add)


class TestDerivativeOperator:
    """Tests for DerivativeOperator class."""

    def test_create_derivative(self):
        """Test creating derivative operators."""
        T = BasisOperator("T", bosonic=True)
        dT = d(T)

        assert dT.base == T
        assert dT.order == 1
        assert dT.parity == T.parity

    def test_higher_order_derivative(self):
        """Test higher order derivatives."""
        T = BasisOperator("T", bosonic=True)
        d2T = d(T, 2)
        d3T = dn(3, T)

        assert d2T.order == 2
        assert d3T.order == 3

    def test_derivative_equality(self):
        """Test derivative equality."""
        T = BasisOperator("T", bosonic=True)
        dT1 = d(T)
        dT2 = d(T)

        assert dT1 == dT2

    def test_derivative_hash(self):
        """Test derivative hashing."""
        T = BasisOperator("T", bosonic=True)
        dT1 = d(T)
        dT2 = d(T)

        assert hash(dT1) == hash(dT2)


class TestNormalOrderedOperator:
    """Tests for NormalOrderedOperator class."""

    def test_create_normal_ordered(self):
        """Test creating normal ordered operators."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        NO_TJ = NormalOrderedOperator(T, J)

        assert NO_TJ.left == T
        assert NO_TJ.right == J
        assert NO_TJ.factors == (T, J)

    def test_normal_ordered_parity(self):
        """Test parity of normal ordered operators."""
        T = BasisOperator("T", bosonic=True)
        psi = BasisOperator("ψ", bosonic=False)

        NO_TT = NormalOrderedOperator(T, T)
        NO_Tpsi = NormalOrderedOperator(T, psi)
        NO_psipsi = NormalOrderedOperator(psi, psi)

        assert NO_TT.parity == 0  # bosonic
        assert NO_Tpsi.parity == 1  # fermionic
        assert NO_psipsi.parity == 0  # bosonic (1+1=2, 2%2=0)

    def test_normal_ordered_equality(self):
        """Test normal ordered equality."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        NO1 = NormalOrderedOperator(T, J)
        NO2 = NormalOrderedOperator(T, J)
        NO3 = NormalOrderedOperator(J, T)

        assert NO1 == NO2
        assert NO1 != NO3
