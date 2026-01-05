"""
Unit tests for OPEData class.
"""

import pytest
import sympy as sp
from pyope.ope_data import OPEData
from pyope.operators import BasisOperator


class TestOPEData:
    """Tests for OPEData class."""

    def test_create_empty_ope(self):
        """Test creating an empty OPE."""
        ope = OPEData()
        assert ope.is_zero()
        assert ope.max_pole == 0
        assert len(ope.poles) == 0

    def test_create_from_dict(self):
        """Test creating OPE from poles dictionary."""
        T = BasisOperator("T", bosonic=True)
        poles = {
            4: sp.Rational(1, 2),
            2: 2 * T,
            1: BasisOperator("T'", bosonic=True),
        }
        ope = OPEData(poles)

        assert ope.max_pole == 4
        assert ope.pole(4) == sp.Rational(1, 2)
        assert ope.pole(2) == 2 * T
        assert ope.pole(3) == 0  # Missing pole

    def test_create_from_list(self):
        """Test creating OPE from list (Mathematica style)."""
        T = BasisOperator("T", bosonic=True)
        dT = BasisOperator("T'", bosonic=True)
        c = sp.Symbol("c")

        # List: [pole_4, pole_3, pole_2, pole_1]
        pole_list = [c / 2, 0, 2 * T, dT]
        ope = OPEData.from_list(pole_list)

        assert ope.max_pole == 4
        assert ope.pole(4) == c / 2
        assert ope.pole(3) == 0
        assert ope.pole(2) == 2 * T
        assert ope.pole(1) == dT

    def test_pole_access(self):
        """Test accessing poles."""
        T = BasisOperator("T", bosonic=True)
        ope = OPEData({2: T, 1: 2 * T})

        assert ope.pole(2) == T
        assert ope.pole(1) == 2 * T
        assert ope.pole(3) == 0  # Non-existent pole
        assert ope.pole(0) == 0  # Regular term (not implemented)

    def test_set_pole(self):
        """Test setting poles."""
        ope = OPEData()
        T = BasisOperator("T", bosonic=True)

        ope.set_pole(2, T)
        assert ope.pole(2) == T
        assert ope.max_pole == 2

        # Setting zero removes the pole
        ope.set_pole(2, 0)
        assert ope.pole(2) == 0
        assert ope.is_zero()

    def test_addition(self):
        """Test adding OPEData objects."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        ope1 = OPEData({2: T, 1: J})
        ope2 = OPEData({2: T, 1: -J})

        result = ope1 + ope2
        assert result.pole(2) == 2 * T
        assert result.pole(1) == 0  # J - J = 0

    def test_scalar_multiplication(self):
        """Test multiplying OPE by scalar."""
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        ope = OPEData({2: T, 1: 2 * T})
        result = c * ope

        assert result.pole(2) == c * T
        assert result.pole(1) == 2 * c * T

    def test_negation(self):
        """Test negating OPE."""
        T = BasisOperator("T", bosonic=True)
        ope = OPEData({2: T})

        result = -ope
        assert result.pole(2) == -T

    def test_subtraction(self):
        """Test subtracting OPEData objects."""
        T = BasisOperator("T", bosonic=True)

        ope1 = OPEData({2: 3 * T})
        ope2 = OPEData({2: T})

        result = ope1 - ope2
        assert result.pole(2) == 2 * T

    def test_equality(self):
        """Test OPE equality."""
        T = BasisOperator("T", bosonic=True)

        ope1 = OPEData({2: T, 1: 2 * T})
        ope2 = OPEData({2: T, 1: 2 * T})
        ope3 = OPEData({2: T})

        assert ope1 == ope2
        assert ope1 != ope3

    def test_simplify(self):
        """Test simplifying OPE."""
        T = BasisOperator("T", bosonic=True)
        x = sp.Symbol("x")

        ope = OPEData({2: x**2 + 2 * x + 1})
        simplified = ope.simplify(sp.factor)

        # (x+1)^2 after factoring
        assert simplified.pole(2) == (x + 1) ** 2

    def test_sum_support(self):
        """Test that OPEData works with sum()."""
        T = BasisOperator("T", bosonic=True)

        opes = [
            OPEData({2: T}),
            OPEData({2: T}),
            OPEData({2: T}),
        ]

        result = sum(opes)
        assert result.pole(2) == 3 * T
