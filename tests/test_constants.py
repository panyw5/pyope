"""
Unit tests for constants (One, Delta).
"""

import pytest
import sympy as sp
from pyope.constants import One, Delta, ConstantOperator


class TestConstantOperator:
    """Tests for ConstantOperator class."""

    def test_one_operator(self):
        """Test the One constant operator."""
        assert One.name == "One"
        assert One.is_bosonic
        assert One.parity == 0

    def test_constant_equality(self):
        """Test constant operator equality."""
        one1 = ConstantOperator("One")
        one2 = ConstantOperator("One")
        two = ConstantOperator("Two")

        assert one1 == one2
        assert one1 != two

    def test_constant_hash(self):
        """Test constant operator hashing."""
        one1 = ConstantOperator("One")
        one2 = ConstantOperator("One")

        assert hash(one1) == hash(one2)


class TestDelta:
    """Tests for Delta function."""

    def test_delta_equal_indices(self):
        """Test Delta with equal indices."""
        i = sp.Symbol("i")
        result = Delta(i, i)
        assert result == 1

    def test_delta_different_numbers(self):
        """Test Delta with different numeric indices."""
        result = Delta(1, 2)
        assert result == 0

    def test_delta_same_numbers(self):
        """Test Delta with same numeric indices."""
        result = Delta(1, 1)
        assert result == 1

    def test_delta_symbolic(self):
        """Test Delta with symbolic indices."""
        i = sp.Symbol("i")
        j = sp.Symbol("j")

        result = Delta(i, j)
        # Should remain symbolic
        assert isinstance(result, sp.Function)

    def test_delta_simplify(self):
        """Test Delta simplification."""
        i = sp.Symbol("i")

        expr = Delta(i, i)
        assert expr == 1

        expr = Delta(1, 2)
        assert expr == 0
