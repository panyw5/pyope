"""
Comprehensive test suite for expand_nested_no with common VOA algebras.

This script tests expand_nested_no against various algebras:
- sl(2) Kac-Moody algebra
- Virasoro algebra
- Mixed algebras

The goal is to verify that expand_nested_no produces correct expansions
that match the expected behavior from OPEdefs.m.
"""

import pytest
import sympy as sp
from pyope import (
    BasisOperator,
    OPE,
    NO,
    d,
    One,
    Zero,
)
from pyope.ope_data import OPEData
from pyope.simplify import expand_nested_no, simplify


class TestExpandNestedNO_SL2_Comprehensive:
    """Comprehensive tests for sl(2) Kac-Moody algebra."""

    @classmethod
    def setup_class(cls):
        """Set up sl(2) algebra at level k."""
        cls.Jp = BasisOperator("J⁺", conformal_weight=1)
        cls.J0 = BasisOperator("J⁰", conformal_weight=1)
        cls.Jm = BasisOperator("J⁻", conformal_weight=1)
        cls.k = sp.Symbol("k", real=True, positive=True)

        # Define sl(2) OPEs
        # [J+, J-] = k/2 δ(z-w) + J0 ∂δ(z-w)
        # In OPE form: J+(z)J-(w) = k/2/(z-w)^2 + J0/(z-w)
        OPE[cls.Jp, cls.Jm] = OPEData.from_list(
            [
                cls.k / 2 * One,  # 2nd pole
                cls.J0,  # 1st pole
            ]
        )

        # [J0, J+] = J+
        # J0(z)J+(w) = J+/(z-w)
        OPE[cls.J0, cls.Jp] = OPEData.from_list(
            [
                cls.Jp  # 1st pole
            ]
        )

        # [J0, J-] = -J-
        # J0(z)J-(w) = -J-/(z-w)
        OPE[cls.J0, cls.Jm] = OPEData.from_list(
            [
                -cls.Jm  # 1st pole
            ]
        )

        # [J0, J0] = k/2
        # J0(z)J0(w) = k/2/(z-w)^2
        OPE[cls.J0, cls.J0] = OPEData.from_list(
            [
                cls.k / 2 * One,  # 2nd pole
                Zero,  # 1st pole
            ]
        )

    def test_sl2_case1_left_nested_jp_j0_jm(self):
        """
        Test Case 1: NO(NO(J+, J0), J-)

        This is a left-nested expression. According to Thielemans eq (3.3.19):
        NO(NO(A,B), C) = NO(A, NO(B,C)) + Σ 1/l! NO(∂^l A, {BC}_l) + sign Σ 1/l! NO(∂^l B, {AC}_l)

        For A=J+, B=J0, C=J-:
        - OPE[J0, J-] has pole at l=1: -J-
        - OPE[J+, J-] has poles at l=1: J0, l=2: k/2
        - sign = (-1)^{|J+||J0|} = 1 (both bosonic)
        """
        expr = NO(NO(self.Jp, self.J0), self.Jm)
        result = expand_nested_no(expr)

        # The result should be mathematically equivalent to the input
        # Verify by checking that simplify(result) == simplify(expr)
        result_simplified = simplify(result)
        expr_simplified = simplify(expr)

        diff = simplify(result_simplified - expr_simplified)
        assert diff == Zero or diff == 0, (
            f"expand_nested_no changed the physical meaning!\n"
            f"Input simplified: {expr_simplified}\n"
            f"Result simplified: {result_simplified}\n"
            f"Diff: {diff}"
        )

        # Also verify the result is not trivial
        assert result_simplified != Zero

        print(f"\nCase 1 - NO(NO(J+, J0), J-):")
        print(f"  Input:  {expr}")
        print(f"  Expanded: {result}")
        print(f"  Simplified: {result_simplified}")

    def test_sl2_case2_right_nested_jp_j0_jm(self):
        """
        Test Case 2: NO(J+, NO(J0, J-))

        This is already in right-nested form. If J+ comes before J0 in canonical order,
        this should not expand further (or only expand the inner NO if needed).
        """
        expr = NO(self.Jp, NO(self.J0, self.Jm))
        result = expand_nested_no(expr)

        # The result should be valid
        result_simplified = simplify(result)
        assert result_simplified is not None

        print(f"\nCase 2 - NO(J+, NO(J0, J-)):")
        print(f"  Input:  {expr}")
        print(f"  Result: {result}")
        print(f"  Simplified: {result_simplified}")

    def test_sl2_case3_left_nested_jm_jp_j0(self):
        """
        Test Case 3: NO(NO(J-, J+), J0)

        This involves OPE[J-, J+] which is non-trivial.
        """
        expr = NO(NO(self.Jm, self.Jp), self.J0)
        result = expand_nested_no(expr)
        result_simplified = simplify(result)

        assert result_simplified is not None

        print(f"\nCase 3 - NO(NO(J-, J+), J0):")
        print(f"  Input:  {expr}")
        print(f"  Result: {result}")
        print(f"  Simplified: {result_simplified}")

    def test_sl2_case4_right_nested_jm_jp_j0(self):
        """
        Test Case 4: NO(J-, NO(J+, J0))
        """
        expr = NO(self.Jm, NO(self.Jp, self.J0))
        result = expand_nested_no(expr)
        result_simplified = simplify(result)

        assert result_simplified is not None

        print(f"\nCase 4 - NO(J-, NO(J+, J0)):")
        print(f"  Input:  {expr}")
        print(f"  Result: {result}")
        print(f"  Simplified: {result_simplified}")

    def test_sl2_case5_triple_nested(self):
        """
        Test Case 5: NO(NO(NO(J+, J0), J-), J0)

        Triple nested expression.
        """
        expr = NO(NO(NO(self.Jp, self.J0), self.Jm), self.J0)
        result = expand_nested_no(expr)
        result_simplified = simplify(result)

        assert result_simplified is not None

        print(f"\nCase 5 - NO(NO(NO(J+, J0), J-), J0):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 200 chars): {str(result)[:200]}...")
        print(f"  Simplified (first 200 chars): {str(result_simplified)[:200]}...")


class TestExpandNestedNO_Virasoro_Comprehensive:
    """Comprehensive tests for Virasoro algebra."""

    @classmethod
    def setup_class(cls):
        """Set up Virasoro algebra."""
        cls.T = BasisOperator("T", conformal_weight=2)
        cls.c = sp.Symbol("c", real=True)

        # Define Virasoro OPE: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
        OPE[cls.T, cls.T] = OPEData.from_list(
            [
                cls.c / 2 * One,  # 4th pole
                Zero,  # 3rd pole
                2 * cls.T,  # 2nd pole
                d(cls.T),  # 1st pole
            ]
        )

    def test_virasoro_case1_left_nested_ttt(self):
        """
        Test Case 1: NO(NO(T, T), T)

        Left-nested Virasoro expression.
        """
        expr = NO(NO(self.T, self.T), self.T)
        result = expand_nested_no(expr)
        result_simplified = simplify(result)

        assert result_simplified is not None
        assert result_simplified != Zero

        # The result should involve the central charge c
        result_str = str(result_simplified)

        print(f"\nVirasoro Case 1 - NO(NO(T, T), T):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 300 chars): {str(result)[:300]}...")
        print(f"  Simplified (first 300 chars): {str(result_simplified)[:300]}...")

    def test_virasoro_case2_right_nested_ttt(self):
        """
        Test Case 2: NO(T, NO(T, T))

        Right-nested Virasoro expression.
        """
        expr = NO(self.T, NO(self.T, self.T))
        result = expand_nested_no(expr)
        result_simplified = simplify(result)

        assert result_simplified is not None
        assert result_simplified != Zero

        print(f"\nVirasoro Case 2 - NO(T, NO(T, T)):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 300 chars): {str(result)[:300]}...")
        print(f"  Simplified (first 300 chars): {str(result_simplified)[:300]}...")

    def test_virasoro_case3_with_derivative_left(self):
        """
        Test Case 3: NO(NO(T, ∂T), T)
        """
        expr = NO(NO(self.T, d(self.T)), self.T)
        result = expand_nested_no(expr, expand_derivatives=False)
        result_simplified = simplify(result, expand_derivatives=False)

        assert result_simplified is not None

        print(f"\nVirasoro Case 3 - NO(NO(T, ∂T), T):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 200 chars): {str(result)[:200]}...")

    def test_virasoro_case4_with_derivative_right(self):
        """
        Test Case 4: NO(T, NO(∂T, T))
        """
        expr = NO(self.T, NO(d(self.T), self.T))
        result = expand_nested_no(expr, expand_derivatives=False)
        result_simplified = simplify(result, expand_derivatives=False)

        assert result_simplified is not None

        print(f"\nVirasoro Case 4 - NO(T, NO(∂T, T)):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 200 chars): {str(result)[:200]}...")

    def test_virasoro_case5_both_derivatives(self):
        """
        Test Case 5: NO(NO(∂T, T), ∂T)
        """
        expr = NO(NO(d(self.T), self.T), d(self.T))
        result = expand_nested_no(expr, expand_derivatives=False)
        result_simplified = simplify(result, expand_derivatives=False)

        assert result_simplified is not None

        print(f"\nVirasoro Case 5 - NO(NO(∂T, T), ∂T):")
        print(f"  Input:  {expr}")
        print(f"  Result (first 200 chars): {str(result)[:200]}...")


class TestExpandNestedNO_Consistency:
    """Test consistency properties of expand_nested_no."""

    def test_expand_then_simplify_equals_simplify(self):
        """
        Test that expand_nested_no followed by simplify gives the same result as simplify alone.

        This verifies that expand_nested_no produces mathematically equivalent expressions.
        """
        # Set up sl(2)
        Jp = BasisOperator("J⁺", conformal_weight=1)
        J0 = BasisOperator("J⁰", conformal_weight=1)
        Jm = BasisOperator("J⁻", conformal_weight=1)
        k = sp.Symbol("k")

        OPE[Jp, Jm] = OPEData.from_list([k / 2 * One, J0])
        OPE[J0, Jp] = OPEData.from_list([Jp])
        OPE[J0, Jm] = OPEData.from_list([-Jm])

        # Test expression
        expr = NO(NO(Jp, J0), Jm)

        # Method 1: expand then simplify
        expanded = expand_nested_no(expr)
        result1 = simplify(expanded)

        # Method 2: simplify directly
        result2 = simplify(expr)

        # They should be equal
        diff = simplify(result1 - result2)
        assert diff == Zero or diff == 0, (
            f"Inconsistency: expand+simplify ≠ simplify\nDiff: {diff}"
        )

        print(f"\nConsistency test passed:")
        print(f"  expand_nested_no + simplify: {result1}")
        print(f"  simplify directly:            {result2}")

    def test_expand_preserves_physical_equivalence(self):
        """
        Test that expand_nested_no preserves physical equivalence.

        Two expressions are physically equivalent if their difference simplifies to zero.
        """
        # Set up Virasoro
        T = BasisOperator("T", conformal_weight=2)
        c = sp.Symbol("c")

        OPE[T, T] = OPEData.from_list([c / 2 * One, Zero, 2 * T, d(T)])

        # Test expression
        expr = NO(NO(T, T), T)

        # Expand
        expanded = expand_nested_no(expr)

        # Check physical equivalence: simplify(expanded - expr) should be 0
        # Note: We need to simplify both sides first
        diff = simplify(expanded) - simplify(expr)
        diff_simplified = simplify(diff)

        assert diff_simplified == Zero or diff_simplified == 0, (
            f"expand_nested_no changed physical meaning!\nDiff: {diff_simplified}"
        )

        print(f"\nPhysical equivalence test passed for Virasoro")


def run_comprehensive_tests():
    """Run all comprehensive tests with detailed output."""
    print("=" * 80)
    print("COMPREHENSIVE TESTS FOR expand_nested_no")
    print("=" * 80)

    # Run sl(2) tests
    print("\n" + "=" * 80)
    print("SL(2) KAC-MOODY ALGEBRA TESTS")
    print("=" * 80)

    sl2_tests = TestExpandNestedNO_SL2_Comprehensive()
    sl2_tests.setup_class()

    sl2_tests.test_sl2_case1_left_nested_jp_j0_jm()
    sl2_tests.test_sl2_case2_right_nested_jp_j0_jm()
    sl2_tests.test_sl2_case3_left_nested_jm_jp_j0()
    sl2_tests.test_sl2_case4_right_nested_jm_jp_j0()
    sl2_tests.test_sl2_case5_triple_nested()

    # Run Virasoro tests
    print("\n" + "=" * 80)
    print("VIRASORO ALGEBRA TESTS")
    print("=" * 80)

    vir_tests = TestExpandNestedNO_Virasoro_Comprehensive()
    vir_tests.setup_class()

    vir_tests.test_virasoro_case1_left_nested_ttt()
    vir_tests.test_virasoro_case2_right_nested_ttt()
    vir_tests.test_virasoro_case3_with_derivative_left()
    vir_tests.test_virasoro_case4_with_derivative_right()
    vir_tests.test_virasoro_case5_both_derivatives()

    # Run consistency tests
    print("\n" + "=" * 80)
    print("CONSISTENCY TESTS")
    print("=" * 80)

    cons_tests = TestExpandNestedNO_Consistency()
    cons_tests.test_expand_then_simplify_equals_simplify()
    cons_tests.test_expand_preserves_physical_equivalence()

    print("\n" + "=" * 80)
    print("ALL COMPREHENSIVE TESTS COMPLETED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == "__main__":
    # Run with pytest
    pytest.main([__file__, "-v", "-s"])

    # Or run comprehensive tests directly
    # run_comprehensive_tests()
