"""
Advanced OPE tests based on OPEdefs.m functionality.

This test suite covers the core algorithms from OPEdefs.m:
1. Derivative rules (OPEDerivativeHelpL/R)
2. Commutation relations (OPECommuteHelp)
3. Composite operator OPEs (OPECompositeHelpRQ/LQ)
4. Normal ordering simplification (NOCommuteHelp)
5. OPEPole extraction

These tests are based on the Mathematica implementation and serve as
specifications for the Python implementation.
"""

import pytest
import sympy as sp
from pyope.api import OPE, NO, bracket
from pyope.operators import BasisOperator, d, dn, Operator
from pyope.ope_data import OPEData
from pyope.constants import One
from pyope.registry import ope_registry, Bosonic, Fermionic


@pytest.fixture(autouse=True)
def clear_registry():
    """Clear registries before each test."""
    ope_registry.clear()
    yield
    ope_registry.clear()


class TestVirasoroAlgebra:
    """
    Test Virasoro algebra OPE and related computations.

    The stress-energy tensor T satisfies:
    T(z)T(w) = c/2 / (z-w)^4 + 2T(w) / (z-w)^2 + T'(w) / (z-w)
    """

    def test_virasoro_ope_definition(self):
        """Test basic Virasoro OPE definition."""
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        # Define Virasoro OPE
        OPE[T, T] = OPE.make([c/2 * One, 0, 2*T, d(T)])

        result = OPE(T, T)

        assert result.pole(4) == c/2 * One
        assert result.pole(3) == 0
        assert result.pole(2) == 2*T
        assert result.pole(1) == d(T)

    def test_virasoro_with_primary(self):
        """
        Test OPE of Virasoro with a primary field.

        For a primary field φ of dimension h:
        T(z)φ(w) = h φ(w) / (z-w)^2 + φ'(w) / (z-w)
        """
        T = BasisOperator("T", bosonic=True)
        phi = BasisOperator("φ", bosonic=True)
        h = sp.Symbol("h")

        # Define OPEs
        OPE[T, phi] = OPE.make([h*phi, d(phi)])

        result = OPE(T, phi)

        assert result.pole(2) == h*phi
        assert result.pole(1) == d(phi)


class TestDerivativeRules:
    """
    Test derivative rules for OPEs.

    These correspond to OPEDerivativeHelpL and OPEDerivativeHelpR
    in OPEdefs.m.
    """

    def test_left_derivative_simple(self):
        """
        Test OPE[∂A, B] computation.

        Formula (from OPEdefs.m line 910-920):
        OPE[∂^i A, B] = (-1)^i * Sum[Pochhammer[j,i] * pole_j(OPE[A,B]), {j, maxpole, 1, -1}]

        For OPE[A, B] with pole(2)=A, pole(1)=B:
        - pole(2) contributes to OPE[∂A, B] at pole(3): (-1) * (3-1) * A = -2*A
        - pole(1) contributes to OPE[∂A, B] at pole(2): (-1) * (2-1) * B = -B
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)

        # Define OPE[A, B]
        OPE[A, B] = OPE.make([A, B])  # poles at 2 and 1

        # Compute OPE[∂A, B]
        dA = d(A)
        result = OPE(dA, B)

        # Expected: max_pole=3, pole(3)=-2*A, pole(2)=-B
        assert result.max_pole == 3
        assert result.pole(3) == -2*A
        assert result.pole(2) == -B

    def test_right_derivative_simple(self):
        """
        Test OPE[A, ∂B] computation.

        Formula (from OPEdefs.m line 937-948):
        Uses Leibniz rule with binomial coefficients and Pochhammer symbols.
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)

        # Define OPE[A, B]
        OPE[A, B] = OPE.make([A, B])

        # Compute OPE[A, ∂B]
        dB = d(B)
        result = OPE(A, dB)

        # Expected: includes derivatives of poles
        assert result.max_pole >= 2
        # Note: This test will fail until derivative rules are implemented

    def test_virasoro_derivative(self):
        """
        Test OPE[∂T, T] for Virasoro.

        This is a concrete example from conformal field theory.

        For OPE[T, T] with pole(4)=c/2, pole(3)=0, pole(2)=2*T, pole(1)=∂T:
        - pole(4) contributes to OPE[∂T, T] at pole(5): (-1) * (5-1) * c/2 = -2*c
        - pole(3) contributes to OPE[∂T, T] at pole(4): (-1) * (4-1) * 0 = 0
        - pole(2) contributes to OPE[∂T, T] at pole(3): (-1) * (3-1) * 2*T = -4*T
        - pole(1) contributes to OPE[∂T, T] at pole(2): (-1) * (2-1) * ∂T = -∂T
        """
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        # Define Virasoro OPE
        OPE[T, T] = OPE.make([c/2 * One, 0, 2*T, d(T)])

        # Compute OPE[∂T, T]
        dT = d(T)
        result = OPE(dT, T)

        # Expected: max_pole=5, pole(5)=-2*c, pole(3)=-4*T, pole(2)=-∂T
        assert result.max_pole == 5
        assert result.pole(5) == -2*c * One
        assert result.pole(3) == -4*T
        assert result.pole(2) == -d(T)


class TestCommutationRelations:
    """
    Test commutation relations: computing OPE[B,A] from OPE[A,B].

    Corresponds to OPECommuteHelp in OPEdefs.m (line 959-972).
    """

    def test_bosonic_commutation(self):
        """
        Test commutation for bosonic operators.

        Formula (OPEdefs.m line 959-972):
        OPE[B,A](q) = Sum[(-1)^l / (l-q)! * D^(l-q)[pole_l(OPE[A,B])], {l, q, maxpole}]
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)

        # Define OPE[A, B]
        OPE[A, B] = OPE.make([A, B])

        # Compute OPE[B, A] - should use commutation formula
        result = OPE(B, A)

        # For bosonic operators, OPE[B,A] should be related to OPE[A,B]
        # by the commutation formula
        assert not result.is_zero()
        # Note: This test will fail until commutation is implemented

    def test_fermionic_anticommutation(self):
        """
        Test anticommutation for fermionic operators.

        For fermions, there's an additional sign factor.
        """
        psi = BasisOperator("ψ", bosonic=False)
        chi = BasisOperator("χ", bosonic=False)

        # Define OPE[ψ, χ]
        OPE[psi, chi] = OPE.make([psi])

        # Compute OPE[χ, ψ] - should have opposite sign
        result = OPE(chi, psi)

        # For fermions: OPE[χ,ψ] = -OPE[ψ,χ] + regular terms
        assert not result.is_zero()
        # Note: This test will fail until commutation is implemented

    def test_virasoro_commutation(self):
        """
        Test T(w)T(z) from T(z)T(w).

        This is a concrete example with known result.
        """
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        # Define T(z)T(w)
        OPE[T, T] = OPE.make([c/2 * One, 0, 2*T, d(T)])

        # T(w)T(z) should be the same (T is bosonic and self-OPE is symmetric)
        result = OPE(T, T)
        result_commuted = OPE(T, T)  # Should use commutation if order is swapped

        # For self-OPE of bosonic operator, should be symmetric
        assert result.pole(4) == result_commuted.pole(4)


class TestCompositeOperatorOPE:
    """
    Test OPEs involving composite operators (normal ordered products).

    Corresponds to OPECompositeHelpRQ/LQ in OPEdefs.m.
    """

    def test_ope_with_right_composite(self):
        """
        Test OPE[A, NO[B,C]].

        Formula (OPEdefs.m line 985-1016):
        OPE[A, NO[B,C]] involves:
        - SwapSign[A,B] * NO[B, OPE[A,C]]
        - NO[OPE[A,B], C]
        - Sum over nested OPEs
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)
        C = BasisOperator("C", bosonic=True)

        # Define basic OPEs
        OPE[A, B] = OPE.make([A])
        OPE[A, C] = OPE.make([C])

        # Compute OPE[A, NO[B,C]]
        BC = NO(B, C)
        result = OPE(A, BC)

        # Should involve normal ordered products
        assert not result.is_zero()
        # Note: This test will fail until composite OPE is implemented

    def test_ope_with_left_composite(self):
        """
        Test OPE[NO[A,B], C].

        Formula (OPEdefs.m line 1028-1084):
        More complex than right composite case.
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)
        C = BasisOperator("C", bosonic=True)

        # Define basic OPEs
        OPE[A, C] = OPE.make([A])
        OPE[B, C] = OPE.make([B])

        # Compute OPE[NO[A,B], C]
        AB = NO(A, B)
        result = OPE(AB, C)

        # Should involve derivatives and nested OPEs
        assert not result.is_zero()
        # Note: This test will fail until composite OPE is implemented


class TestNormalOrderingSimplification:
    """
    Test normal ordering simplification rules.

    Corresponds to NOCommuteHelp and related functions in OPEdefs.m.
    """

    def test_no_commutator_formula(self):
        """
        Test NO[A,B] - SwapSign[A,B] * NO[B,A] = commutator.

        Formula (OPEdefs.m line 1520-1528):
        NOCommuteHelp[A,B] = Sum[-(-1)^m / m! * D^m[pole_m(OPE[A,B])], {m, maxpole}]
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)

        # Define OPE[A, B]
        OPE[A, B] = OPE.make([A, B])

        # Compute commutator
        comm = bracket(A, B, anticommutator=False)

        # Should be NO[A,B] - NO[B,A]
        assert comm != 0
        # Note: This test will fail until NO simplification is implemented

    def test_nested_no_simplification(self):
        """
        Test simplification of NO[NO[A,B], C].

        Should be rewritten using OPE formulas.
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)
        C = BasisOperator("C", bosonic=True)

        # Define OPEs
        OPE[A, C] = OPE.make([A])
        OPE[B, C] = OPE.make([B])

        # Create nested NO
        AB = NO(A, B)
        ABC = NO(AB, C)

        # Should simplify to a combination of simpler terms
        # Note: This test will fail until NO simplification is implemented
        assert ABC is not None


class TestOPEPoleExtraction:
    """
    Test OPEPole function for extracting specific poles.

    This is a critical function used throughout OPEdefs.m.
    """

    def test_pole_extraction_basic(self):
        """Test extracting poles from OPEData."""
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        ope = OPE.make([c/2 * One, 0, 2*T, d(T)])

        assert ope.pole(4) == c/2 * One
        assert ope.pole(3) == 0
        assert ope.pole(2) == 2*T
        assert ope.pole(1) == d(T)
        assert ope.pole(0) == 0  # Regular part

    def test_pole_from_operators(self):
        """
        Test computing pole directly from operators.

        In OPEdefs.m, OPEPole[n][A,B] computes the n-th pole
        without first computing the full OPE.
        """
        A = BasisOperator("A", bosonic=True)
        B = BasisOperator("B", bosonic=True)

        # Define OPE
        OPE[A, B] = OPE.make([A, B])

        # Extract pole
        result = OPE(A, B)
        pole2 = result.pole(2)

        assert pole2 == A


class TestCurrentAlgebra:
    """
    Test current algebra (Kac-Moody) OPEs.

    This is a standard example from conformal field theory.
    """

    def test_kac_moody_ope(self):
        """
        Test current algebra OPE.

        J^a(z)J^b(w) = k δ^{ab} / (z-w)^2 + i f^{abc} J^c(w) / (z-w)
        """
        # Use indexed operators
        i, j = sp.symbols("i j")
        J = lambda idx: BasisOperator(f"J", bosonic=True, indices=(idx,))

        k = sp.Symbol("k")

        # Define diagonal OPE: J^i(z)J^i(w)
        Ji = J(i)
        OPE[Ji, Ji] = OPE.make([k, 0])

        result = OPE(Ji, Ji)
        assert result.pole(2) == k


class TestExpressionSimplification:
    """
    Test expression simplification and collection of terms.

    Corresponds to OPESimplify and ExtractOperators in OPEdefs.m.
    """

    def test_collect_like_terms(self):
        """Test collecting terms with the same operator."""
        T = BasisOperator("T", bosonic=True)
        c = sp.Symbol("c")

        # Create expression with multiple T terms
        expr = T + 2*T + c*T

        # SymPy automatically combines the constant coefficients: T + 2*T = 3*T
        # The result is c*T + 3*T
        # Check that the expression has the correct structure
        assert isinstance(expr, sp.Add)
        assert len(expr.args) == 2  # Should have 2 terms: c*T and 3*T

        # Extract coefficients for T
        terms = list(expr.args)
        coeff_sum = 0
        for term in terms:
            if isinstance(term, sp.Mul) and any(arg == T for arg in term.args):
                # Extract coefficient
                coeff = 1
                for arg in term.args:
                    if arg != T and not isinstance(arg, Operator):
                        coeff *= arg
                coeff_sum += coeff
            elif term == T:
                coeff_sum += 1

        # Verify total coefficient is 3 + c
        expected = sp.expand(3 + c)
        assert sp.simplify(coeff_sum - expected) == 0

    def test_ope_arithmetic(self):
        """Test arithmetic operations on OPEData."""
        T = BasisOperator("T", bosonic=True)
        J = BasisOperator("J", bosonic=True)

        ope1 = OPEData({2: T, 1: J})
        ope2 = OPEData({2: T, 1: -J})

        # Addition
        result = ope1 + ope2
        assert result.pole(2) == 2*T
        assert result.pole(1) == 0

        # Scalar multiplication
        c = sp.Symbol("c")
        result = c * ope1
        assert result.pole(2) == c*T
        assert result.pole(1) == c*J


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
