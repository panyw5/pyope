"""
Simple example demonstrating PyOPE basic functionality.

This example shows how to:
1. Declare operators
2. Define OPEs
3. Compute OPEs with linearity
4. Create normal ordered products
"""

import sympy as sp
from pyope import BasisOperator, OPE, NO, d
from pyope.ope_data import OPEData
from pyope.constants import One


def example_virasoro():
    """
    Example: Virasoro algebra OPE.

    The stress-energy tensor T has OPE:
    T(z)T(w) = c/2 / (z-w)^4 + 2T(w) / (z-w)^2 + T'(w) / (z-w)
    """
    print("=" * 60)
    print("Example: Virasoro Algebra")
    print("=" * 60)

    # Declare the stress-energy tensor T (bosonic)
    T = BasisOperator("T", bosonic=True)
    print(f"\nDeclared operator: {T}")
    print(f"  - Bosonic: {T.is_bosonic}")
    print(f"  - Parity: {T.parity}")

    # Define the central charge
    c = sp.Symbol("c")

    # Define the OPE T(z)T(w)
    # Using the list format: [pole_4, pole_3, pole_2, pole_1]
    OPE[T, T] = OPE.make([
        c / 2 * One,  # 4th order pole
        0,            # 3rd order pole (zero)
        2 * T,        # 2nd order pole
        d(T)          # 1st order pole (derivative)
    ])

    print(f"\nDefined OPE[T, T]:")
    result = OPE(T, T)
    print(f"  {result}")

    # Access individual poles
    print(f"\nIndividual poles:")
    print(f"  4th order: {result.pole(4)}")
    print(f"  3rd order: {result.pole(3)}")
    print(f"  2nd order: {result.pole(2)}")
    print(f"  1st order: {result.pole(1)}")

    # Test linearity: OPE(T, 2T) = 2 * OPE(T, T)
    print(f"\nTesting linearity:")
    result_2T = OPE(T, 2 * T)
    print(f"  OPE(T, 2T) pole(4) = {result_2T.pole(4)}")
    print(f"  Expected: {2 * result.pole(4)}")
    print(f"  Match: {result_2T.pole(4) == 2 * result.pole(4)}")


def example_current_algebra():
    """
    Example: Current algebra (Kac-Moody).

    For currents J^a, the OPE is:
    J^a(z)J^b(w) = k δ^{ab} / (z-w)^2 + i f^{abc} J^c(w) / (z-w)
    """
    print("\n" + "=" * 60)
    print("Example: Current Algebra (simplified)")
    print("=" * 60)

    # Declare current operators
    J1 = BasisOperator("J1", bosonic=True)
    J2 = BasisOperator("J2", bosonic=True)

    print(f"\nDeclared operators: {J1}, {J2}")

    # Define level k
    k = sp.Symbol("k")

    # Define OPE J1(z)J1(w) (diagonal case)
    OPE[J1, J1] = OPE.make([
        k,    # 2nd order pole
        0     # 1st order pole (no structure constants for diagonal)
    ])

    # Define OPE J1(z)J2(w) (off-diagonal, regular)
    OPE[J1, J2] = OPEData()  # Regular (no poles)

    print(f"\nDefined OPE[J1, J1]:")
    result = OPE(J1, J1)
    print(f"  {result}")

    print(f"\nDefined OPE[J1, J2]:")
    result = OPE(J1, J2)
    print(f"  {result}")
    print(f"  Is zero: {result.is_zero()}")


def example_normal_ordering():
    """
    Example: Normal ordered products.
    """
    print("\n" + "=" * 60)
    print("Example: Normal Ordered Products")
    print("=" * 60)

    # Declare operators
    T = BasisOperator("T", bosonic=True)
    J = BasisOperator("J", bosonic=True)

    print(f"\nDeclared operators: {T}, {J}")

    # Create normal ordered products
    TJ = NO(T, J)
    print(f"\nNormal ordered product: {TJ}")
    print(f"  Left factor: {TJ.left}")
    print(f"  Right factor: {TJ.right}")

    # Test linearity of NO
    print(f"\nTesting NO linearity:")
    result = NO(T, J + 2 * J)
    print(f"  NO(T, J + 2J) = {result}")

    # Nested normal ordering
    TT = NO(T, T)
    TTJ = NO(TT, J)
    print(f"\nNested normal ordering:")
    print(f"  NO(T, T) = {TT}")
    print(f"  NO(NO(T, T), J) = {TTJ}")


def example_derivatives():
    """
    Example: Derivative operators.
    """
    print("\n" + "=" * 60)
    print("Example: Derivative Operators")
    print("=" * 60)

    # Declare operator
    T = BasisOperator("T", bosonic=True)

    print(f"\nDeclared operator: {T}")

    # Create derivatives
    dT = d(T)
    d2T = d(T, 2)
    d3T = d(T, 3)

    print(f"\nDerivatives:")
    print(f"  ∂T = {dT}")
    print(f"  ∂²T = {d2T}")
    print(f"  ∂³T = {d3T}")

    print(f"\nDerivative properties:")
    print(f"  Base of ∂T: {dT.base}")
    print(f"  Order of ∂²T: {d2T.order}")
    print(f"  Parity of ∂T: {dT.parity} (same as T)")


def example_ope_arithmetic():
    """
    Example: OPE arithmetic operations.
    """
    print("\n" + "=" * 60)
    print("Example: OPE Arithmetic")
    print("=" * 60)

    T = BasisOperator("T", bosonic=True)
    J = BasisOperator("J", bosonic=True)

    # Create some OPEs
    ope1 = OPEData({2: T, 1: J})
    ope2 = OPEData({2: T, 1: -J})

    print(f"\nOPE 1: {ope1}")
    print(f"OPE 2: {ope2}")

    # Addition
    result = ope1 + ope2
    print(f"\nOPE1 + OPE2 = {result}")
    print(f"  pole(2) = {result.pole(2)}")
    print(f"  pole(1) = {result.pole(1)}")

    # Scalar multiplication
    c = sp.Symbol("c")
    result = c * ope1
    print(f"\nc * OPE1 = {result}")
    print(f"  pole(2) = {result.pole(2)}")

    # Subtraction
    result = ope1 - ope2
    print(f"\nOPE1 - OPE2 = {result}")
    print(f"  pole(1) = {result.pole(1)}")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("PyOPE - Python Operator Product Expansion Library")
    print("Basic Functionality Examples")
    print("=" * 60)

    # Run examples
    example_virasoro()
    example_current_algebra()
    example_normal_ordering()
    example_derivatives()
    example_ope_arithmetic()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
