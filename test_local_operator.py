"""
Test script for LocalOperator architecture.

This script tests the new unified LocalOperator type system.
"""

import sys
sys.path.insert(0, '/Users/lelouch/Nutstore Files/Math and Physics/NoteBooks/opepackage/pyope')

from pyope import (
    BasisOperator,
    LocalOperator,
    OperatorSum,
    OperatorProduct,
    is_local_operator,
    OPE,
    NO,
    d,
    One,
)
from pyope.ope_data import OPEData
import sympy as sp

print("=" * 60)
print("Testing LocalOperator Architecture")
print("=" * 60)

# Test 1: Create basic operators
print("\n1. Creating basic operators:")
T = BasisOperator("T", bosonic=True)
J = BasisOperator("J", bosonic=True)
print(f"T = {T}, type = {type(T)}")
print(f"J = {J}, type = {type(J)}")
print(f"T is LocalOperator: {is_local_operator(T)}")
print(f"J is LocalOperator: {is_local_operator(J)}")

# Test 2: Arithmetic operations
print("\n2. Testing arithmetic operations:")
sum_op = T + J
print(f"T + J = {sum_op}, type = {type(sum_op)}")
print(f"(T + J) is LocalOperator: {is_local_operator(sum_op)}")

prod_op = 2 * T
print(f"2 * T = {prod_op}, type = {type(prod_op)}")
print(f"(2 * T) is LocalOperator: {is_local_operator(prod_op)}")

complex_op = 2 * T + 3 * J
print(f"2*T + 3*J = {complex_op}, type = {type(complex_op)}")
print(f"(2*T + 3*J) is LocalOperator: {is_local_operator(complex_op)}")

# Test 3: OPE with new types
print("\n3. Testing OPE with new types:")
c = sp.Symbol('c')

# Define T(z)T(w) OPE
OPE[T, T] = OPEData.from_list([
    c/2 * One,  # 4阶极点
    0,          # 3阶极点
    2*T,        # 2阶极点
    d(T)        # 1阶极点
])

print(f"OPE[T, T] = {OPE(T, T)}")

# Test OPE with sum
print(f"\nOPE[T, T+J]:")
ope_sum = OPE(T, T + J)
print(f"  Result = {ope_sum}")

# Test OPE with product
print(f"\nOPE[T, 2*T]:")
ope_prod = OPE(T, 2*T)
print(f"  Result = {ope_prod}")

# Test 4: NO with new types
print("\n4. Testing NO with new types:")
no_TJ = NO(T, J)
print(f"NO[T, J] = {no_TJ}, type = {type(no_TJ)}")

no_sum = NO(T, T + J)
print(f"NO[T, T+J] = {no_sum}")

no_prod = NO(2*T, J)
print(f"NO[2*T, J] = {no_prod}")

# Test 5: Parity
print("\n5. Testing parity:")
psi = BasisOperator("ψ", bosonic=False)
print(f"ψ parity = {psi.parity}, is_fermionic = {psi.is_fermionic}")

sum_mixed = T + psi
print(f"(T + ψ) parity = {sum_mixed.parity if hasattr(sum_mixed, 'parity') else 'N/A'}")

prod_fermion = 2 * psi
print(f"(2*ψ) parity = {prod_fermion.parity if hasattr(prod_fermion, 'parity') else 'N/A'}")

print("\n" + "=" * 60)
print("All tests completed!")
print("=" * 60)
