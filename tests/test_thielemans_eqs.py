"""
Test eq 3.3.1 and 3.3.2 from Thielemans paper.
"""

import sys
sys.path.insert(0, '/Users/lelouch/Nutstore Files/Math and Physics/NoteBooks/opepackage/pyope')

from pyope import (
    BasisOperator,
    OPE,
    d,
    One,
)
from pyope.ope_data import OPEData
import sympy as sp

print("=" * 60)
print("验证 Thielemans 论文 eq 3.3.1 和 3.3.2")
print("=" * 60)

# 创建算符
T = BasisOperator("T", )
c = sp.Symbol('c')

# 定义 Virasoro OPE: T(z)T(w)
OPE[T, T] = OPEData.from_list([
    c/2 * One,  # 4阶极点
    0,          # 3阶极点
    2*T,        # 2阶极点
    d(T)        # 1阶极点
])

print("\n定义的 OPE:")
print(f"OPE[T, T] = {OPE(T, T)}")

# 验证 eq 3.3.1: [∂T, T]_q = -(q-1)[T, T]_{q-1}
print("\n" + "=" * 60)
print("验证 eq 3.3.1: [∂T, T]_q = -(q-1)[T, T]_{q-1}")
print("=" * 60)

dT_T = OPE(d(T), T)
T_T = OPE(T, T)

print(f"\nOPE[∂T, T] = {dT_T}")
print(f"最高极点: {dT_T.max_pole}")

all_passed_1 = True
for q in range(1, dT_T.max_pole + 1):
    lhs = dT_T.pole(q)
    rhs = -(q - 1) * T_T.pole(q - 1) if q > 1 else sp.S.Zero
    diff = sp.expand(lhs - rhs)

    passed = (diff == 0)
    all_passed_1 = all_passed_1 and passed

    print(f"\nq = {q}:")
    print(f"  [∂T, T]_{q} = {lhs}")
    print(f"  -(q-1)[T, T]_{{{q-1}}} = {rhs}")
    print(f"  差值 = {diff}")
    print(f"  {'✓ 通过' if passed else '✗ 失败'}")

print(f"\neq 3.3.1 总体: {'✓ 全部通过' if all_passed_1 else '✗ 有失败'}")

# 验证 eq 3.3.2: [T, ∂T]_q = (q-1)[T, T]_{q-1} + ∂[T, T]_q
print("\n" + "=" * 60)
print("验证 eq 3.3.2: [T, ∂T]_q = (q-1)[T, T]_{q-1} + ∂[T, T]_q")
print("=" * 60)

T_dT = OPE(T, d(T))

print(f"\nOPE[T, ∂T] = {T_dT}")
print(f"最高极点: {T_dT.max_pole}")

all_passed_2 = True
for q in range(1, T_dT.max_pole + 1):
    lhs = T_dT.pole(q)

    term1 = (q - 1) * T_T.pole(q - 1) if q > 1 else sp.S.Zero
    term2 = d(T_T.pole(q)) if T_T.pole(q) != 0 else sp.S.Zero
    rhs = term1 + term2

    diff = sp.expand(lhs - rhs)

    passed = (diff == 0)
    all_passed_2 = all_passed_2 and passed

    print(f"\nq = {q}:")
    print(f"  [T, ∂T]_{q} = {lhs}")
    print(f"  (q-1)[T, T]_{{{q-1}}} = {term1}")
    print(f"  ∂[T, T]_{q} = {term2}")
    print(f"  右边总和 = {rhs}")
    print(f"  差值 = {diff}")
    print(f"  {'✓ 通过' if passed else '✗ 失败'}")

print(f"\neq 3.3.2 总体: {'✓ 全部通过' if all_passed_2 else '✗ 有失败'}")

print("\n" + "=" * 60)
if all_passed_1 and all_passed_2:
    print("🎉 所有验证通过！")
else:
    print("⚠️  有验证失败")
print("=" * 60)
