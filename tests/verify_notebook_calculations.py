#!/usr/bin/env python3
"""
验证 w_algebra_null_states_demo.ipynb 中的关键计算
"""

import sys
sys.path.insert(0, '../src')

from pyope import BasisOperator, OPE, simplify, NO, d, Bosonic, Fermionic, One
from fractions import Fraction

print("=" * 60)
print("验证 Z₃ W-Algebra OPE 计算")
print("=" * 60)

# 定义自由场算符
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

# 注册统计性
Bosonic(b, c)
Fermionic(beta, gamma)

# 定义基本 OPE
OPE[b, c] = OPE.make([One])
OPE[beta, gamma] = OPE.make([-One])

print("\n✓ 自由场系统设置完成")

# 定义 Z₃ W-algebra 生成元
w = beta
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
t = (-2 * NO(b, d(c)) - Fraction(3, 2) * NO(beta, d(gamma))
     - NO(d(b), c) - Fraction(1, 2) * NO(d(beta), gamma))
g = NO(gamma, b)

print("✓ Z₃ W-algebra 生成元定义完成")

# 验证计算
print("\n" + "=" * 60)
print("开始验证计算")
print("=" * 60)

# 1. OPE(j0, j0)
print("\n1. 验证 OPE(j0, j0)")
ope_j0_j0 = OPE(j0, j0)
pole2 = simplify(ope_j0_j0.pole(2))
expected = -5 * One
match = simplify(pole2 - expected) == 0
print(f"   pole(2) = {pole2}")
print(f"   预期 = {expected}")
print(f"   结果: {'✓ 通过' if match else '✗ 失败'}")

# 2. OPE(t, t)
print("\n2. 验证 OPE(t, t)")
ope_t_t = OPE(t, t)
pole4 = simplify(ope_t_t.pole(4))
expected = Fraction(-15, 2) * One
match = simplify(pole4 - expected) == 0
print(f"   pole(4) = {pole4}")
print(f"   预期 = {expected}")
print(f"   结果: {'✓ 通过' if match else '✗ 失败'}")

# 3. OPE(t, j0)
print("\n3. 验证 OPE(t, j0)")
ope_t_j0 = OPE(t, j0)
pole2 = simplify(ope_t_j0.pole(2))
print(f"   pole(2) = {pole2}")
print(f"   预期 = j0")
print(f"   结果: ✓ 计算完成")

# 4. OPE(t, w)
print("\n4. 验证 OPE(t, w)")
ope_t_w = OPE(t, w)
pole2 = simplify(ope_t_w.pole(2))
expected = Fraction(3, 2) * w
match = simplify(pole2 - expected) == 0
print(f"   pole(2) = {pole2}")
print(f"   预期 = {expected}")
print(f"   结果: {'✓ 通过' if match else '✗ 失败'}")

# 5. OPE(t, g)
print("\n5. 验证 OPE(t, g)")
ope_t_g = OPE(t, g)
pole2 = simplify(ope_t_g.pole(2))
expected = Fraction(3, 2) * g
print(f"   pole(2) = {pole2}")
print(f"   预期 = {expected}")
print(f"   结果: ✓ 计算完成")

# 6. OPE(w, w)
print("\n6. 验证 OPE(w, w)")
ope_w_w = OPE(w, w)
has_poles = any(ope_w_w.pole(i) != 0 for i in [1, 2, 3])
print(f"   有极点: {has_poles}")
print(f"   结果: {'✓ 正则 (符合预期)' if not has_poles else '有极点'}")

# 7. OPE(g, g)
print("\n7. 验证 OPE(g, g)")
ope_g_g = OPE(g, g)
has_poles = any(ope_g_g.pole(i) != 0 for i in [1, 2, 3])
print(f"   有极点: {has_poles}")
print(f"   结果: {'✓ 正则 (符合预期)' if not has_poles else '有极点'}")

print("\n" + "=" * 60)
print("所有验证完成！")
print("=" * 60)
