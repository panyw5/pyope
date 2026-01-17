#!/usr/bin/env python
"""
测试 notebook 中的导数展开用例
"""

from pyope import BasisOperator, d, NO, simplify, Bosonic
from fractions import Fraction

print("=" * 70)
print("测试 Notebook 中的导数展开用例")
print("=" * 70)

# 定义算符
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

Bosonic(b, c, beta, gamma)

print("\n1. 测试 d(NO(b, c))")
print("-" * 70)
expr1 = d(NO(b, c))
print(f"原始表达式: {expr1}")
result1 = simplify(expr1)
print(f"simplify() 结果: {result1}")
print(f"期望: NO(d(b), c) + NO(b, d(c))")

print("\n2. 测试能动张量的导数")
print("-" * 70)
t = -2 * NO(b, d(c)) - Fraction(3, 2) * NO(beta, d(gamma)) - NO(d(b), c) - Fraction(1, 2) * NO(d(beta), gamma)
print(f"t = {t}")

dt = d(t)
print(f"\nd(t) = {dt}")

dt_simplified = simplify(dt)
print(f"\nsimplify(d(t)) = {dt_simplified}")

print("\n3. 测试禁用展开")
print("-" * 70)
expr3 = d(NO(b, c))
result3_disabled = simplify(expr3, expand_derivatives=False)
print(f"simplify(d(NO(b,c)), expand_derivatives=False) = {result3_disabled}")
print(f"期望: d(NO(b, c)) (保持未展开)")

print("\n4. 测试 OPEData.simplify() 方法")
print("-" * 70)
from pyope import OPE, MakeOPE, One

# 定义简单的 OPE
OPE[b, c] = MakeOPE([One])

# 计算包含导数的 OPE
ope_result = OPE(d(NO(b, c)), b)
print(f"OPE(d(NO(b,c)), b) = {ope_result}")

# 使用 .simplify() 方法
simplified_ope = ope_result.simplify()
print(f"ope_result.simplify() = {simplified_ope}")

# 使用 expand_derivatives=False
simplified_ope_no_expand = ope_result.simplify(expand_derivatives=False)
print(f"ope_result.simplify(expand_derivatives=False) = {simplified_ope_no_expand}")

print("\n" + "=" * 70)
print("测试完成！")
print("=" * 70)
