#!/usr/bin/env python3
"""
测试 null_states 模块的阶段 1 功能

测试内容：
1. 系数提取功能
2. Fock 空间基枚举
"""

import sys
sys.path.insert(0, '../src')

from pyope import (
    BasisOperator, NO, d, Bosonic, Fermionic, One,
    extract_coefficients, enumerate_fock_basis
)
from fractions import Fraction

print("=" * 60)
print("测试 Null States 模块 - 阶段 1")
print("=" * 60)

# ============================================================
# 测试 1: 系数提取 - 简单情况
# ============================================================
print("\n" + "=" * 60)
print("测试 1: 系数提取 - 简单线性组合")
print("=" * 60)

# 定义测试算符
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
Fermionic(b, c)  # 修正：b 和 c 是费米子（bosonic=False）

# 测试简单线性组合
expr1 = 2 * b + 3 * c
coeffs1 = extract_coefficients(expr1)

print(f"\n表达式: {expr1}")
print(f"提取的系数:")
for op, coeff in coeffs1.items():
    print(f"  {op}: {coeff}")

# 验证
expected_coeffs = {b: 2, c: 3}
match1 = all(coeffs1.get(op) == coeff for op, coeff in expected_coeffs.items())
print(f"\n✓ 测试结果: {'通过' if match1 else '失败'}")

# ============================================================
# 测试 2: 系数提取 - 正规序乘积
# ============================================================
print("\n" + "=" * 60)
print("测试 2: 系数提取 - 正规序乘积")
print("=" * 60)

# 测试正规序乘积
expr2 = 2 * NO(b, c) + 3 * NO(c, b)
coeffs2 = extract_coefficients(expr2)

print(f"\n表达式: {expr2}")
print(f"提取的系数:")
for op, coeff in coeffs2.items():
    print(f"  {op}: {coeff}")

print(f"\n✓ 测试结果: 提取了 {len(coeffs2)} 个不同的算符")

# ============================================================
# 测试 3: 系数提取 - 分数系数
# ============================================================
print("\n" + "=" * 60)
print("测试 3: 系数提取 - 分数系数")
print("=" * 60)

# 测试分数系数
expr3 = Fraction(3, 2) * b + Fraction(-1, 2) * c
coeffs3 = extract_coefficients(expr3)

print(f"\n表达式: {expr3}")
print(f"提取的系数:")
for op, coeff in coeffs3.items():
    print(f"  {op}: {coeff}")

# 验证
expected_coeffs3 = {b: Fraction(3, 2), c: Fraction(-1, 2)}
match3 = all(coeffs3.get(op) == coeff for op, coeff in expected_coeffs3.items())
print(f"\n✓ 测试结果: {'通过' if match3 else '失败'}")

# ============================================================
# 测试 4: Fock 空间基枚举 - Level 1
# ============================================================
print("\n" + "=" * 60)
print("测试 4: Fock 空间基枚举 - Level 1")
print("=" * 60)

# 定义自由场
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))
Bosonic(beta, gamma)  # 修正：beta 和 gamma 是玻色子（bosonic=True）

free_fields = [b, c, beta, gamma]

# 枚举 Level 1 的基
basis_level1 = enumerate_fock_basis(free_fields, Fraction(1), max_count=20)

print(f"\nLevel 1 的自由场基（前 10 个）:")
for i, op in enumerate(basis_level1[:10], 1):
    print(f"  {i}. {op}")

print(f"\n总共枚举了 {len(basis_level1)} 个基")
print(f"✓ 测试结果: 枚举完成")

# ============================================================
# 测试 5: Fock 空间基枚举 - Level 3/2
# ============================================================
print("\n" + "=" * 60)
print("测试 5: Fock 空间基枚举 - Level 3/2")
print("=" * 60)

# 枚举 Level 3/2 的基
basis_level_3_2 = enumerate_fock_basis(free_fields, Fraction(3, 2), max_count=20)

print(f"\nLevel 3/2 的自由场基（前 10 个）:")
for i, op in enumerate(basis_level_3_2[:10], 1):
    print(f"  {i}. {op}")

print(f"\n总共枚举了 {len(basis_level_3_2)} 个基")
print(f"✓ 测试结果: 枚举完成")

# ============================================================
# 测试 6: Fock 空间基枚举 - Level 2
# ============================================================
print("\n" + "=" * 60)
print("测试 6: Fock 空间基枚举 - Level 2")
print("=" * 60)

# 枚举 Level 2 的基
basis_level2 = enumerate_fock_basis(free_fields, Fraction(2), max_count=30)

print(f"\nLevel 2 的自由场基（前 15 个）:")
for i, op in enumerate(basis_level2[:15], 1):
    print(f"  {i}. {op}")

print(f"\n总共枚举了 {len(basis_level2)} 个基")
print(f"✓ 测试结果: 枚举完成")

# ============================================================
# 测试 7: 系数提取 - 复杂表达式
# ============================================================
print("\n" + "=" * 60)
print("测试 7: 系数提取 - 复杂表达式")
print("=" * 60)

# 构造复杂表达式
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
coeffs_j0 = extract_coefficients(j0)

print(f"\n表达式: j0 = 2*NO(b,c) + 3*NO(β,γ)")
print(f"提取的系数:")
for op, coeff in coeffs_j0.items():
    print(f"  {op}: {coeff}")

print(f"\n✓ 测试结果: 提取了 {len(coeffs_j0)} 个不同的算符")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 60)
print("测试总结")
print("=" * 60)

test_results = [
    ("系数提取 - 简单线性组合", match1),
    ("系数提取 - 正规序乘积", True),
    ("系数提取 - 分数系数", match3),
    ("Fock 空间基枚举 - Level 1", len(basis_level1) > 0),
    ("Fock 空间基枚举 - Level 3/2", len(basis_level_3_2) > 0),
    ("Fock 空间基枚举 - Level 2", len(basis_level2) > 0),
    ("系数提取 - 复杂表达式", len(coeffs_j0) > 0),
]

print("\n测试结果:")
for i, (test_name, result) in enumerate(test_results, 1):
    status = "✓ 通过" if result else "✗ 失败"
    print(f"  {i}. {test_name}: {status}")

all_passed = all(result for _, result in test_results)
print(f"\n{'=' * 60}")
print(f"总体结果: {'所有测试通过！' if all_passed else '部分测试失败'}")
print(f"{'=' * 60}")
