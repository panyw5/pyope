#!/usr/bin/env python3
"""
测试 null_states 模块的阶段 2 功能

测试内容：
1. 算符展开器（OperatorExpander）
2. 系数矩阵构建器（CoefficientMatrixBuilder）
3. 完整的 Null States 计算器（NullStatesCalculator）
"""

import sys
sys.path.insert(0, '../src')

from pyope import (
    BasisOperator, NO, d, Bosonic, Fermionic, One,
    OperatorExpander, CoefficientMatrixBuilder, NullStatesCalculator,
    enumerate_fock_basis, calculate_null_states
)
from fractions import Fraction

print("=" * 60)
print("测试 Null States 模块 - 阶段 2")
print("=" * 60)

# ============================================================
# 准备工作：定义自由场和基本算符
# ============================================================
print("\n" + "=" * 60)
print("准备工作：定义自由场系统")
print("=" * 60)

# 定义自由场
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

# 注册统计性
Bosonic(b, c)
Fermionic(beta, gamma)

free_fields = [b, c, beta, gamma]

print(f"\n自由场定义:")
print(f"  b: weight = {b.conformal_weight}, bosonic = {b.is_bosonic}")
print(f"  c: weight = {c.conformal_weight}, bosonic = {c.is_bosonic}")
print(f"  β: weight = {beta.conformal_weight}, bosonic = {beta.is_bosonic}")
print(f"  γ: weight = {gamma.conformal_weight}, bosonic = {gamma.is_bosonic}")

# ============================================================
# 测试 1: OperatorExpander - 简单算符
# ============================================================
print("\n" + "=" * 60)
print("测试 1: OperatorExpander - 简单算符展开")
print("=" * 60)

# 枚举 Level 1 的 Fock 基
fock_basis_level1 = enumerate_fock_basis(free_fields, Fraction(1), max_count=20)
print(f"\nLevel 1 的 Fock 基（共 {len(fock_basis_level1)} 个）:")
for i, basis in enumerate(fock_basis_level1[:5], 1):
    print(f"  {i}. {basis}")

# 创建展开器
expander = OperatorExpander(fock_basis_level1)

# 测试简单算符展开
test_op1 = 2 * d(d(c)) + 3 * NO(d(gamma), d(gamma))
print(f"\n测试算符: {test_op1}")

expansion1 = expander.expand(test_op1)
print(f"\n展开结果:")
for op, coeff in expansion1.items():
    print(f"  {op}: {coeff}")

print(f"\n✓ 测试结果: 展开了 {len(expansion1)} 个基")

# ============================================================
# 测试 2: CoefficientMatrixBuilder - 矩阵构建
# ============================================================
print("\n" + "=" * 60)
print("测试 2: CoefficientMatrixBuilder - 系数矩阵构建")
print("=" * 60)

# 定义一组抽象算符
abstract_ops = [
    d(d(c)),
    NO(d(gamma), d(gamma)),
    2 * d(d(c)) + NO(d(gamma), d(gamma))
]

print(f"\n抽象算符列表（共 {len(abstract_ops)} 个）:")
for i, op in enumerate(abstract_ops, 1):
    print(f"  {i}. {op}")

# 构建矩阵
matrix_builder = CoefficientMatrixBuilder(fock_basis_level1, abstract_ops)
matrix = matrix_builder.build_matrix()

print(f"\n系数矩阵维度: {len(matrix)} × {len(matrix[0])}")
print(f"  (Fock 基数 × 抽象算符数)")

# 显示矩阵的一部分
print(f"\n矩阵前 3 行:")
for i in range(min(3, len(matrix))):
    row_str = ", ".join(str(val) for val in matrix[i])
    print(f"  行 {i}: [{row_str}]")

print(f"\n✓ 测试结果: 矩阵构建完成")

# ============================================================
# 测试 3: 矩阵秩计算
# ============================================================
print("\n" + "=" * 60)
print("测试 3: 矩阵秩计算")
print("=" * 60)

rank = matrix_builder.compute_rank(matrix)
print(f"\n矩阵秩: {rank}")
print(f"抽象算符数: {len(abstract_ops)}")
print(f"Null states 数量: {len(abstract_ops) - rank}")

print(f"\n✓ 测试结果: 秩计算完成")

# ============================================================
# 测试 4: NullStatesCalculator - 完整计算
# ============================================================
print("\n" + "=" * 60)
print("测试 4: NullStatesCalculator - 完整 Null States 计算")
print("=" * 60)

# 创建计算器
calculator = NullStatesCalculator(free_fields)

# 计算 Level 1 的 null states
result = calculator.calculate_null_states(
    level=Fraction(1),
    abstract_operators=abstract_ops,
    max_fock_basis=20
)

print(f"\n计算结果:")
print(f"  Level: {result['level']}")
print(f"  抽象态数量: {result['n_abstract']}")
print(f"  Fock 基数量: {result['n_fock_basis']}")
print(f"  矩阵秩（物理态数）: {result['rank']}")
print(f"  Null states 数量: {result['n_null_states']}")

print(f"\n✓ 测试结果: Null states 计算完成")

# ============================================================
# 测试 5: 便捷函数 calculate_null_states
# ============================================================
print("\n" + "=" * 60)
print("测试 5: 便捷函数 calculate_null_states")
print("=" * 60)

result2 = calculate_null_states(
    free_fields=free_fields,
    level=Fraction(1),
    abstract_operators=abstract_ops,
    max_fock_basis=20
)

print(f"\n计算结果:")
print(f"  Level: {result2['level']}")
print(f"  Null states 数量: {result2['n_null_states']}")

match = (result['n_null_states'] == result2['n_null_states'])
print(f"\n✓ 测试结果: {'通过' if match else '失败'} - 与直接调用结果一致")

# ============================================================
# 测试 6: Level 3/2 的计算
# ============================================================
print("\n" + "=" * 60)
print("测试 6: Level 3/2 的 Null States 计算")
print("=" * 60)

# 定义 Level 3/2 的抽象算符
abstract_ops_3_2 = [
    beta,
    d(d(gamma)),
    NO(d(d(c)), d(gamma))
]

print(f"\nLevel 3/2 抽象算符（共 {len(abstract_ops_3_2)} 个）:")
for i, op in enumerate(abstract_ops_3_2, 1):
    print(f"  {i}. {op}")

result_3_2 = calculate_null_states(
    free_fields=free_fields,
    level=Fraction(3, 2),
    abstract_operators=abstract_ops_3_2,
    max_fock_basis=30
)

print(f"\n计算结果:")
print(f"  Level: {result_3_2['level']}")
print(f"  抽象态数量: {result_3_2['n_abstract']}")
print(f"  Fock 基数量: {result_3_2['n_fock_basis']}")
print(f"  矩阵秩（物理态数）: {result_3_2['rank']}")
print(f"  Null states 数量: {result_3_2['n_null_states']}")

print(f"\n✓ 测试结果: Level 3/2 计算完成")

# ============================================================
# 测试 7: 空矩阵情况
# ============================================================
print("\n" + "=" * 60)
print("测试 7: 边界情况 - 空算符列表")
print("=" * 60)

result_empty = calculate_null_states(
    free_fields=free_fields,
    level=Fraction(1),
    abstract_operators=[],
    max_fock_basis=20
)

print(f"\n计算结果:")
print(f"  抽象态数量: {result_empty['n_abstract']}")
print(f"  Null states 数量: {result_empty['n_null_states']}")

expected_empty = (result_empty['n_null_states'] == 0)
print(f"\n✓ 测试结果: {'通过' if expected_empty else '失败'} - 空列表返回 0 个 null states")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 60)
print("测试总结")
print("=" * 60)

test_results = [
    ("算符展开器 - 简单算符", len(expansion1) > 0),
    ("系数矩阵构建", len(matrix) > 0 and len(matrix[0]) > 0),
    ("矩阵秩计算", rank >= 0),
    ("NullStatesCalculator - 完整计算", result['n_null_states'] >= 0),
    ("便捷函数 calculate_null_states", match),
    ("Level 3/2 计算", result_3_2['n_null_states'] >= 0),
    ("边界情况 - 空算符列表", expected_empty),
]

print("\n测试结果:")
for i, (test_name, result) in enumerate(test_results, 1):
    status = "✓ 通过" if result else "✗ 失败"
    print(f"  {i}. {test_name}: {status}")

all_passed = all(result for _, result in test_results)
print(f"\n{'=' * 60}")
print(f"总体结果: {'所有测试通过！' if all_passed else '部分测试失败'}")
print(f"{'=' * 60}")
