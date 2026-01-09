#!/usr/bin/env python3
"""
Level 4 Null States 验证测试

与 Mathematica 参考实现进行完整对比
"""

import sys
sys.path.insert(0, '../src')

from pyope import (
    BasisOperator, NO, d, Bosonic, Fermionic, One, OPE,
    GroupedNullStatesCalculator, OperatorEnumerator
)
from fractions import Fraction

print("=" * 70)
print("Level 4 Null States 验证测试")
print("=" * 70)

# ============================================================
# 第 1 部分：定义自由场系统
# ============================================================
print("\n" + "=" * 70)
print("第 1 部分：定义自由场系统")
print("=" * 70)

# 定义自由场
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

Bosonic(b, c)
Fermionic(beta, gamma)

# 定义基本 OPE
OPE[b, c] = OPE.make([One])
OPE[beta, gamma] = OPE.make([-One])

free_fields = [b, c, beta, gamma]

print(f"\n自由场定义:")
print(f"  bc 系统: b (权重 2), c (权重 -1)")
print(f"  βγ 系统: β (权重 3/2), γ (权重 -1/2)")
print(f"\n基本 OPE:")
print(f"  OPE[b, c] = 1/(z-w)")
print(f"  OPE[β, γ] = -1/(z-w)")

# ============================================================
# 第 2 部分：定义 Z₃ W-algebra 生成元
# ============================================================
print("\n" + "=" * 70)
print("第 2 部分：定义 Z₃ W-algebra 生成元")
print("=" * 70)

# 简单生成元
w = beta  # 权重 3/2
gw = b    # 权重 2

# U(1) 流
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)  # 权重 1

# 能动张量
t = (-2 * NO(b, d(c)) - Fraction(3, 2) * NO(beta, d(gamma))
     - NO(d(b), c) - Fraction(1, 2) * NO(d(beta), gamma))  # 权重 2

# 费米生成元
g = NO(gamma, b)  # 权重 3/2
gt = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))  # 权重 3/2

# 复杂 W-场 wb（权重 3/2）
print("\n构造复杂生成元 wb 和 gwb...")
wb = (NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
      + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c))))
      - 4 * NO(beta, NO(d(gamma), gamma))
      - Fraction(4, 3) * NO(gamma, NO(b, d(c)))
      + Fraction(2, 3) * NO(gamma, NO(d(b), c))
      + Fraction(2, 3) * NO(d(beta), NO(gamma, gamma))
      - Fraction(8, 3) * NO(d(gamma), NO(b, c))
      + Fraction(10, 9) * d(d(gamma)))

# 复杂费米场 gwb（权重 2）
gwb = (Fraction(8, 3) * NO(b, NO(d(d(c)), c))
       + 3 * NO(beta, NO(beta, NO(gamma, NO(gamma, d(c)))))
       - 4 * NO(beta, NO(gamma, NO(b, NO(d(c), c))))
       - 4 * NO(beta, NO(gamma, d(d(c))))
       - 4 * NO(beta, NO(d(gamma), d(c)))
       - Fraction(2, 3) * NO(d(b), NO(d(c), c))
       + 2 * NO(d(beta), NO(beta, NO(gamma, NO(gamma, c))))
       - Fraction(8, 3) * NO(d(beta), NO(d(gamma), c))
       + Fraction(2, 3) * NO(d(d(beta)), NO(gamma, c))
       + Fraction(10, 9) * d(d(d(c))))

print(f"✓ 生成元定义完成")

# 生成元字典
generators = {
    'w': {'op': w, 'weight': Fraction(3, 2), 'm': Fraction(3, 2), 'r': Fraction(0)},
    'j0': {'op': j0, 'weight': Fraction(1), 'm': Fraction(0), 'r': Fraction(0)},
    'wb': {'op': wb, 'weight': Fraction(3, 2), 'm': Fraction(-3, 2), 'r': Fraction(0)},
    't': {'op': t, 'weight': Fraction(2), 'm': Fraction(0), 'r': Fraction(0)},
    'g': {'op': g, 'weight': Fraction(3, 2), 'm': Fraction(-1, 2), 'r': Fraction(1, 2)},
    'gt': {'op': gt, 'weight': Fraction(3, 2), 'm': Fraction(1, 2), 'r': Fraction(-1, 2)},
    'gw': {'op': gw, 'weight': Fraction(2), 'm': Fraction(1), 'r': Fraction(1, 2)},
    'gwb': {'op': gwb, 'weight': Fraction(2), 'm': Fraction(-1), 'r': Fraction(-1, 2)}
}

print(f"\n生成元列表（共 {len(generators)} 个）:")
for name, info in generators.items():
    print(f"  {name:4s}: 权重={info['weight']}, (m={info['m']}, r={info['r']})")

# ============================================================
# 第 3 部分：定义量子数映射
# ============================================================
print("\n" + "=" * 70)
print("第 3 部分：定义量子数映射")
print("=" * 70)

quantum_number_map = {
    # 自由场
    b: (Fraction(1), Fraction(1, 2)),
    c: (Fraction(-1), Fraction(-1, 2)),
    beta: (Fraction(3, 2), Fraction(0)),
    gamma: (Fraction(-3, 2), Fraction(0)),
    # W-algebra 生成元
    w: (Fraction(3, 2), Fraction(0)),
    j0: (Fraction(0), Fraction(0)),
    wb: (Fraction(-3, 2), Fraction(0)),
    t: (Fraction(0), Fraction(0)),
    g: (Fraction(-1, 2), Fraction(1, 2)),
    gt: (Fraction(1, 2), Fraction(-1, 2)),
    gw: (Fraction(1), Fraction(1, 2)),
    gwb: (Fraction(-1), Fraction(-1, 2))
}

print(f"✓ 量子数映射定义完成（共 {len(quantum_number_map)} 个算符）")

# ============================================================
# 第 4 部分：枚举 Level 4 的抽象算符
# ============================================================
print("\n" + "=" * 70)
print("第 4 部分：枚举 Level 4 的抽象算符")
print("=" * 70)

enumerator = OperatorEnumerator(generators)
abstract_ops_lv4 = enumerator.enumerate_operators(Fraction(4), max_derivative_order=10)

print(f"\nLevel 4 抽象算符（共 {len(abstract_ops_lv4)} 个）:")
for i, op in enumerate(abstract_ops_lv4[:10], 1):
    print(f"  {i}. {str(op)[:50]}")
if len(abstract_ops_lv4) > 10:
    print(f"  ... 还有 {len(abstract_ops_lv4) - 10} 个算符")

# ============================================================
# 第 5 部分：创建分组计算器并计算 Level 4
# ============================================================
print("\n" + "=" * 70)
print("第 5 部分：Level 4 Null States 计算")
print("=" * 70)

grouped_calc = GroupedNullStatesCalculator(free_fields, quantum_number_map)

print(f"\n开始计算 Level 4 (m≥0 扇区)...")
result_lv4 = grouped_calc.calculate_null_states_grouped(
    level=Fraction(4),
    abstract_operators=abstract_ops_lv4,
    max_fock_basis=200,
    only_non_negative_m=True,
    filter_linearly_independent=False  # 禁用过滤，使用所有枚举的算符
)

print(f"✓ 计算完成")

# ============================================================
# 第 6 部分：输出结果
# ============================================================
print("\n" + "=" * 70)
print("第 6 部分：计算结果")
print("=" * 70)

print(f"\n总体结果:")
print(f"  总抽象态数: {result_lv4['total_n_abstract']}")
print(f"  总秩: {result_lv4['total_rank']}")
print(f"  总 Null states: {result_lv4['total_n_null_states']}")

print(f"\n各量子数扇区详细结果:")
print(f"{'量子数 (m, r)':20s} {'抽象':>6s} {'秩':>6s} {'Null':>6s}")
print("-" * 45)

for qn in sorted(result_lv4['groups'].keys()):
    group = result_lv4['groups'][qn]
    qn_str = f"({qn[0]}, {qn[1]})"
    print(f"{qn_str:20s} {group['n_abstract']:6d} {group['rank']:6d} {group['n_null_states']:6d}")

# ============================================================
# 第 7 部分：与 Mathematica 参考结果对比
# ============================================================
print("\n" + "=" * 70)
print("第 7 部分：与 Mathematica 参考结果对比")
print("=" * 70)

# Mathematica Level 4 参考结果（m≥0 扇区）
mathematica_results = {
    (Fraction(3), Fraction(0)): {'n_abstract': 2, 'rank': 1, 'n_null': 1},
    (Fraction(2), Fraction(1)): {'n_abstract': 1, 'rank': 1, 'n_null': 0},
    (Fraction(2), Fraction(-1, 2)): {'n_abstract': 4, 'rank': 1, 'n_null': 3},
    (Fraction(1), Fraction(1, 2)): {'n_abstract': 11, 'rank': 1, 'n_null': 10},
    (Fraction(1), Fraction(-1)): {'n_abstract': 2, 'rank': 1, 'n_null': 1},
    (Fraction(0), Fraction(0)): {'n_abstract': 21, 'rank': 1, 'n_null': 20}
}

mathematica_total_abstract = sum(v['n_abstract'] for v in mathematica_results.values())
mathematica_total_rank = sum(v['rank'] for v in mathematica_results.values())
mathematica_total_null = sum(v['n_null'] for v in mathematica_results.values())

print(f"\nMathematica 参考结果:")
print(f"  总抽象态数: {mathematica_total_abstract}")
print(f"  总秩: {mathematica_total_rank}")
print(f"  总 Null states: {mathematica_total_null}")

print(f"\n对比分析:")
print(f"{'量子数 (m, r)':20s} {'pyope':>15s} {'Mathematica':>15s} {'状态':>10s}")
print("-" * 70)

all_quantum_numbers = set(result_lv4['groups'].keys()) | set(mathematica_results.keys())

for qn in sorted(all_quantum_numbers):
    qn_str = f"({qn[0]}, {qn[1]})"

    if qn in result_lv4['groups']:
        pyope_str = f"{result_lv4['groups'][qn]['n_abstract']}/{result_lv4['groups'][qn]['rank']}/{result_lv4['groups'][qn]['n_null_states']}"
    else:
        pyope_str = "N/A"

    if qn in mathematica_results:
        math_str = f"{mathematica_results[qn]['n_abstract']}/{mathematica_results[qn]['rank']}/{mathematica_results[qn]['n_null']}"
    else:
        math_str = "N/A"

    # 判断是否匹配
    if qn in result_lv4['groups'] and qn in mathematica_results:
        pyope_data = result_lv4['groups'][qn]
        math_data = mathematica_results[qn]
        match = (pyope_data['n_abstract'] == math_data['n_abstract'] and
                 pyope_data['rank'] == math_data['rank'] and
                 pyope_data['n_null_states'] == math_data['n_null'])
        status = "✓ 匹配" if match else "✗ 不匹配"
    else:
        status = "? 缺失"

    print(f"{qn_str:20s} {pyope_str:>15s} {math_str:>15s} {status:>10s}")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 70)
print("总结")
print("=" * 70)

print(f"\npyope 计算结果:")
print(f"  总抽象态数: {result_lv4['total_n_abstract']}")
print(f"  总秩: {result_lv4['total_rank']}")
print(f"  总 Null states: {result_lv4['total_n_null_states']}")

print(f"\nMathematica 参考结果:")
print(f"  总抽象态数: {mathematica_total_abstract}")
print(f"  总秩: {mathematica_total_rank}")
print(f"  总 Null states: {mathematica_total_null}")

print(f"\n差异分析:")
print(f"  抽象态数差异: {result_lv4['total_n_abstract'] - mathematica_total_abstract}")
print(f"  秩差异: {result_lv4['total_rank'] - mathematica_total_rank}")
print(f"  Null states 差异: {result_lv4['total_n_null_states'] - mathematica_total_null}")

print(f"\n{'=' * 70}")
print(f"Level 4 验证测试完成")
print(f"{'=' * 70}")
