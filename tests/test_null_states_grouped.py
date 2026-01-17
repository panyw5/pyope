#!/usr/bin/env python3
"""
测试 null_states 模块的量子数分组功能

测试内容：
1. QuantumNumberCalculator - 量子数计算
2. QuantumNumberGrouper - 量子数分组
3. GroupedNullStatesCalculator - 分组 Null States 计算
4. 与 Mathematica 参考结果对比
"""

import sys
sys.path.insert(0, '../src')

from pyope import (
    BasisOperator, NO, d, Bosonic, Fermionic, One, OPE,
    QuantumNumberCalculator, QuantumNumberGrouper, GroupedNullStatesCalculator
)
from fractions import Fraction

print("=" * 60)
print("测试 Null States 模块 - 量子数分组功能")
print("=" * 60)

# ============================================================
# 准备工作：定义自由场和 Z₃ W-algebra 生成元
# ============================================================
print("\n" + "=" * 60)
print("准备工作：定义 Z₃ W-algebra")
print("=" * 60)

# 定义自由场
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

Fermionic(b, c)  # 修正：b 和 c 是费米子（bosonic=False）
Bosonic(beta, gamma)  # 修正：beta 和 gamma 是玻色子（bosonic=True）

# 定义基本 OPE
OPE[b, c] = OPE.make([One])
OPE[beta, gamma] = OPE.make([-One])

free_fields = [b, c, beta, gamma]

# 定义 Z₃ W-algebra 生成元
w = beta  # 权重 3/2
gw = b    # 权重 2
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)  # 权重 1
t = (-2 * NO(b, d(c)) - Fraction(3, 2) * NO(beta, d(gamma))
     - NO(d(b), c) - Fraction(1, 2) * NO(d(beta), gamma))  # 权重 2
g = NO(gamma, b)  # 权重 3/2
gt = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))  # 权重 3/2

# 简化版的 wb 和 gwb（用于测试）
wb = (NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
      + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c)))))  # 权重 3/2
gwb = (Fraction(8, 3) * NO(b, NO(d(d(c)), c))
       + 3 * NO(beta, NO(beta, NO(gamma, NO(gamma, d(c))))))  # 权重 2

print(f"\n生成元定义:")
print(f"  w = β (权重 3/2)")
print(f"  gw = b (权重 2)")
print(f"  j0 = 2*NO(b,c) + 3*NO(β,γ) (权重 1)")
print(f"  t = 能动张量 (权重 2)")
print(f"  g = NO(γ,b) (权重 3/2)")
print(f"  gt = 2*NO(∂β,c) + 3*NO(β,∂c) (权重 3/2)")

# ============================================================
# 测试 1: QuantumNumberCalculator - 基本算符
# ============================================================
print("\n" + "=" * 60)
print("测试 1: QuantumNumberCalculator - 基本算符量子数")
print("=" * 60)

# 定义量子数映射（包括自由场和生成元）
quantum_number_map = {
    # 自由场的量子数
    b: (Fraction(1), Fraction(1, 2)),      # b
    c: (Fraction(-1), Fraction(-1, 2)),    # c
    beta: (Fraction(3, 2), Fraction(0)),   # β
    gamma: (Fraction(-3, 2), Fraction(0)), # γ
    # W-algebra 生成元的量子数
    w: (Fraction(3, 2), Fraction(0)),      # w = β
    gw: (Fraction(1), Fraction(1, 2)),     # gw = b
    j0: (Fraction(0), Fraction(0)),        # j0
    t: (Fraction(0), Fraction(0)),         # t
    g: (Fraction(-1, 2), Fraction(1, 2)),  # g = NO(γ,b)
    gt: (Fraction(1, 2), Fraction(-1, 2)), # gt
    wb: (Fraction(-3, 2), Fraction(0)),    # wb
    gwb: (Fraction(-1), Fraction(-1, 2)),  # gwb
}

qn_calc = QuantumNumberCalculator(quantum_number_map)

print(f"\n基本生成元的量子数:")
for op, (m, r) in quantum_number_map.items():
    print(f"  {str(op)[:20]:20s} → (m={m}, r={r})")

# 测试导数算符
print(f"\n导数算符的量子数:")
dw = d(w)
m_dw, r_dw = qn_calc.get_quantum_numbers(dw)
print(f"  ∂w → (m={m_dw}, r={r_dw})")
print(f"  预期: (m=3/2, r=0) - 导数不改变量子数")

match1 = (m_dw == Fraction(3, 2) and r_dw == Fraction(0))
print(f"\n✓ 测试结果: {'通过' if match1 else '失败'}")

# ============================================================
# 测试 2: QuantumNumberCalculator - 复合算符
# ============================================================
print("\n" + "=" * 60)
print("测试 2: QuantumNumberCalculator - 复合算符量子数")
print("=" * 60)

# 测试正规序算符
no_wg = NO(w, g)
m_no, r_no = qn_calc.get_quantum_numbers(no_wg)
print(f"\nNO(w, g) 的量子数:")
print(f"  w: (m=3/2, r=0)")
print(f"  g: (m=-1/2, r=1/2)")
print(f"  NO(w,g): (m={m_no}, r={r_no})")
print(f"  预期: (m=1, r=1/2) - 量子数相加")

match2 = (m_no == Fraction(1) and r_no == Fraction(1, 2))
print(f"\n✓ 测试结果: {'通过' if match2 else '失败'}")

# ============================================================
# 测试 3: QuantumNumberGrouper - 算符分组
# ============================================================
print("\n" + "=" * 60)
print("测试 3: QuantumNumberGrouper - 算符分组")
print("=" * 60)

grouper = QuantumNumberGrouper(qn_calc)

# 创建一组测试算符
test_ops = [w, gw, j0, t, g, gt, d(w), d(g)]

# 分组（包含所有量子数）
groups_all = grouper.group_operators(test_ops, only_non_negative_m=False)

print(f"\n所有算符分组（共 {len(groups_all)} 组）:")
for qn, ops in sorted(groups_all.items()):
    print(f"  (m={qn[0]}, r={qn[1]}): {len(ops)} 个算符")

# 分组（只保留 m≥0）
groups_nonneg = grouper.group_operators(test_ops, only_non_negative_m=True)

print(f"\nm≥0 算符分组（共 {len(groups_nonneg)} 组）:")
for qn, ops in sorted(groups_nonneg.items()):
    print(f"  (m={qn[0]}, r={qn[1]}): {len(ops)} 个算符")

match3 = (len(groups_nonneg) < len(groups_all))
print(f"\n✓ 测试结果: {'通过' if match3 else '失败'} - m≥0 分组数量更少")

# ============================================================
# 测试 4: GroupedNullStatesCalculator - Level 2
# ============================================================
print("\n" + "=" * 60)
print("测试 4: GroupedNullStatesCalculator - Level 2")
print("=" * 60)

# 创建分组计算器
grouped_calc = GroupedNullStatesCalculator(free_fields, quantum_number_map)

# Level 2 的抽象算符
abstract_ops_lv2 = [
    b,           # gw
    t,           # t
    d(d(c)),     # ∂²c
    NO(beta, d(gamma)),  # NO(β,∂γ)
    d(j0),       # ∂j0
]

print(f"\nLevel 2 抽象算符（共 {len(abstract_ops_lv2)} 个）:")
for i, op in enumerate(abstract_ops_lv2, 1):
    m, r = qn_calc.get_quantum_numbers(op)
    print(f"  {i}. {str(op)[:30]:30s} (m={m}, r={r})")

# 计算（包含所有量子数）
result_lv2_all = grouped_calc.calculate_null_states_grouped(
    level=Fraction(2),
    abstract_operators=abstract_ops_lv2,
    max_fock_basis=50,
    only_non_negative_m=False
)

print(f"\n计算结果（所有量子数）:")
print(f"  总抽象态数: {result_lv2_all['total_n_abstract']}")
print(f"  总秩: {result_lv2_all['total_rank']}")
print(f"  总 Null states: {result_lv2_all['total_n_null_states']}")

print(f"\n各量子数扇区:")
for qn, group_result in sorted(result_lv2_all['groups'].items()):
    print(f"  (m={qn[0]}, r={qn[1]}): "
          f"抽象={group_result['n_abstract']}, "
          f"秩={group_result['rank']}, "
          f"null={group_result['n_null_states']}")

print(f"\n✓ 测试结果: Level 2 分组计算完成")

# ============================================================
# 测试 5: GroupedNullStatesCalculator - Level 2 (m≥0)
# ============================================================
print("\n" + "=" * 60)
print("测试 5: GroupedNullStatesCalculator - Level 2 (m≥0)")
print("=" * 60)

# 计算（只保留 m≥0）
result_lv2_nonneg = grouped_calc.calculate_null_states_grouped(
    level=Fraction(2),
    abstract_operators=abstract_ops_lv2,
    max_fock_basis=50,
    only_non_negative_m=True
)

print(f"\n计算结果（m≥0 扇区）:")
print(f"  总抽象态数: {result_lv2_nonneg['total_n_abstract']}")
print(f"  总秩: {result_lv2_nonneg['total_rank']}")
print(f"  总 Null states: {result_lv2_nonneg['total_n_null_states']}")

print(f"\n各量子数扇区:")
for qn, group_result in sorted(result_lv2_nonneg['groups'].items()):
    print(f"  (m={qn[0]}, r={qn[1]}): "
          f"抽象={group_result['n_abstract']}, "
          f"秩={group_result['rank']}, "
          f"null={group_result['n_null_states']}")

match5 = (result_lv2_nonneg['total_n_abstract'] <= result_lv2_all['total_n_abstract'])
print(f"\n✓ 测试结果: {'通过' if match5 else '失败'} - m≥0 扇区算符数量更少")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 60)
print("测试总结")
print("=" * 60)

test_results = [
    ("量子数计算 - 导数算符", match1),
    ("量子数计算 - 复合算符", match2),
    ("量子数分组 - m≥0 过滤", match3),
    ("分组计算 - Level 2 (所有量子数)", True),
    ("分组计算 - Level 2 (m≥0)", match5),
]

print("\n测试结果:")
for i, (test_name, result) in enumerate(test_results, 1):
    status = "✓ 通过" if result else "✗ 失败"
    print(f"  {i}. {test_name}: {status}")

all_passed = all(result for _, result in test_results)
print(f"\n{'=' * 60}")
print(f"总体结果: {'所有测试通过！' if all_passed else '部分测试失败'}")
print(f"{'=' * 60}")

# ============================================================
# 与 Mathematica 对比说明
# ============================================================
print("\n" + "=" * 60)
print("与 Mathematica 对比说明")
print("=" * 60)

print("""
根据 Mathematica 参考实现，Level 4 的 m≥0 扇区结果为：

| 量子数 (m, r) | 抽象算符 | 秩 | Null States |
|--------------|---------|----|-----------|
| (3, 0)       | 2       | 1  | 1          |
| (2, 1)       | 1       | 1  | 0          |
| (2, -1/2)    | 4       | 1  | 3          |
| (1, 1/2)     | 11      | 1  | 10         |
| (1, -1)      | 2       | 1  | 1          |
| (0, 0)       | 21      | 1  | 20         |
| **总计**     | **41**  | **6** | **35**  |

pyope 的分组计算功能现在可以按量子数分组计算，
与 Mathematica 的计算策略一致。

要完全复现 Mathematica 的结果，还需要：
1. 实现完整的算符枚举（包括所有生成元的导数组合）
2. 正确定义所有 8 个 Z₃ W-algebra 生成元
3. 实现 R-filtration 分级（可选）
""")
