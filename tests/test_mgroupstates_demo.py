#!/usr/bin/env python3
"""
演示 mgroupstates 对应的 Python 实现

对比 Mathematica 和 Python 的分组逻辑：
- Mathematica: mgroupstates -> mgfulllist -> rtset -> fullgrouping
- Python: group_by_m_only() -> group_by_r_within_m()
"""

import sys
sys.path.insert(0, '../src')

from fractions import Fraction
from pyope import (
    BasisOperator, NO, d,
    QuantumNumberCalculator, QuantumNumberGrouper
)

print("=" * 70)
print("演示：mgroupstates 的 Python 实现")
print("=" * 70)

# ============================================================
# 第 1 部分：定义自由场和生成元
# ============================================================
print("\n第 1 部分：定义自由场和生成元")

b = BasisOperator('b', weight=Fraction(1), statistics='fermionic')
c = BasisOperator('c', weight=Fraction(0), statistics='fermionic')
beta = BasisOperator('β', weight=Fraction(1, 2), statistics='bosonic')
gamma = BasisOperator('γ', weight=Fraction(1, 2), statistics='bosonic')

# 定义 Z₃ W-algebra 生成元（简化版）
w = beta  # 权重 3/2
gw = b    # 权重 2
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)  # 权重 1
g = NO(gamma, b)  # 权重 3/2

# ============================================================
# 第 2 部分：定义量子数映射
# ============================================================
print("\n第 2 部分：定义量子数映射")
print("  m: osp(2|2) gl(1) 荷")
print("  r: 外自同构 gl(1) 荷")

quantum_number_map = {
    w: (Fraction(3, 2), Fraction(0)),      # m=3/2, r=0
    gw: (Fraction(1), Fraction(1, 2)),     # m=1, r=1/2
    j0: (Fraction(0), Fraction(0)),        # m=0, r=0
    g: (Fraction(-1, 2), Fraction(1, 2)),  # m=-1/2, r=1/2
}

# ============================================================
# 第 3 部分：创建测试算符
# ============================================================
print("\n第 3 部分：创建测试算符")

# 单个生成元及其导数
test_ops = [
    w,
    d(w),
    gw,
    j0,
    d(j0),
    g,
    d(g),
    # 正规序乘积
    NO(w, w),      # m=3, r=0
    NO(w, g),      # m=1, r=1/2
    NO(g, g),      # m=-1, r=1
]

print(f"  总共 {len(test_ops)} 个测试算符")

# ============================================================
# 第 4 部分：创建分组器并演示 mgroupstates
# ============================================================
print("\n" + "=" * 70)
print("第 4 部分：演示 group_by_m_only() (对应 mgroupstates)")
print("=" * 70)

calculator = QuantumNumberCalculator(quantum_number_map)
grouper = QuantumNumberGrouper(calculator)

# 方法 1: 只按 m 分组（对应 mgroupstates）
states_by_m = grouper.group_by_m_only(test_ops)

print(f"\n按 m 量子数分组（共 {len(states_by_m)} 个扇区）:")
for m, ops in states_by_m.items():
    print(f"\nm = {m}:")
    for i, op in enumerate(ops, 1):
        m_val, r_val = calculator.get_quantum_numbers(op)
        print(f"  [{i}] {op} (m={m_val}, r={r_val})")

# ============================================================
# 第 5 部分：演示 fullgrouping
# ============================================================
print("\n" + "=" * 70)
print("第 5 部分：演示 group_by_r_within_m() (对应 fullgrouping)")
print("=" * 70)

# 方法 2: 在 m 扇区内按 r 分组（对应 fullgrouping）
states_by_m_r = grouper.group_by_r_within_m(states_by_m)

print(f"\n按 (m, r) 扇区分组（共 {len(states_by_m_r)} 个扇区）:")
for (m, r), ops in states_by_m_r.items():
    print(f"\n(m={m}, r={r}): {len(ops)} 个算符")
    for i, op in enumerate(ops, 1):
        print(f"  [{i}] {op}")

# ============================================================
# 第 6 部分：演示 only_non_negative_m（对应 sfullgrouping）
# ============================================================
print("\n" + "=" * 70)
print("第 6 部分：演示 only_non_negative_m=True (对应 sfullgrouping)")
print("=" * 70)

states_by_m_nonneg = grouper.group_by_m_only(test_ops, only_non_negative_m=True)

print(f"\n只保留 m≥0 扇区（共 {len(states_by_m_nonneg)} 个扇区）:")
for m, ops in states_by_m_nonneg.items():
    print(f"  m={m}: {len(ops)} 个算符")

# ============================================================
# 第 7 部分：获取 mparam
# ============================================================
print("\n" + "=" * 70)
print("第 7 部分：演示 get_m_values() (对应 mparam)")
print("=" * 70)

m_values = grouper.get_m_values(test_ops)
m_values_nonneg = grouper.get_m_values(test_ops, only_non_negative_m=True)

print(f"\n所有 m 值: {m_values}")
print(f"非负 m 值: {m_values_nonneg}")

# ============================================================
# 第 8 部分：对比两种方法
# ============================================================
print("\n" + "=" * 70)
print("第 8 部分：对比两种分组方法")
print("=" * 70)

# 方法 A: 两步法（对应 Mathematica）
states_by_m = grouper.group_by_m_only(test_ops)
states_by_m_r_twostep = grouper.group_by_r_within_m(states_by_m)

# 方法 B: 直接法（现有实现）
states_by_m_r_direct = grouper.group_operators(test_ops)

print("\n验证两种方法结果一致:")
print(f"  两步法扇区数: {len(states_by_m_r_twostep)}")
print(f"  直接法扇区数: {len(states_by_m_r_direct)}")

# 验证键集合相同
assert set(states_by_m_r_twostep.keys()) == set(states_by_m_r_direct.keys()), \
    "两种方法的扇区键不一致！"

# 验证每个扇区的算符列表相同
for key in states_by_m_r_twostep.keys():
    ops_twostep = set(str(op) for op in states_by_m_r_twostep[key])
    ops_direct = set(str(op) for op in states_by_m_r_direct[key])
    assert ops_twostep == ops_direct, f"扇区 {key} 的算符不一致！"

print("  ✓ 验证通过：两种方法结果完全一致")

print("\n" + "=" * 70)
print("演示完成")
print("=" * 70)
print("\n总结:")
print("  - group_by_m_only() 对应 Mathematica 的 mgroupstates")
print("  - group_by_r_within_m() 对应 Mathematica 的 fullgrouping")
print("  - 两步法与直接法 group_operators() 结果一致")
print("  - Python 实现使用显式字典，避免了隐式索引对齐")
