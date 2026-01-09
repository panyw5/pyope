"""
检查所有 41 个算符的量子数分配

找出哪些算符被错误过滤
"""

from fractions import Fraction
from pyope.operators import BasisOperator, d
from pyope.api import NO
from pyope.null_states import (
    OperatorEnumerator,
    QuantumNumberCalculator,
    QuantumNumberGrouper
)

# 定义自由场
b = BasisOperator('b', 0, 2)
c = BasisOperator('c', 0, -1)
beta = BasisOperator('β', 0, Fraction(3, 2))
gamma = BasisOperator('γ', 0, Fraction(-1, 2))

# 定义 Z₃ W-algebra 生成元
w = beta
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
wb = NO(beta, NO(beta, NO(gamma, NO(gamma, gamma)))) + \
     2 * NO(beta, NO(gamma, NO(gamma, NO(b, c)))) - \
     4 * NO(beta, NO(d(gamma), gamma)) - \
     Fraction(4, 3) * NO(gamma, NO(b, d(c))) + \
     Fraction(2, 3) * NO(gamma, NO(d(b), c)) + \
     Fraction(2, 3) * NO(d(beta), NO(gamma, gamma)) - \
     Fraction(8, 3) * NO(d(gamma), NO(b, c)) + \
     Fraction(10, 9) * d(gamma, 2)

t = -2 * NO(b, d(c)) - Fraction(3, 2) * NO(beta, d(gamma)) - \
    NO(d(b), c) - Fraction(1, 2) * NO(d(beta), gamma)

g = NO(gamma, b)
gt = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))
gw = b
gwb = Fraction(8, 3) * NO(b, NO(d(c, 2), c)) + \
      3 * NO(beta, NO(beta, NO(gamma, NO(gamma, d(c))))) - \
      4 * NO(beta, NO(gamma, NO(b, NO(d(c), c)))) - \
      4 * NO(beta, NO(gamma, d(c, 2))) - \
      4 * NO(beta, NO(d(gamma), d(c))) - \
      Fraction(2, 3) * NO(d(b), NO(d(c), c)) + \
      2 * NO(d(beta), NO(beta, NO(gamma, NO(gamma, c)))) - \
      Fraction(8, 3) * NO(d(beta), NO(d(gamma), c)) + \
      Fraction(2, 3) * NO(d(beta, 2), NO(gamma, c)) + \
      Fraction(10, 9) * d(c, 3)

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

# 量子数映射
quantum_number_map = {
    b: (Fraction(1), Fraction(1, 2)),
    c: (Fraction(-1), Fraction(-1, 2)),
    beta: (Fraction(3, 2), Fraction(0)),
    gamma: (Fraction(-3, 2), Fraction(0)),
    w: (Fraction(3, 2), Fraction(0)),
    j0: (Fraction(0), Fraction(0)),
    wb: (Fraction(-3, 2), Fraction(0)),
    t: (Fraction(0), Fraction(0)),
    g: (Fraction(-1, 2), Fraction(1, 2)),
    gt: (Fraction(1, 2), Fraction(-1, 2)),
    gw: (Fraction(1), Fraction(1, 2)),
    gwb: (Fraction(-1), Fraction(-1, 2))
}

# 创建枚举器和量子数计算器
enumerator = OperatorEnumerator(generators)
quantum_calc = QuantumNumberCalculator(quantum_number_map)
grouper = QuantumNumberGrouper(quantum_calc)

# 枚举 level=4 的算符
print("=" * 80)
print("检查所有 41 个算符的量子数分配")
print("=" * 80)

ops = enumerator.enumerate_operators(Fraction(4))
print(f"\n总共枚举了 {len(ops)} 个算符\n")

# 检查每个算符的量子数
print(f"{'序号':>4s} {'量子数 (m, r)':>20s} {'算符':s}")
print("-" * 80)

for i, op in enumerate(ops, 1):
    qn = quantum_calc.get_quantum_numbers(op)
    op_str = str(op)
    if len(op_str) > 50:
        op_str = op_str[:47] + "..."
    print(f"{i:4d} ({str(qn[0]):>6s}, {str(qn[1]):>6s}) {op_str}")

# 按量子数分组（包括负 m）
print("\n" + "=" * 80)
print("按量子数分组（包括负 m 扇区）")
print("=" * 80)

operator_groups_all = grouper.group_operators(ops, only_non_negative_m=False)
print(f"\n{'量子数 (m, r)':>20s} {'算符数量':>10s}")
print("-" * 35)
for qn in sorted(operator_groups_all.keys()):
    print(f"{str(qn):>20s} {len(operator_groups_all[qn]):>10d}")

# 按量子数分组（只有 m≥0）
print("\n" + "=" * 80)
print("按量子数分组（只有 m≥0 扇区）")
print("=" * 80)

operator_groups_pos = grouper.group_operators(ops, only_non_negative_m=True)
print(f"\n{'量子数 (m, r)':>20s} {'算符数量':>10s}")
print("-" * 35)
for qn in sorted(operator_groups_pos.keys()):
    print(f"{str(qn):>20s} {len(operator_groups_pos[qn]):>10d}")

# 找出被过滤掉的算符
print("\n" + "=" * 80)
print("被过滤掉的算符（m < 0）")
print("=" * 80)

filtered_count = 0
for i, op in enumerate(ops, 1):
    qn = quantum_calc.get_quantum_numbers(op)
    if qn[0] < 0:
        filtered_count += 1
        op_str = str(op)
        if len(op_str) > 50:
            op_str = op_str[:47] + "..."
        print(f"{i:4d} ({str(qn[0]):>6s}, {str(qn[1]):>6s}) {op_str}")

print(f"\n总共 {filtered_count} 个算符被过滤（m < 0）")
print(f"保留 {len(ops) - filtered_count} 个算符（m ≥ 0）")
