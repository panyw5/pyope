#!/usr/bin/env python3
"""
对比不同方法枚举 Level 4 算符的数量
"""

import sys
sys.path.insert(0, '../src')

from pyope import BasisOperator, NO, d, Bosonic, Fermionic
from pyope.null_states import OperatorEnumerator, integer_partitions
from fractions import Fraction
import itertools

print("=" * 80)
print("Level 4 算符枚举对比")
print("=" * 80)

# 定义算符（按 Mathematica opfields 顺序）
w = BasisOperator('w', bosonic=True, conformal_weight=Fraction(3, 2))
j0 = BasisOperator('j0', bosonic=True, conformal_weight=Fraction(1))
wb = BasisOperator('wb', bosonic=True, conformal_weight=Fraction(3, 2))
t = BasisOperator('t', bosonic=True, conformal_weight=Fraction(2))
g = BasisOperator('g', bosonic=False, conformal_weight=Fraction(3, 2))
gt = BasisOperator('gt', bosonic=False, conformal_weight=Fraction(3, 2))
gw = BasisOperator('gw', bosonic=False, conformal_weight=Fraction(2))
gwb = BasisOperator('gwb', bosonic=False, conformal_weight=Fraction(2))

# 注册
Bosonic(w, j0, wb, t)
Fermionic(g, gt, gw, gwb)

# 创建生成元字典
generators = {
    'w': {'op': w, 'weight': Fraction(3, 2)},
    'j0': {'op': j0, 'weight': Fraction(1)},
    'wb': {'op': wb, 'weight': Fraction(3, 2)},
    't': {'op': t, 'weight': Fraction(2)},
    'g': {'op': g, 'weight': Fraction(3, 2)},
    'gt': {'op': gt, 'weight': Fraction(3, 2)},
    'gw': {'op': gw, 'weight': Fraction(2)},
    'gwb': {'op': gwb, 'weight': Fraction(2)}
}

# ============================================================
# 方法 1: 使用 OperatorEnumerator (use_partition_method=True)
# ============================================================
print("\n" + "=" * 80)
print("方法 1: OperatorEnumerator (use_partition_method=True)")
print("=" * 80)

enumerator1 = OperatorEnumerator(generators)
ops1 = enumerator1.enumerate_operators(Fraction(4), max_derivative_order=10, use_partition_method=True)
print(f"\n生成算符数量: {len(ops1)}")
print(f"\n前 10 个算符:")
for i, op in enumerate(ops1[:10], 1):
    print(f"  {i:3d}. {op}")

# ============================================================
# 方法 2: 使用 OperatorEnumerator (use_partition_method=False)
# ============================================================
print("\n" + "=" * 80)
print("方法 2: OperatorEnumerator (use_partition_method=False)")
print("=" * 80)

enumerator2 = OperatorEnumerator(generators)
ops2 = enumerator2.enumerate_operators(Fraction(4), max_derivative_order=10, use_partition_method=False)
print(f"\n生成算符数量: {len(ops2)}")
print(f"\n前 10 个算符:")
for i, op in enumerate(ops2[:10], 1):
    print(f"  {i:3d}. {op}")

# ============================================================
# 方法 3: 手动枚举（不去重）
# ============================================================
print("\n" + "=" * 80)
print("方法 3: 手动枚举（不去重）")
print("=" * 80)

def manual_enumerate(level, generators, max_deriv=10):
    """手动枚举，不去重"""
    operators = []

    # 生成所有整数分拆
    partitions = list(integer_partitions(level))
    print(f"\nLevel {level} 的分拆数量: {len(partitions)}")
    print(f"前 10 个分拆:")
    for i, p in enumerate(partitions[:10], 1):
        print(f"  {i}. {p}")

    for partition in partitions:
        # 对于每个分拆，生成所有可能的排列
        unique_permutations = set()
        for perm in itertools.permutations(partition):
            unique_permutations.add(perm)

        for perm in unique_permutations:
            # 为排列中的每个权重生成算符列表
            operator_lists = []
            for weight in perm:
                ops_at_weight = []
                for name, gen_info in generators.items():
                    base_weight = gen_info['weight']
                    deriv_order = weight - base_weight

                    if deriv_order >= 0 and deriv_order <= max_deriv:
                        if deriv_order.denominator == 1:
                            deriv_order_int = int(deriv_order)
                            if deriv_order_int == 0:
                                ops_at_weight.append(gen_info['op'])
                            else:
                                ops_at_weight.append(d(gen_info['op'], deriv_order_int))

                if ops_at_weight:
                    operator_lists.append(ops_at_weight)

            # 如果所有权重都有对应的算符，生成笛卡尔积
            if len(operator_lists) == len(perm):
                for combination in itertools.product(*operator_lists):
                    if len(combination) == 1:
                        operators.append(combination[0])
                    else:
                        # 构造嵌套 NO
                        result = combination[-1]
                        for op in reversed(combination[:-1]):
                            result = NO(op, result)
                        operators.append(result)

    return operators

ops3 = manual_enumerate(Fraction(4), generators)
print(f"\n生成算符数量（去重前）: {len(ops3)}")

# 去重
unique_ops3 = []
seen = set()
for op in ops3:
    op_str = str(op)
    if op_str not in seen:
        seen.add(op_str)
        unique_ops3.append(op)

print(f"生成算符数量（去重后）: {len(unique_ops3)}")

# ============================================================
# 对比分析
# ============================================================
print("\n" + "=" * 80)
print("对比分析")
print("=" * 80)

print(f"\n方法 1 (partition=True):  {len(ops1)} 个算符")
print(f"方法 2 (partition=False): {len(ops2)} 个算符")
print(f"方法 3 (手动枚举):        {len(unique_ops3)} 个算符")

# 检查差异
set1 = set(str(op) for op in ops1)
set2 = set(str(op) for op in ops2)
set3 = set(str(op) for op in unique_ops3)

print(f"\n集合差异:")
print(f"  方法1 - 方法2: {len(set1 - set2)} 个")
print(f"  方法2 - 方法1: {len(set2 - set1)} 个")
print(f"  方法1 - 方法3: {len(set1 - set3)} 个")
print(f"  方法3 - 方法1: {len(set3 - set1)} 个")

if set1 - set3:
    print(f"\n方法1 有但方法3 没有的算符:")
    for op_str in sorted(set1 - set3)[:10]:
        print(f"  {op_str}")

if set3 - set1:
    print(f"\n方法3 有但方法1 没有的算符:")
    for op_str in sorted(set3 - set1)[:10]:
        print(f"  {op_str}")

print("\n" + "=" * 80)
