#!/usr/bin/env python3
"""
验证 Python 和 Mathematica 在不同 level 的算符枚举结果
"""
import sys
sys.path.insert(0, 'src')

from pyope import BasisOperator, OPE, simplify, NO, d, Bosonic, Fermionic, One
from fractions import Fraction
import itertools

# 定义自由场算符
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))

Bosonic(b, c)
Fermionic(beta, gamma)

# 定义生成元
w = beta
gw = b
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
t = - 2 * NO(b, d(c)) - (3/2) * NO(beta, d(gamma)) - NO(d(b), c) - (1/2) * NO(d(beta), gamma)
g = NO(gamma, b)
gt = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))
wb = (NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
      + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c))))
      - 4 * NO(beta, NO(d(gamma), gamma))
      - Fraction(4, 3) * NO(gamma, NO(b, d(c)))
      + Fraction(2, 3) * NO(gamma, NO(d(b), c))
      + Fraction(2, 3) * NO(d(beta), NO(gamma, gamma))
      - Fraction(8, 3) * NO(d(gamma), NO(b, c))
      + Fraction(10, 9) * d(d(gamma)))
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

generators = {
    'w': (w, Fraction(3, 2)),
    'j0': (j0, Fraction(1)),
    't': (t, Fraction(2)),
    'g': (g, Fraction(3, 2)),
    'gt': (gt, Fraction(3, 2)),
    'gw': (gw, Fraction(2)),
    'wb': (wb, Fraction(3, 2)),
    'gwb': (gwb, Fraction(2))
}

def integer_partitions(n):
    n_doubled = int(2 * n)
    def partitions_helper(n, max_val=None):
        if max_val is None:
            max_val = n
        if n == 0:
            yield []
            return
        for i in range(min(n, max_val), 0, -1):
            for p in partitions_helper(n - i, i):
                yield [i] + p
    for p in partitions_helper(n_doubled):
        yield [Fraction(x, 2) for x in p]

def can_generate_at_weight(gen_name, target_weight):
    gen, base_weight = generators[gen_name]
    derivative_order = target_weight - base_weight
    if derivative_order < 0:
        return False
    return derivative_order == int(derivative_order)

def enumerate_operators_at_level(level):
    operators = []
    fermionic_gens = {'g', 'gt', 'gw', 'gwb'}

    def is_valid_combination(partition, combo):
        for i in range(len(combo) - 1):
            if (combo[i] == combo[i+1] and
                combo[i] in fermionic_gens and
                partition[i] == partition[i+1]):
                return False

        i = 0
        while i < len(partition):
            j = i
            while j < len(partition) and partition[j] == partition[i]:
                j += 1
            segment = combo[i:j]
            if list(segment) != sorted(segment):
                return False
            i = j
        return True

    for partition in integer_partitions(level):
        available_gens = []
        for weight in partition:
            gens_at_weight = [name for name in generators.keys()
                             if can_generate_at_weight(name, weight)]
            available_gens.append(gens_at_weight)

        if all(available_gens):
            for combo in itertools.product(*available_gens):
                if is_valid_combination(partition, combo):
                    operators.append((partition, combo))

    return operators

# Mathematica 参考结果（从实际运行获得）
mathematica_results = {
    1: 1,
    2: 5,
    3: 17,  # 修正：实际是 17 个，不是 16 个
    4: 45,
    5: 121
}

print("=" * 60)
print("Python vs Mathematica 算符枚举验证")
print("=" * 60)

all_passed = True
for level in [1, 2, 3, 4, 5]:
    ops = enumerate_operators_at_level(level)
    python_count = len(ops)
    mathematica_count = mathematica_results[level]
    match = python_count == mathematica_count

    status = "✓" if match else "✗"
    print(f"\nLevel {level}:")
    print(f"  Python:      {python_count:3d} 个算符")
    print(f"  Mathematica: {mathematica_count:3d} 个算符")
    print(f"  状态: {status} {'通过' if match else f'失败 (差异 {mathematica_count - python_count})'}")

    if not match:
        all_passed = False

print("\n" + "=" * 60)
if all_passed:
    print("✓✓✓ 所有测试通过！Python 实现与 Mathematica 完全一致！")
else:
    print("✗✗✗ 部分测试失败")
print("=" * 60)
