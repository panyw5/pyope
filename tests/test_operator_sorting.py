#!/usr/bin/env python3
"""
测试算符排序功能

验证算符枚举时的排序是否与 Mathematica 的 opfields 顺序一致
"""

import sys
sys.path.insert(0, '../src')

from pyope import BasisOperator, d, Bosonic, Fermionic, ope_registry
from pyope.null_states import OperatorEnumerator
from fractions import Fraction

print("=" * 80)
print("算符排序测试")
print("=" * 80)

# ============================================================
# 第 1 部分：定义 Z₃ W-algebra 生成元（按 Mathematica opfields 顺序）
# ============================================================
print("\n" + "=" * 80)
print("第 1 部分：定义生成元（按 Mathematica opfields 顺序）")
print("=" * 80)

# 按照 Mathematica 的 opfields 顺序定义算符
# opfields={ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]}

w = BasisOperator('w', bosonic=True, conformal_weight=Fraction(3, 2))
j0 = BasisOperator('j0', bosonic=True, conformal_weight=Fraction(1))
wb = BasisOperator('wb', bosonic=True, conformal_weight=Fraction(3, 2))
t = BasisOperator('t', bosonic=True, conformal_weight=Fraction(2))
g = BasisOperator('g', bosonic=False, conformal_weight=Fraction(3, 2))
gt = BasisOperator('gt', bosonic=False, conformal_weight=Fraction(3, 2))
gw = BasisOperator('gw', bosonic=False, conformal_weight=Fraction(2))
gwb = BasisOperator('gwb', bosonic=False, conformal_weight=Fraction(2))

# 按照 opfields 顺序注册算符
# 这个顺序很关键！必须与 Mathematica 一致
Bosonic(w, j0, wb, t)
Fermionic(g, gt, gw, gwb)

print("\n算符注册顺序（对应 Mathematica opfields）:")
print("  1. w   (weight=3/2, bosonic)")
print("  2. j0  (weight=1, bosonic)")
print("  3. wb  (weight=3/2, bosonic)")
print("  4. t   (weight=2, bosonic)")
print("  5. g   (weight=3/2, fermionic)")
print("  6. gt  (weight=3/2, fermionic)")
print("  7. gw  (weight=2, fermionic)")
print("  8. gwb (weight=2, fermionic)")

# 验证注册顺序
print("\n验证注册位置:")
for op in [w, j0, wb, t, g, gt, gw, gwb]:
    pos = ope_registry.get_position(op)
    print(f"  {op.name}: position = {pos}")

# ============================================================
# 第 2 部分：测试 ope_registry.compare_operators()
# ============================================================
print("\n" + "=" * 80)
print("第 2 部分：测试 ope_registry.compare_operators()")
print("=" * 80)

# 测试基本算符比较
print("\n测试 1: 基本算符比较")
test_pairs = [
    (w, j0, "w < j0"),
    (j0, wb, "j0 < wb"),
    (wb, g, "wb < g"),
    (g, gt, "g < gt"),
]

for left, right, expected in test_pairs:
    result = ope_registry.compare_operators(left, right)
    status = "✓" if result > 0 else "✗"
    print(f"  {status} compare({left.name}, {right.name}) = {result:2d}  (期望: {expected})")

# 测试导数算符比较
print("\n测试 2: 导数算符比较")
test_pairs_deriv = [
    (w, d(w), "w < ∂w"),
    (d(w), d(d(w)), "∂w < ∂²w"),
    (d(j0), wb, "∂j0 < wb (因为 j0 < wb)"),
    (w, d(j0), "w < ∂j0 (因为 w < j0)"),
]

for left, right, expected in test_pairs_deriv:
    result = ope_registry.compare_operators(left, right)
    status = "✓" if result > 0 else "✗"
    left_str = str(left).replace('d(', '∂').replace(')', '')
    right_str = str(right).replace('d(', '∂').replace(')', '')
    print(f"  {status} compare({left_str}, {right_str}) = {result:2d}  (期望: {expected})")

# ============================================================
# 第 3 部分：测试 OperatorEnumerator
# ============================================================
print("\n" + "=" * 80)
print("第 3 部分：测试 OperatorEnumerator 的排序")
print("=" * 80)

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

enumerator = OperatorEnumerator(generators)

# 测试权重 3/2 的算符排序
print("\n测试 3: 权重 3/2 的算符排序")
ops_3_2 = enumerator._generate_operators_at_weight(Fraction(3, 2), max_derivative_order=5)
print(f"权重 3/2 的算符（共 {len(ops_3_2)} 个）:")
for i, op in enumerate(ops_3_2, 1):
    op_str = str(op).replace('d(', '∂').replace(')', '')
    print(f"  {i}. {op_str}")

print("\n期望顺序（按 opfields）:")
print("  1. w   (基础权重 3/2)")
print("  2. wb  (基础权重 3/2)")
print("  3. g   (基础权重 3/2)")
print("  4. gt  (基础权重 3/2)")
print("  注: j0 的基础权重是 1，∂j0 的权重是 2，不在这里")

# 验证顺序
expected_order = [w, wb, g, gt]
match = all(str(ops_3_2[i]) == str(expected_order[i]) for i in range(len(expected_order)))
print(f"\n排序验证: {'✓ 通过' if match else '✗ 失败'}")

# 测试权重 5/2 的算符排序
print("\n测试 4: 权重 5/2 的算符排序")
ops_5_2 = enumerator._generate_operators_at_weight(Fraction(5, 2), max_derivative_order=5)
print(f"权重 5/2 的算符（共 {len(ops_5_2)} 个）:")
for i, op in enumerate(ops_5_2, 1):
    op_str = str(op).replace('d(', '∂').replace(')', '').replace('d^2', '∂²').replace('d^3', '∂³')
    print(f"  {i}. {op_str}")

print("\n期望顺序（按 opfields）:")
print("  1. ∂w   (3/2 + 1 = 5/2)")
print("  2. ∂wb  (3/2 + 1 = 5/2)")
print("  3. ∂g   (3/2 + 1 = 5/2)")
print("  4. ∂gt  (3/2 + 1 = 5/2)")
print("  注: ∂²j0 的权重是 1+2=3，∂t 的权重是 2+1=3，都不是 5/2")

# 验证顺序
expected_order_5_2 = [d(w), d(wb), d(g), d(gt)]
match_5_2 = all(str(ops_5_2[i]) == str(expected_order_5_2[i]) for i in range(len(expected_order_5_2)))
print(f"\n排序验证: {'✓ 通过' if match_5_2 else '✗ 失败'}")

# ============================================================
# 第 4 部分：测试完整枚举
# ============================================================
print("\n" + "=" * 80)
print("第 4 部分：测试完整算符枚举")
print("=" * 80)

# 枚举 level 2 的算符
print("\n测试 5: Level 2 的算符枚举")
ops_level_2 = enumerator.enumerate_operators(Fraction(2), max_derivative_order=5)
print(f"Level 2 的算符（共 {len(ops_level_2)} 个）:")
for i, op in enumerate(ops_level_2[:20], 1):  # 只显示前 20 个
    op_str = str(op).replace('d(', '∂').replace(')', '').replace('d^2', '∂²')
    print(f"  {i}. {op_str}")

if len(ops_level_2) > 20:
    print(f"  ... (还有 {len(ops_level_2) - 20} 个)")

# ============================================================
# 第 5 部分：与 Mathematica 结果对比
# ============================================================
print("\n" + "=" * 80)
print("第 5 部分：与 Mathematica 结果对比")
print("=" * 80)

print("\nMathematica 的示例输出（来自 finalgen）:")
print("  NO[w[1,1,z],g[0,1,z]]  -> NO(∂w, g)")
print("  NO[w[0,1,z],g[1,1,z]]  -> NO(w, ∂g)")
print("  NO[g[0,1,z],g[1,1,z]]  -> NO(g, ∂g)")
print("  NO[j0[0,1,z],j0[2,1,z]] -> NO(j0, ∂²j0)")

print("\n这些来自不同的分拆:")
print("  NO(∂w, g)   来自分拆 [5/2, 3/2]")
print("  NO(w, ∂g)   来自分拆 [3/2, 3/2]")
print("  NO(g, ∂g)   来自分拆 [3/2, 3/2]")
print("  NO(j0, ∂²j0) 来自分拆 [1, 3]")

print("\n关键点:")
print("  1. 分拆按降序排列（权重大的在前）")
print("  2. 同一权重位置的算符按 opfields 顺序枚举")
print("  3. 导数不改变基本算符的相对顺序")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("测试总结")
print("=" * 80)

all_passed = match and match_5_2
print(f"\n总体结果: {'✓ 所有测试通过' if all_passed else '✗ 部分测试失败'}")

if all_passed:
    print("\n算符排序实现正确！")
    print("- 使用 ope_registry.compare_operators() 进行排序")
    print("- 排序顺序与 Mathematica 的 opfields 一致")
    print("- 导数算符正确继承基本算符的顺序")
else:
    print("\n需要检查:")
    print("- 算符注册顺序是否正确")
    print("- compare_operators() 的实现是否正确")
    print("- _generate_operators_at_weight() 的排序逻辑")

print("\n" + "=" * 80)
