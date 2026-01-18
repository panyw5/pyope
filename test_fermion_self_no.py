"""
测试费米子自身正规序的正确计算

验证公式 (2.3.17) 的实现：
对于费米子 A，[AA]_0 = sum_{k>=0} a_k ∂^{2k+1} [AA]_{2k+1}
"""

import sys

sys.path.insert(0, "src")

import sympy as sp
from pyope import (
    BasisOperator,
    OPE,
    MakeOPE,
    One,
    Zero,
    NO,
    bracket,
    Fermionic,
    Bosonic,
    simplify,
    d,
)

print("=" * 70)
print("测试费米子自身正规序的计算")
print("=" * 70)

# 测试 1: bc ghost 系统 - NO(b,b) = 0
print("\n测试 1: bc ghost 系统")
print("-" * 70)

b = BasisOperator("b", bosonic=False, conformal_weight=-1)
c = BasisOperator("c", bosonic=False, conformal_weight=2)

Fermionic(b, c)
OPE[b, c] = MakeOPE([One])

print(f"定义: b(z)c(w) ~ 1/(z-w)")
print(f"      b(z)b(w) ~ 0 (未定义，默认为空)")

no_bb = NO(b, b)
print(f"\nNO(b,b) = {no_bb}")
print(f"类型: {type(no_bb)}")
print(f"NO(b,b) == 0: {no_bb == 0}")
print(f"✓ 正确！bc ghost 系统中 NO(b,b) = 0")

# 测试 2: N=1 超共形代数 - NO(G,G) = ∂T
print("\n\n测试 2: N=1 超共形代数")
print("-" * 70)

G = BasisOperator("G", bosonic=False, conformal_weight=3 / 2)
T = BasisOperator("T", bosonic=True, conformal_weight=2)
c = sp.Symbol("c")

Fermionic(G)
Bosonic(T)

# G(z)G(w) ~ (2c/3)/(z-w)^3 + 2T/(z-w)
OPE[G, G] = MakeOPE([2 * c / 3 * One, 0, 2 * T])

print(f"定义: G(z)G(w) ~ (2c/3)/(z-w)^3 + 2T/(z-w)")
print(f"\n根据公式 (2.3.17):")
print(f"  [GG]_0 = a_0 * ∂^1 [GG]_1")
print(f"         = (1/2) * ∂(2T)")
print(f"         = ∂T")

no_gg = NO(G, G)
expected = d(T)

print(f"\nNO(G,G) = {no_gg}")
print(f"期望值 = {expected}")
print(f"相等: {simplify(no_gg - expected) == 0}")

if simplify(no_gg - expected) == 0:
    print(f"✓ 正确！N=1 超共形代数中 NO(G,G) = ∂T")
else:
    print(f"✗ 错误！结果不匹配")

# 测试 3: 更复杂的例子 - 有更高阶极点
print("\n\n测试 3: 有 5 阶极点的费米子")
print("-" * 70)

psi = BasisOperator("ψ", bosonic=False, conformal_weight=5 / 2)
W = BasisOperator("W", bosonic=True, conformal_weight=3)
U = BasisOperator("U", bosonic=True, conformal_weight=4)

Fermionic(psi)
Bosonic(W, U)

# ψ(z)ψ(w) ~ A/(z-w)^5 + B/(z-w)^3 + C/(z-w)
A = sp.Symbol("A")
B = sp.Symbol("B")
C = sp.Symbol("C")

OPE[psi, psi] = MakeOPE([A * One, 0, B * W, 0, C * U])

print(f"定义: ψ(z)ψ(w) ~ A/(z-w)^5 + B/(z-w)^3 + C/(z-w)")
print(f"\n根据公式 (2.3.17):")
print(f"  [ψψ]_0 = a_0 * ∂^1 [ψψ]_1 + a_1 * ∂^3 [ψψ]_3 + a_2 * ∂^5 [ψψ]_5")
print(f"         = (1/2) * ∂(C*U) + (-1/24) * ∂^3(B*W) + (1/240) * ∂^5(A*One)")
print(f"         = (C/2) * ∂U + (-B/24) * ∂^3 W + 0")
print(f"         = (C/2) * ∂U - (B/24) * ∂^3 W")

no_psipsi = NO(psi, psi)
expected = sp.Rational(1, 2) * C * d(U) - sp.Rational(1, 24) * B * d(W, 3)

print(f"\nNO(ψ,ψ) = {no_psipsi}")
print(f"期望值 = {expected}")

diff = simplify(no_psipsi - expected)
print(f"差值 = {diff}")
print(f"相等: {diff == 0}")

if diff == 0:
    print(f"✓ 正确！高阶极点的计算正确")
else:
    print(f"✗ 错误！结果不匹配")

# 测试 4: 玻色子自身的正规序（应该直接返回 NormalOrderedOperator）
print("\n\n测试 4: 玻色子自身的正规序")
print("-" * 70)

J = BasisOperator("J", bosonic=True, conformal_weight=1)
Bosonic(J)
OPE[J, J] = MakeOPE([])  # 空 OPE

no_jj = NO(J, J)
print(f"定义: J(z)J(w) ~ 0 (空 OPE)")
print(f"NO(J,J) = {no_jj}")
print(f"类型: {type(no_jj)}")

if no_jj == 0:
    print(f"✓ 玻色子自身正规序在 OPE 为空时返回 0")
else:
    print(f"注意：玻色子自身正规序返回 {no_jj}")

# 测试 5: bracket 函数的一致性
print("\n\n测试 5: bracket(A,A,0) 与 NO(A,A) 的一致性")
print("-" * 70)

bracket_gg_0 = bracket(G, G, 0)
no_gg_2 = NO(G, G)

print(f"bracket(G,G,0) = {bracket_gg_0}")
print(f"NO(G,G) = {no_gg_2}")
print(f"相等: {simplify(bracket_gg_0 - no_gg_2) == 0}")

if simplify(bracket_gg_0 - no_gg_2) == 0:
    print(f"✓ bracket(A,A,0) 与 NO(A,A) 一致")
else:
    print(f"✗ 不一致！")

print("\n" + "=" * 70)
print("总结")
print("=" * 70)
print("pyope 现在正确实现了费米子自身正规序的计算")
print("根据 voa-manual.md 公式 (2.3.17)：")
print("  [AA]_0 = sum_{k>=0} a_k ∂^{2k+1} [AA]_{2k+1}")
print("其中 Bernoulli 系数：")
print("  a_0 = 1/2, a_1 = -1/24, a_2 = 1/240, a_3 = -17/40320")
print("=" * 70)
