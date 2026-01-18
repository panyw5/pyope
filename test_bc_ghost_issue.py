"""
测试 bc ghost 系统中 NO(b,b) 的化简问题

在 OPEdefs.m 中，NO(b,b) 会自动化简为 0（因为费米子的自身正规序为零）
但在 pyope 中，NO(b,b) 不会自动化简，即使使用 simplify 也不会。

这个测试文件用于重现和分析这个问题。
"""

import sys

sys.path.insert(0, "src")

from pyope import (
    BasisOperator,
    OPE,
    MakeOPE,
    One,
    Zero,
    NO,
    bracket,
    Fermionic,
    simplify,
)

# 设置 bc ghost 系统
b = BasisOperator("b", bosonic=False, conformal_weight=-1)
c = BasisOperator("c", bosonic=False, conformal_weight=2)

# 声明它们是费米子
Fermionic(b, c)

# 定义 OPE: b(z)c(w) ~ 1/(z-w)
OPE[b, c] = MakeOPE([One])

print("=" * 60)
print("测试 bc ghost 系统中的 NO(b,b)")
print("=" * 60)

# 测试 NO(b,b)
print("\n1. 计算 NO(b,b):")
no_bb = NO(b, b)
print(f"   NO(b,b) = {no_bb}")
print(f"   类型: {type(no_bb)}")

# 测试 simplify
print("\n2. 使用 simplify:")
simplified = simplify(no_bb)
print(f"   simplify(NO(b,b)) = {simplified}")
print(f"   类型: {type(simplified)}")

# 检查是否为零
print("\n3. 检查是否为零:")
print(f"   NO(b,b) == 0: {no_bb == 0}")
print(f"   NO(b,b) == Zero: {no_bb == Zero}")
print(f"   simplified == 0: {simplified == 0}")
print(f"   simplified == Zero: {simplified == Zero}")

# 测试 NO(b,b) 的 OPE
print("\n4. 测试 NO(b,b) 与其他算符的 OPE:")
try:
    ope_result = OPE(no_bb, c)
    print(f"   OPE(NO(b,b), c) = {ope_result}")
    print(f"   最高极点: {ope_result.max_pole}")
    if ope_result.max_pole > 0:
        print(f"   pole(1) = {ope_result.pole(1)}")
except Exception as e:
    print(f"   错误: {e}")

# 分析问题
print("\n" + "=" * 60)
print("问题分析:")
print("=" * 60)

print("\n在 OPEdefs.m 中:")
print("- 费米子算符 b 与自身的 OPE 是 b(z)b(w) ~ 0")
print("- 因此 NO(b,b) = {bb}_0 = 0")
print("- OPEdefs.m 会自动识别并化简为 0")

print("\n在 pyope 中:")
print("- NO(b,b) 返回一个 NormalOrderedOperator 对象")
print("- 这个对象不会自动化简为 0")
print("- 即使使用 simplify() 也不会化简")

print("\n可能的原因:")
print("1. pyope 没有定义 b(z)b(w) 的 OPE（因为用户只定义了 b(z)c(w)）")
print("2. NO() 函数没有检查费米子自身的正规序应该为零")
print("3. simplify() 函数没有识别这种特殊情况")

# 验证：检查 b(z)b(w) 的 OPE
print("\n5. 验证 b(z)b(w) 的 OPE:")
try:
    ope_bb = OPE(b, b)
    print(f"   OPE(b,b) = {ope_bb}")
    print(f"   最高极点: {ope_bb.max_pole}")
    if ope_bb.max_pole > 0:
        for i in range(1, ope_bb.max_pole + 1):
            print(f"   pole({i}) = {ope_bb.pole(i)}")

    # 计算 bracket
    bracket_0 = bracket(b, b, 0)
    print(f"   bracket(b,b,0) = {bracket_0}")
    print(f"   这应该等于 NO(b,b)")
except Exception as e:
    print(f"   错误: {e}")

print("\n" + "=" * 60)
print("结论:")
print("=" * 60)
print("对于费米子算符 b，理论上应该有：")
print("- b(z)b(w) ~ 0（费米子的反对易性）")
print("- NO(b,b) = 0")
print("\n但 pyope 需要显式处理这种情况。")
