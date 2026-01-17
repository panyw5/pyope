#!/usr/bin/env python
"""
演示 simplify() 的导数自动展开功能（莱布尼茨法则）

这个脚本展示了如何使用 simplify() 自动展开正规序算符的导数。
"""

from pyope import BasisOperator, d, NO, simplify, Bosonic

print("=" * 70)
print("演示：simplify() 自动展开导数（莱布尼茨法则）")
print("=" * 70)

# 定义基本算符
T = BasisOperator("T", bosonic=True)
J = BasisOperator("J", bosonic=True)
Bosonic(T, J)

print("\n1. 一阶导数展开")
print("-" * 70)
expr1 = d(NO(T, J))
print(f"原始表达式: {expr1}")
result1 = simplify(expr1)
print(f"展开结果:   {result1}")
print(f"公式: d(NO(A,B)) = NO(d(A), B) + NO(A, d(B))")

print("\n2. 二阶导数展开")
print("-" * 70)
expr2 = d(NO(T, J), 2)
print(f"原始表达式: {expr2}")
result2 = simplify(expr2)
print(f"展开结果:   {result2}")
print(f"公式: d²(NO(A,B)) = NO(d²(A), B) + 2*NO(d(A), d(B)) + NO(A, d²(B))")

print("\n3. 三阶导数展开")
print("-" * 70)
expr3 = d(NO(T, J), 3)
print(f"原始表达式: {expr3}")
result3 = simplify(expr3)
print(f"展开结果:   {result3}")
print(f"公式: d³(NO(A,B)) = NO(d³(A), B) + 3*NO(d²(A), d(B)) + 3*NO(d(A), d²(B)) + NO(A, d³(B))")

print("\n4. 禁用导数展开")
print("-" * 70)
expr4 = d(NO(T, J))
print(f"原始表达式: {expr4}")
result4 = simplify(expr4, expand_derivatives=False)
print(f"禁用展开:   {result4}")
print(f"说明: 使用 expand_derivatives=False 可以保持导数未展开")

print("\n5. 嵌套正规序的导数")
print("-" * 70)
W = BasisOperator("W", bosonic=True)
Bosonic(W)
expr5 = d(NO(NO(T, J), W))
print(f"原始表达式: {expr5}")
result5 = simplify(expr5)
print(f"展开结果:   {result5}")
print(f"说明: 递归展开嵌套的 NO")

print("\n6. 带标量系数的导数")
print("-" * 70)
expr6 = 3 * d(NO(T, J))
print(f"原始表达式: {expr6}")
result6 = simplify(expr6)
print(f"展开结果:   {result6}")
print(f"说明: 标量系数会分配到每一项")

print("\n" + "=" * 70)
print("演示完成！")
print("=" * 70)
