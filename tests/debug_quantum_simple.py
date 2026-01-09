"""
简单的量子数计算测试

直接测试基本算符的量子数计算
"""

from fractions import Fraction
from pyope.operators import BasisOperator, d
from pyope.api import NO
from pyope.null_states import QuantumNumberCalculator

# 定义自由场
b = BasisOperator('b', 0, 2)
c = BasisOperator('c', 0, -1)
beta = BasisOperator('β', 0, Fraction(3, 2))
gamma = BasisOperator('γ', 0, Fraction(-1, 2))

# 定义简单的生成元
w = beta
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
gw = b

# 量子数映射
quantum_number_map = {
    b: (Fraction(1), Fraction(1, 2)),
    c: (Fraction(-1), Fraction(-1, 2)),
    beta: (Fraction(3, 2), Fraction(0)),
    gamma: (Fraction(-3, 2), Fraction(0)),
    w: (Fraction(3, 2), Fraction(0)),
    j0: (Fraction(0), Fraction(0)),
    gw: (Fraction(1), Fraction(1, 2))
}

# 创建量子数计算器
quantum_calc = QuantumNumberCalculator(quantum_number_map)

print("=" * 70)
print("测试基本算符的量子数计算")
print("=" * 70)

# 测试 1: 简单的 BasisOperator
print("\n测试 1: 简单的 BasisOperator")
print(f"b 的量子数: {quantum_calc.get_quantum_numbers(b)}")
print(f"beta 的量子数: {quantum_calc.get_quantum_numbers(beta)}")

# 测试 2: DerivativeOperator
print("\n测试 2: DerivativeOperator")
print(f"∂b 的量子数: {quantum_calc.get_quantum_numbers(d(b))}")
print(f"∂^2b 的量子数: {quantum_calc.get_quantum_numbers(d(b, 2))}")

# 测试 3: NormalOrderedOperator
print("\n测试 3: NormalOrderedOperator")
print(f"NO(b,c) 的量子数: {quantum_calc.get_quantum_numbers(NO(b, c))}")
print(f"NO(beta,gamma) 的量子数: {quantum_calc.get_quantum_numbers(NO(beta, gamma))}")

# 测试 4: OperatorSum (j0)
print("\n测试 4: OperatorSum (j0)")
print(f"j0 的类型: {type(j0)}")
print(f"j0 的字符串表示: {str(j0)[:50]}...")
print(f"j0 在映射表中吗? {j0 in quantum_number_map}")
print(f"j0 的量子数: {quantum_calc.get_quantum_numbers(j0)}")

# 测试 5: 包含 j0 的复合算符
print("\n测试 5: 包含 j0 的复合算符")
print(f"∂j0 的量子数: {quantum_calc.get_quantum_numbers(d(j0))}")
print(f"NO(j0, j0) 的量子数: {quantum_calc.get_quantum_numbers(NO(j0, j0))}")

# 测试 6: 检查对象相等性
print("\n测试 6: 检查对象相等性")
j0_copy = 2 * NO(b, c) + 3 * NO(beta, gamma)
print(f"j0 == j0_copy: {j0 == j0_copy}")
print(f"j0 is j0_copy: {j0 is j0_copy}")
print(f"id(j0): {id(j0)}")
print(f"id(j0_copy): {id(j0_copy)}")

# 测试 7: 检查 OperatorSum 的 terms
print("\n测试 7: 检查 OperatorSum 的 terms")
if hasattr(j0, 'terms'):
    print(f"j0.terms 的类型: {type(j0.terms)}")
    print(f"j0.terms 的长度: {len(j0.terms)}")
    print(f"j0.terms 的键:")
    for i, key in enumerate(list(j0.terms.keys())[:3]):
        print(f"  {i+1}. {key} (类型: {type(key)})")
