"""
测试 bracket(A, B, n=0) 的正确行为

验证：
1. bracket(A, B, 0) 等同于 NO(A, B)
2. bracket(A, NO(B,C), 0) 等同于 NO(A, NO(B,C))
3. 确保不会错误地使用 _ope_composite_right 计算 q=0
"""

import sympy as sp
from pyope import BasisOperator, bracket, NO, OPE, MakeOPE, d
from pyope.constants import One

def test_bracket_n0_basic():
    """测试基本情况：bracket(A, B, 0) = NO(A, B)"""
    # 创建基础算符
    T = BasisOperator("T", bosonic=True)
    J = BasisOperator("J", bosonic=True)

    # 定义一个简单的 OPE（不重要，因为 n=0 不应该使用它）
    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c/2*One, 0, 2*T, d(T)])

    # bracket(T, J, 0) 应该等同于 NO(T, J)
    result1 = bracket(T, J, 0)
    result2 = NO(T, J)

    # 两者应该相等（都是 NormalOrderedOperator 或简化后的相同表达式）
    assert str(result1) == str(result2), f"Expected {result2}, got {result1}"
    print(f"✓ bracket(T, J, 0) = {result1} = NO(T, J)")


def test_bracket_n0_composite():
    """测试复合算符：bracket(A, NO(B,C), 0) = NO(A, NO(B,C))"""
    # 创建基础算符
    T = BasisOperator("T", bosonic=True)
    J = BasisOperator("J", bosonic=True)

    # 定义 OPE
    c = sp.Symbol("c")
    OPE[T, T] = MakeOPE([c/2*One, 0, 2*T, d(T)])
    OPE[T, J] = MakeOPE([J, d(J)])  # 假设的 OPE

    # 创建复合算符 NO(T, J)
    composite = NO(T, J)

    # bracket(J, NO(T, J), 0) 应该等同于 NO(J, NO(T, J))
    result1 = bracket(J, composite, 0)
    result2 = NO(J, composite)

    # 两者应该相等
    assert str(result1) == str(result2), f"Expected {result2}, got {result1}"
    print(f"✓ bracket(J, NO(T,J), 0) = {result1} = NO(J, NO(T,J))")


def test_bracket_n0_not_from_ope():
    """
    验证 bracket(A, B, 0) 不会从 OPE 中提取第 0 阶极点

    即使 OPE 中没有定义第 0 阶极点，bracket(A, B, 0) 也应该返回 NO(A, B)
    """
    # 创建基础算符
    A = BasisOperator("A", bosonic=True)
    B = BasisOperator("B", bosonic=True)

    # 定义一个只有奇异部分的 OPE（没有 q=0 项）
    OPE[A, B] = MakeOPE([A, B])  # 只有 q=2 和 q=1

    # bracket(A, B, 0) 应该返回 NO(A, B)，而不是 0
    result = bracket(A, B, 0)
    expected = NO(A, B)

    assert str(result) == str(expected), f"Expected {expected}, got {result}"
    assert result != 0, "bracket(A, B, 0) should not return 0 even if OPE has no q=0 pole"
    print(f"✓ bracket(A, B, 0) = {result} (not extracted from OPE)")


if __name__ == "__main__":
    print("Testing bracket(A, B, n=0) behavior...\n")
    test_bracket_n0_basic()
    test_bracket_n0_composite()
    test_bracket_n0_not_from_ope()
    print("\n✅ All tests passed!")
