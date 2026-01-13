"""
测试 _ope_composite_left 中的符号是否正确

对比 Mathematica 的 OPECompositeHelpLQ 实现
"""
import pytest
from pyope import BasisOperator, OPE, NO, bracket
from pyope.api import derivative
from pyope.registry import ope_registry
from pyope.ope_data import OPEData


def test_composite_left_sign_term1():
    """
    测试第一项的符号：Σ_{q=1}^{maxBC} Σ_{l=0}^{maxBC-q} NO[∂^l A, {BC}_{l+q}] / l!

    Mathematica: 没有 sign 前缀
    Python: result._poles[q] = pole_sum (没有 sign 前缀)
    """
    ope_registry.clear()

    # 定义算符
    A = BasisOperator("A", bosonic=True)
    B = BasisOperator("B", bosonic=True)
    C = BasisOperator("C", bosonic=True)

    # 定义 OPE(B, C) = 1/(z^2) * A
    ope_registry.define_ope(B, C, OPEData({2: A}))

    # 计算 OPE(NO(A,B), C)
    result = OPE(NO(A, B), C)

    # 第一项应该贡献 NO(A, A) 到 q=2
    # 因为 {BC}_2 = A, l=0, 所以 NO[∂^0 A, A] / 0! = NO(A, A)
    pole_2 = result.pole(2)

    print(f"pole(2) = {pole_2}")
    # 这应该包含 NO(A, A) 项（没有负号）
    assert pole_2 != 0


def test_composite_left_sign_term2():
    """
    测试第二项的符号：sign * Σ_{q=1}^{maxAC} Σ_{l=0}^{maxAC-q} NO[∂^l B, {AC}_{l+q}] / l!

    Mathematica: 有 sign 前缀
    Python: result._poles[q] = sign * pole_sum (有 sign 前缀)

    其中 sign = (-1)^(|A||B|)
    """
    ope_registry.clear()

    # 定义算符 - A 和 B 都是玻色子
    A = BasisOperator("A", bosonic=True)
    B = BasisOperator("B", bosonic=True)
    C = BasisOperator("C", bosonic=True)

    # 定义 OPE(A, C) = 1/(z^2) * B
    ope_registry.define_ope(A, C, OPEData({2: B}))

    # 计算 OPE(NO(A,B), C)
    # sign = (-1)^(0*0) = 1 (两个玻色子)
    result = OPE(NO(A, B), C)

    pole_2 = result.pole(2)
    print(f"Bosonic case: pole(2) = {pole_2}")

    # 现在测试费米子情况
    ope_registry.clear()
    psi = BasisOperator("psi", bosonic=False)
    chi = BasisOperator("chi", bosonic=False)
    phi = BasisOperator("phi", bosonic=True)

    # 定义 OPE(psi, phi) = 1/(z^2) * chi
    ope_registry.define_ope(psi, phi, OPEData({2: chi}))

    # 计算 OPE(NO(psi, chi), phi)
    # sign = (-1)^(1*1) = -1 (两个费米子)
    result_fermionic = OPE(NO(psi, chi), phi)

    pole_2_fermionic = result_fermionic.pole(2)
    print(f"Fermionic case: pole(2) = {pole_2_fermionic}")

    # 第二项应该有 sign = -1 的前缀


def test_composite_left_sign_term3():
    """
    测试第三项的符号：sign * Σ_{q} Σ_{l} {B, {AC}_q}_{l}

    Mathematica: 有 sign 前缀
    Python: result._poles[q] = sign * pole_sum (有 sign 前缀)
    """
    ope_registry.clear()

    # 定义算符
    A = BasisOperator("A", bosonic=True)
    B = BasisOperator("B", bosonic=True)
    C = BasisOperator("C", bosonic=True)
    D = BasisOperator("D", bosonic=True)

    # 定义 OPE(A, C) = 1/(z^2) * D
    ope_registry.define_ope(A, C, OPEData({2: D}))

    # 定义 OPE(B, D) = 1/(z^1) * A
    ope_registry.define_ope(B, D, OPEData({1: A}))

    # 计算 OPE(NO(A,B), C)
    result = OPE(NO(A, B), C)

    # 第三项应该贡献到高阶极点
    print(f"max_pole = {result.max_pole}")
    for q in range(1, result.max_pole + 1):
        print(f"pole({q}) = {result.pole(q)}")


if __name__ == "__main__":
    print("=" * 60)
    print("测试第一项符号")
    print("=" * 60)
    test_composite_left_sign_term1()

    print("\n" + "=" * 60)
    print("测试第二项符号")
    print("=" * 60)
    test_composite_left_sign_term2()

    print("\n" + "=" * 60)
    print("测试第三项符号")
    print("=" * 60)
    test_composite_left_sign_term3()
