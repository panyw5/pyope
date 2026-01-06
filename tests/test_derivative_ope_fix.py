"""
测试右侧导数算符 OPE 的修复

验证 OPE(T, ∂J) 的正确性，其中 J 是 Virasoro primary
"""

import pytest
import sympy as sp
from sympy import Symbol, simplify

from pyope import (
    BasisOperator,
    d,
    One,
    OPE,
    MakeOPE,
    Bosonic,
)


class TestDerivativeOPEFix:
    """测试右侧导数算符 OPE 的修复"""

    def test_virasoro_primary_derivative_ope(self):
        """
        测试 Virasoro primary 的导数 OPE

        如果 J 是 conformal weight = 1 的 Virasoro primary：
        - OPE(T, J) = J/(z-w)² + ∂J/(z-w)
        - OPE(T, ∂J) = 2J/(z-w)³ + 2∂J/(z-w)² + ∂²J/(z-w)
        """
        # 定义算符
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(T, J)

        # 定义 OPE(T, J) - J 是 Virasoro primary
        OPE[T, J] = MakeOPE([
            J,      # pole_2
            d(J)    # pole_1
        ])

        # 计算 OPE(T, ∂J)
        dJ = d(J)
        ope_result = OPE(T, dJ)

        print(f"\nOPE(T, ∂J) 的结果:")
        print(f"  最高极点: {ope_result.max_pole}")

        # 验证最高极点
        assert ope_result.max_pole == 3, f"Expected max_pole=3, got {ope_result.max_pole}"

        # 验证各个极点
        pole_3 = ope_result.pole(3)
        pole_2 = ope_result.pole(2)
        pole_1 = ope_result.pole(1)

        print(f"  pole(3): {pole_3}")
        print(f"  pole(2): {pole_2}")
        print(f"  pole(1): {pole_1}")

        # pole_3 应该是 2J
        assert pole_3 == 2*J, f"Expected pole(3)=2*J, got {pole_3}"

        # pole_2 应该是 2∂J
        assert pole_2 == 2*d(J), f"Expected pole(2)=2*∂J, got {pole_2}"

        # pole_1 应该是 ∂²J
        assert pole_1 == d(J, 2), f"Expected pole(1)=∂²J, got {pole_1}"

        print("✓ OPE(T, ∂J) 计算正确！")

    def test_second_derivative_ope(self):
        """
        测试二阶导数的 OPE

        OPE(T, ∂²J) 应该有更高阶的极点
        """
        # 定义算符
        T = BasisOperator("T", bosonic=True, conformal_weight=2)
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(T, J)

        # 定义 OPE(T, J)
        OPE[T, J] = MakeOPE([J, d(J)])

        # 计算 OPE(T, ∂²J)
        d2J = d(J, 2)
        ope_result = OPE(T, d2J)

        print(f"\nOPE(T, ∂²J) 的结果:")
        print(f"  最高极点: {ope_result.max_pole}")

        # 验证最高极点应该是 4
        assert ope_result.max_pole == 4, f"Expected max_pole=4, got {ope_result.max_pole}"

        # 打印各个极点
        for q in range(ope_result.max_pole, 0, -1):
            pole_q = ope_result.pole(q)
            if pole_q != 0:
                print(f"  pole({q}): {pole_q}")

        print("✓ OPE(T, ∂²J) 计算正确！")

    def test_general_derivative_formula(self):
        """
        测试一般的导数公式

        验证 [A, ∂B]_q = (q-1)[A,B]_{q-1} + ∂[A,B]_q
        """
        # 定义算符
        A = BasisOperator("A", bosonic=True, conformal_weight=2)
        B = BasisOperator("B", bosonic=True, conformal_weight=1)
        Bosonic(A, B)

        # 定义一个简单的 OPE
        OPE[A, B] = MakeOPE([
            3*B,    # pole_2: {AB}_2 = 3B
            2*d(B)  # pole_1: {AB}_1 = 2∂B
        ])

        # 计算 OPE(A, ∂B)
        dB = d(B)
        ope_result = OPE(A, dB)

        # 根据公式：
        # {A∂B}_3 = (3-1)*{AB}_2 = 2*3B = 6B
        # {A∂B}_2 = (2-1)*{AB}_1 + ∂{AB}_2 = 2∂B + 3∂B = 5∂B
        # {A∂B}_1 = ∂{AB}_1 = 2∂²B

        pole_3 = ope_result.pole(3)
        pole_2 = ope_result.pole(2)
        pole_1 = ope_result.pole(1)

        print(f"\nOPE(A, ∂B) 的结果:")
        print(f"  pole(3): {pole_3} (expected: 6B)")
        print(f"  pole(2): {pole_2} (expected: 5∂B)")
        print(f"  pole(1): {pole_1} (expected: 2∂²B)")

        # 验证
        assert pole_3 == 6*B, f"Expected 6*B, got {pole_3}"

        # pole_2 的验证需要展开
        expected_pole_2 = 5*d(B)
        diff = sp.expand(pole_2 - expected_pole_2)
        assert diff == 0, f"Expected 5*∂B, got {pole_2}"

        assert pole_1 == 2*d(B, 2), f"Expected 2*∂²B, got {pole_1}"

        print("✓ 一般导数公式验证通过！")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
