"""
Jacobi 恒等式测试 - Virasoro 代数

本测试使用 pyope 计算 Virasoro 算符 T 的 Jacobi 恒等式，
并与 Mathematica OPEdefs.m 的参考结果进行对比。

测试内容：
1. 定义 Virasoro 算符 T，共形权重为 2
2. 使用 check_jacobi_identity(T, T, T) 检查 Jacobi 恒等式
3. 验证结果全为 0（Jacobi 恒等式成立）
4. 与 Mathematica 参考结果对比

参考：
- OPEdefs.m 第 1598-1637 行：OPEJacobi 实现
- ref_jacobi_virasoro.wls：Mathematica 参考测试
"""

import pytest
import sympy as sp
from pyope import (
    BasisOperator,
    OPE,
    d,
    One,
    Bosonic,
    check_jacobi_identity,
    verify_jacobi_identity,
)


class TestJacobiIdentityVirasoro:
    """测试 Virasoro 代数的 Jacobi 恒等式"""

    def setup_method(self):
        """设置测试环境"""
        # 定义符号
        self.c = sp.Symbol('c')

        # 定义 Virasoro 算符 T
        self.T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(self.T)

        # 定义 T 的 OPE: T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + T'(w)/(z-w)
        OPE[self.T, self.T] = OPE.make([self.c/2 * One, 0, 2*self.T, d(self.T)])

    def test_virasoro_ope_structure(self):
        """测试 Virasoro OPE 的结构"""
        ope_TT = OPE(self.T, self.T)

        # 验证各阶极点
        assert ope_TT.pole(4) == self.c/2 * One, "Pole 4 should be c/2*One"
        assert ope_TT.pole(3) == 0, "Pole 3 should be 0"
        assert ope_TT.pole(2) == 2*self.T, "Pole 2 should be 2*T"
        assert ope_TT.pole(1) == d(self.T), "Pole 1 should be d(T)"

        print("Virasoro OPE structure verified:")
        print(f"  Pole 4: {ope_TT.pole(4)}")
        print(f"  Pole 3: {ope_TT.pole(3)}")
        print(f"  Pole 2: {ope_TT.pole(2)}")
        print(f"  Pole 1: {ope_TT.pole(1)}")

    def test_jacobi_identity_TTT(self):
        """测试 Jacobi 恒等式 check_jacobi_identity(T, T, T)"""
        print("\nComputing Jacobi identity check_jacobi_identity(T, T, T)...")

        # 计算 Jacobi 恒等式
        jacobi_result = check_jacobi_identity(self.T, self.T, self.T, simplify_func=sp.expand)

        print(f"Result dimensions: {len(jacobi_result)} x {len(jacobi_result[0]) if jacobi_result else 0}")
        print("Jacobi identity result:")

        # 检查所有元素是否为 0
        all_zero = True
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                if value != 0:
                    print(f"  Non-zero entry at [{i+1},{j+1}]: {value}")
                    all_zero = False

        if all_zero:
            print("  All entries are zero!")
        else:
            print("  Some entries are non-zero!")

        # 断言所有元素为 0
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                assert value == 0, f"Jacobi identity violated at position [{i+1},{j+1}]: {value}"

        print("SUCCESS: Jacobi identity holds for Virasoro algebra!")

    def test_verify_jacobi_identity_TTT(self):
        """测试 verify_jacobi_identity(T, T, T) 便捷函数"""
        print("\nTesting verify_jacobi_identity(T, T, T)...")

        # 使用便捷函数验证
        result = verify_jacobi_identity(self.T, self.T, self.T, simplify_func=sp.expand)

        print(f"verify_jacobi_identity result: {result}")

        assert result is True, "Jacobi identity should hold for Virasoro algebra"

        print("SUCCESS: verify_jacobi_identity confirms Jacobi identity holds!")

    def test_jacobi_identity_dimensions(self):
        """测试 Jacobi 恒等式结果的维度"""
        print("\nTesting Jacobi identity result dimensions...")

        jacobi_result = check_jacobi_identity(self.T, self.T, self.T)

        # 根据 Mathematica 参考结果，应该是 5x5 矩阵
        # 这是因为 T(z)T(w) 的最高极点是 4，所以需要检查的范围是 1 到 5
        expected_rows = 5
        expected_cols = 5

        actual_rows = len(jacobi_result)
        actual_cols = len(jacobi_result[0]) if jacobi_result else 0

        print(f"Expected dimensions: {expected_rows} x {expected_cols}")
        print(f"Actual dimensions: {actual_rows} x {actual_cols}")

        assert actual_rows == expected_rows, f"Expected {expected_rows} rows, got {actual_rows}"
        assert actual_cols == expected_cols, f"Expected {expected_cols} columns, got {actual_cols}"

        print("SUCCESS: Jacobi identity result has correct dimensions!")

    def test_comparison_with_mathematica(self):
        """与 Mathematica 参考结果对比"""
        print("\nComparing with Mathematica reference results...")

        # Mathematica 参考结果：5x5 全零矩阵
        mathematica_result = [[0]*5 for _ in range(5)]

        # pyope 计算结果
        pyope_result = check_jacobi_identity(self.T, self.T, self.T, simplify_func=sp.expand)

        # 对比维度
        assert len(pyope_result) == len(mathematica_result), "Dimension mismatch"
        assert len(pyope_result[0]) == len(mathematica_result[0]), "Dimension mismatch"

        # 对比每个元素
        for i in range(len(pyope_result)):
            for j in range(len(pyope_result[i])):
                pyope_value = pyope_result[i][j]
                mathematica_value = mathematica_result[i][j]

                assert pyope_value == mathematica_value, \
                    f"Mismatch at [{i+1},{j+1}]: pyope={pyope_value}, mathematica={mathematica_value}"

        print("SUCCESS: pyope results match Mathematica reference!")


class TestJacobiIdentityProperties:
    """测试 Jacobi 恒等式的一般性质"""

    def test_jacobi_identity_for_simple_operators(self):
        """测试简单算符的 Jacobi 恒等式"""
        # 定义一个简单的流算符 J，共形权重为 1
        J = BasisOperator("J", bosonic=True, conformal_weight=1)
        Bosonic(J)

        # 定义 J 的 OPE: J(z)J(w) = k/(z-w)^2
        k = sp.Symbol('k')
        OPE[J, J] = OPE.make([k * One, 0])

        # 验证 Jacobi 恒等式
        result = verify_jacobi_identity(J, J, J, simplify_func=sp.expand)

        assert result is True, "Jacobi identity should hold for current algebra"

        print("SUCCESS: Jacobi identity holds for simple current algebra!")


if __name__ == "__main__":
    # 运行测试
    print("=" * 70)
    print("Jacobi Identity Test for Virasoro Algebra")
    print("=" * 70)

    test = TestJacobiIdentityVirasoro()
    test.setup_method()

    print("\n" + "=" * 70)
    print("Test 1: Virasoro OPE Structure")
    print("=" * 70)
    test.test_virasoro_ope_structure()

    print("\n" + "=" * 70)
    print("Test 2: Jacobi Identity check_jacobi_identity(T, T, T)")
    print("=" * 70)
    test.test_jacobi_identity_TTT()

    print("\n" + "=" * 70)
    print("Test 3: verify_jacobi_identity(T, T, T)")
    print("=" * 70)
    test.test_verify_jacobi_identity_TTT()

    print("\n" + "=" * 70)
    print("Test 4: Jacobi Identity Dimensions")
    print("=" * 70)
    test.test_jacobi_identity_dimensions()

    print("\n" + "=" * 70)
    print("Test 5: Comparison with Mathematica")
    print("=" * 70)
    test.test_comparison_with_mathematica()

    print("\n" + "=" * 70)
    print("Test 6: Simple Current Algebra")
    print("=" * 70)
    test2 = TestJacobiIdentityProperties()
    test2.test_jacobi_identity_for_simple_operators()

    print("\n" + "=" * 70)
    print("All tests passed!")
    print("=" * 70)
