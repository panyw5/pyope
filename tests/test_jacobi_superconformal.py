"""
Jacobi 恒等式测试 - N=1 超共形代数

本测试使用 pyope 计算 N=1 超共形代数的 Jacobi 恒等式，
并与 Mathematica OPEdefs.m 的参考结果进行对比。

N=1 超共形代数包含：
- T(z): 能动张量，h=2，玻色子（parity=0）
- G(z): 超流，h=3/2，费米子（parity=1）

OPE 公式：
1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)
2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)
3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)

测试内容：
1. check_jacobi_identity(G, G, G) - 最关键，测试费米子符号
2. check_jacobi_identity(T, G, G) - 混合玻色-费米
3. check_jacobi_identity(T, T, G) - 验证 G 的变换性质

关键点：
- G 是费米子，parity=1
- G(z)G(w) 中的符号因子: (-1)^(1*1) = -1
- 这是测试费米子处理的关键用例

参考：
- OPEdefs.m 第 1601-1637 行：OPEJacobi 实现
- OPEdefs.m 第 740-743 行：Fermionic 声明
- OPEdefs.m 第 1604 行：费米子符号因子 sign = (-1)^(OPEParity[A] OPEParity[B])
- ref_jacobi_superconformal.wls：Mathematica 参考测试
"""

import pytest
import sympy as sp
import json
from pathlib import Path
from pyope import (
    BasisOperator,
    OPE,
    d,
    One,
    Bosonic,
    Fermionic,
    check_jacobi_identity,
    verify_jacobi_identity,
)


class TestJacobiIdentitySuperconformal:
    """测试 N=1 超共形代数的 Jacobi 恒等式"""

    def setup_method(self):
        """设置测试环境"""
        # 定义符号
        self.c = sp.Symbol('c')

        # 定义算符
        # T: 能动张量，玻色子，共形权重 h=2
        self.T = BasisOperator("T", bosonic=True, conformal_weight=2)
        Bosonic(self.T)

        # G: 超流，费米子，共形权重 h=3/2
        self.G = BasisOperator("G", bosonic=False, conformal_weight=sp.Rational(3, 2))
        Fermionic(self.G)

        # 定义 OPE
        # 1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)
        OPE[self.T, self.T] = OPE.make([self.c/2 * One, 0, 2*self.T, d(self.T)])

        # 2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)
        OPE[self.T, self.G] = OPE.make([0, sp.Rational(3, 2)*self.G, d(self.G)])

        # 3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)
        OPE[self.G, self.G] = OPE.make([2*self.c/3 * One, 0, 2*self.T])

        print("\n" + "="*70)
        print("N=1 Superconformal Algebra Setup")
        print("="*70)
        print(f"T: Energy-momentum tensor (bosonic, h=2)")
        print(f"   parity = {self.T.parity}")
        print(f"G: Supercurrent (fermionic, h=3/2)")
        print(f"   parity = {self.G.parity}")
        print("\nOPEs defined:")
        print("1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)")
        print("2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)")
        print("3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)")
        print("="*70)

    def test_operator_properties(self):
        """测试算符的基本性质"""
        print("\nTest: Operator Properties")
        print("-"*70)

        # 检查 T 是玻色子
        assert self.T.parity == 0, "T should be bosonic (parity=0)"
        print(f"✓ T is bosonic: parity = {self.T.parity}")

        # 检查 G 是费米子
        assert self.G.parity == 1, "G should be fermionic (parity=1)"
        print(f"✓ G is fermionic: parity = {self.G.parity}")

        # 检查共形权重
        assert self.T.conformal_weight == 2, "T should have h=2"
        print(f"✓ T has conformal weight h = {self.T.conformal_weight}")

        assert self.G.conformal_weight == sp.Rational(3, 2), "G should have h=3/2"
        print(f"✓ G has conformal weight h = {self.G.conformal_weight}")

        print("SUCCESS: All operator properties verified!")

    def test_ope_structure(self):
        """测试 OPE 的结构"""
        print("\nTest: OPE Structure")
        print("-"*70)

        # 测试 T(z)T(w)
        ope_TT = OPE(self.T, self.T)
        print("T(z)T(w):")
        print(f"  Pole 4: {ope_TT.pole(4)}")
        print(f"  Pole 3: {ope_TT.pole(3)}")
        print(f"  Pole 2: {ope_TT.pole(2)}")
        print(f"  Pole 1: {ope_TT.pole(1)}")
        assert ope_TT.pole(4) == self.c/2 * One
        assert ope_TT.pole(3) == 0
        assert ope_TT.pole(2) == 2*self.T
        assert ope_TT.pole(1) == d(self.T)

        # 测试 T(z)G(w)
        ope_TG = OPE(self.T, self.G)
        print("\nT(z)G(w):")
        print(f"  Pole 2: {ope_TG.pole(2)}")
        print(f"  Pole 1: {ope_TG.pole(1)}")
        assert ope_TG.pole(2) == sp.Rational(3, 2)*self.G
        assert ope_TG.pole(1) == d(self.G)

        # 测试 G(z)G(w)
        ope_GG = OPE(self.G, self.G)
        print("\nG(z)G(w):")
        print(f"  Pole 3: {ope_GG.pole(3)}")
        print(f"  Pole 2: {ope_GG.pole(2)}")
        print(f"  Pole 1: {ope_GG.pole(1)}")
        assert ope_GG.pole(3) == 2*self.c/3 * One
        assert ope_GG.pole(2) == 0
        assert ope_GG.pole(1) == 2*self.T

        print("\nSUCCESS: All OPE structures verified!")

    def test_jacobi_GGG(self):
        """测试 Jacobi 恒等式 check_jacobi_identity(G, G, G)

        这是最关键的测试，验证费米子符号处理。
        符号因子: (-1)^(parity[G] * parity[G]) = (-1)^(1*1) = -1
        """
        print("\nTest: Jacobi Identity check_jacobi_identity(G, G, G)")
        print("-"*70)
        print("This is the most critical test for fermionic sign factors")
        print("Sign factor: (-1)^(1*1) = -1")
        print("\nComputing...")

        # 计算 Jacobi 恒等式
        jacobi_result = check_jacobi_identity(self.G, self.G, self.G, simplify_func=sp.expand)

        print(f"Result dimensions: {len(jacobi_result)} x {len(jacobi_result[0]) if jacobi_result else 0}")
        print("Jacobi identity result:")

        # 打印结果矩阵
        for i, row in enumerate(jacobi_result):
            row_str = "  [" + ", ".join(str(val) for val in row) + "]"
            print(row_str)

        # 检查所有元素是否为 0
        all_zero = True
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                if value != 0:
                    print(f"  Non-zero entry at [{i+1},{j+1}]: {value}")
                    all_zero = False

        if all_zero:
            print("  All entries are zero!")
            print("SUCCESS: Jacobi identity holds for G-G-G!")
        else:
            print("  Some entries are non-zero!")
            print("FAILURE: Jacobi identity violated for G-G-G!")

        # 断言所有元素为 0
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                assert value == 0, f"Jacobi identity violated at position [{i+1},{j+1}]: {value}"

    def test_jacobi_TGG(self):
        """测试 Jacobi 恒等式 check_jacobi_identity(T, G, G)

        混合玻色-费米测试。
        符号因子: (-1)^(parity[T] * parity[G]) = (-1)^(0*1) = 1
        """
        print("\nTest: Jacobi Identity check_jacobi_identity(T, G, G)")
        print("-"*70)
        print("Mixed bosonic-fermionic test")
        print("Sign factor: (-1)^(0*1) = 1")
        print("\nComputing...")

        # 计算 Jacobi 恒等式
        jacobi_result = check_jacobi_identity(self.T, self.G, self.G, simplify_func=sp.expand)

        print(f"Result dimensions: {len(jacobi_result)} x {len(jacobi_result[0]) if jacobi_result else 0}")
        print("Jacobi identity result:")

        # 打印结果矩阵
        for i, row in enumerate(jacobi_result):
            row_str = "  [" + ", ".join(str(val) for val in row) + "]"
            print(row_str)

        # 检查所有元素是否为 0
        all_zero = True
        non_zero_entries = []
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                if value != 0:
                    print(f"  Non-zero entry at [{i+1},{j+1}]: {value}")
                    non_zero_entries.append((i+1, j+1, value))
                    all_zero = False

        if all_zero:
            print("  All entries are zero!")
            print("SUCCESS: Jacobi identity holds for T-G-G!")
        else:
            print("  Some entries are non-zero!")
            print("FAILURE: Jacobi identity violated for T-G-G!")

        # 断言所有元素为 0
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                assert value == 0, f"Jacobi identity violated at position [{i+1},{j+1}]: {value}"

    def test_jacobi_TTG(self):
        """测试 Jacobi 恒等式 check_jacobi_identity(T, T, G)

        验证 G 的变换性质。
        符号因子: (-1)^(parity[T] * parity[T]) = (-1)^(0*0) = 1
        """
        print("\nTest: Jacobi Identity check_jacobi_identity(T, T, G)")
        print("-"*70)
        print("Verifies transformation properties of G")
        print("Sign factor: (-1)^(0*0) = 1")
        print("\nComputing...")

        # 计算 Jacobi 恒等式
        jacobi_result = check_jacobi_identity(self.T, self.T, self.G, simplify_func=sp.expand)

        print(f"Result dimensions: {len(jacobi_result)} x {len(jacobi_result[0]) if jacobi_result else 0}")
        print("Jacobi identity result:")

        # 打印结果矩阵
        for i, row in enumerate(jacobi_result):
            row_str = "  [" + ", ".join(str(val) for val in row) + "]"
            print(row_str)

        # 检查所有元素是否为 0
        all_zero = True
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                if value != 0:
                    print(f"  Non-zero entry at [{i+1},{j+1}]: {value}")
                    all_zero = False

        if all_zero:
            print("  All entries are zero!")
            print("SUCCESS: Jacobi identity holds for T-T-G!")
        else:
            print("  Some entries are non-zero!")
            print("FAILURE: Jacobi identity violated for T-T-G!")

        # 断言所有元素为 0
        for i, row in enumerate(jacobi_result):
            for j, value in enumerate(row):
                assert value == 0, f"Jacobi identity violated at position [{i+1},{j+1}]: {value}"

    def test_comparison_with_mathematica(self):
        """与 Mathematica 参考结果对比"""
        print("\nTest: Comparison with Mathematica Reference")
        print("-"*70)

        # 尝试加载 Mathematica 参考结果
        ref_file = Path(__file__).parent / "ref_jacobi_superconformal_output.json"

        if not ref_file.exists():
            print(f"WARNING: Reference file not found: {ref_file}")
            print("Skipping comparison test.")
            pytest.skip("Mathematica reference file not available")
            return

        with open(ref_file, 'r') as f:
            ref_data = json.load(f)

        print("Loaded Mathematica reference results")
        print(f"  GGG dimensions: {ref_data['dimensions_GGG']}")
        print(f"  TGG dimensions: {ref_data['dimensions_TGG']}")
        print(f"  TTG dimensions: {ref_data['dimensions_TTG']}")
        print(f"  GGG all zero: {ref_data['all_zero_GGG']}")
        print(f"  TGG all zero: {ref_data['all_zero_TGG']}")
        print(f"  TTG all zero: {ref_data['all_zero_TTG']}")

        # 计算 pyope 结果
        pyope_GGG = check_jacobi_identity(self.G, self.G, self.G, simplify_func=sp.expand)
        pyope_TGG = check_jacobi_identity(self.T, self.G, self.G, simplify_func=sp.expand)
        pyope_TTG = check_jacobi_identity(self.T, self.T, self.G, simplify_func=sp.expand)

        print("\nPyope results:")
        print(f"  GGG dimensions: {len(pyope_GGG)} x {len(pyope_GGG[0])}")
        print(f"  TGG dimensions: {len(pyope_TGG)} x {len(pyope_TGG[0])}")
        print(f"  TTG dimensions: {len(pyope_TTG)} x {len(pyope_TTG[0])}")

        # 对比 all_zero 状态（这是最重要的）
        print("\nComparing all_zero status (most important):")
        pyope_all_zero_GGG = all(val == 0 for row in pyope_GGG for val in row)
        pyope_all_zero_TGG = all(val == 0 for row in pyope_TGG for val in row)
        pyope_all_zero_TTG = all(val == 0 for row in pyope_TTG for val in row)

        assert pyope_all_zero_GGG == ref_data['all_zero_GGG'], "GGG all_zero status mismatch"
        print(f"  ✓ GGG all_zero matches: {pyope_all_zero_GGG}")

        assert pyope_all_zero_TGG == ref_data['all_zero_TGG'], "TGG all_zero status mismatch"
        print(f"  ✓ TGG all_zero matches: {pyope_all_zero_TGG}")

        assert pyope_all_zero_TTG == ref_data['all_zero_TTG'], "TTG all_zero status mismatch"
        print(f"  ✓ TTG all_zero matches: {pyope_all_zero_TTG}")

        print("\nNote: Dimensions may differ between Mathematica and pyope")
        print("      due to different range calculations, but the all_zero")
        print("      status (which is what matters for Jacobi identity)")
        print("      matches perfectly!")

        print("\nSUCCESS: pyope results match Mathematica reference!")
        print("  All Jacobi identities hold with the corrected OPE formula!")


class TestJacobiIdentityVerifyFunction:
    """测试 verify_jacobi_identity 便捷函数"""

    def setup_method(self):
        """设置测试环境"""
        self.c = sp.Symbol('c')
        self.T = BasisOperator("T", bosonic=True, conformal_weight=2)
        self.G = BasisOperator("G", bosonic=False, conformal_weight=sp.Rational(3, 2))
        Bosonic(self.T)
        Fermionic(self.G)

        OPE[self.T, self.T] = OPE.make([self.c/2 * One, 0, 2*self.T, d(self.T)])
        OPE[self.T, self.G] = OPE.make([0, sp.Rational(3, 2)*self.G, d(self.G)])
        OPE[self.G, self.G] = OPE.make([2*self.c/3 * One, 0, 2*self.T])

    def test_verify_jacobi_GGG(self):
        """测试 verify_jacobi_identity(G, G, G)"""
        print("\nTest: verify_jacobi_identity(G, G, G)")
        print("-"*70)

        result = verify_jacobi_identity(self.G, self.G, self.G, simplify_func=sp.expand)
        print(f"Result: {result}")

        assert result is True, "Jacobi identity should hold for G-G-G"
        print("SUCCESS: verify_jacobi_identity confirms Jacobi identity holds!")

    def test_verify_jacobi_TTG(self):
        """测试 verify_jacobi_identity(T, T, G)"""
        print("\nTest: verify_jacobi_identity(T, T, G)")
        print("-"*70)

        result = verify_jacobi_identity(self.T, self.T, self.G, simplify_func=sp.expand)
        print(f"Result: {result}")

        assert result is True, "Jacobi identity should hold for T-T-G"
        print("SUCCESS: verify_jacobi_identity confirms Jacobi identity holds!")


if __name__ == "__main__":
    # 运行测试
    print("="*70)
    print("Jacobi Identity Test for N=1 Superconformal Algebra")
    print("="*70)

    test = TestJacobiIdentitySuperconformal()
    test.setup_method()

    print("\n" + "="*70)
    print("Test 1: Operator Properties")
    print("="*70)
    test.test_operator_properties()

    print("\n" + "="*70)
    print("Test 2: OPE Structure")
    print("="*70)
    test.test_ope_structure()

    print("\n" + "="*70)
    print("Test 3: Jacobi Identity G-G-G")
    print("="*70)
    test.test_jacobi_GGG()

    print("\n" + "="*70)
    print("Test 4: Jacobi Identity T-G-G")
    print("="*70)
    test.test_jacobi_TGG()

    print("\n" + "="*70)
    print("Test 5: Jacobi Identity T-T-G")
    print("="*70)
    test.test_jacobi_TTG()

    print("\n" + "="*70)
    print("Test 6: Comparison with Mathematica")
    print("="*70)
    test.test_comparison_with_mathematica()

    print("\n" + "="*70)
    print("Test 7: verify_jacobi_identity Functions")
    print("="*70)
    test2 = TestJacobiIdentityVerifyFunction()
    test2.setup_method()
    test2.test_verify_jacobi_GGG()
    test2.test_verify_jacobi_TTG()

    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print("✓ G-G-G: PASS (fermionic sign factor works correctly)")
    print("✓ T-G-G: PASS (corrected OPE formula with 2c/3)")
    print("✓ T-T-G: PASS")
    print("\nNote: The correct G(z)G(w) OPE formula is:")
    print("      G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)")
    print("      NOT (c/3)/(z-w)^3 + 2T(w)/(z-w)")
    print("      This is confirmed by Thielemans' paper (eq. 2.4.13)")
    print("="*70)
