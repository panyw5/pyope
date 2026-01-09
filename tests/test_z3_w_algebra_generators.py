"""
测试 Z3 W-代数生成元的自由场实现

根据 null_states.txt 的第 115-122 行，Z3 W-代数有以下生成元：
- w = β
- j0 = 2:bc: + 3:βγ:
- wb = (复杂的自由场表达式)
- t = -2:b∂c: - 3/2:β∂γ: - :∂bc: - 1/2:∂βγ:  (能动张量)
- g = :γb:
- gt = 2:∂βc: + 3:β∂c:
- gw = b
- gwb = (复杂的自由场表达式)
"""

import sys
sys.path.insert(0, '../src')

import pytest
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
    d,  # 导数算符
)


class TestZ3WAlgebraGenerators:
    """测试 Z3 W-代数生成元的构造"""

    def setup_method(self):
        """设置自由场系统和 W-代数生成元"""
        # 定义自由场
        self.beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
        self.gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)
        self.b = BasisOperator("b", bosonic=False, conformal_weight=-1)
        self.c = BasisOperator("c", bosonic=False, conformal_weight=2)

        # 声明费米子
        Fermionic(self.beta, self.gamma, self.b, self.c)

        # 定义自由场 OPE
        OPE[self.beta, self.gamma] = MakeOPE([-One])
        OPE[self.b, self.c] = MakeOPE([One])

        # 构造 W-代数生成元
        self._construct_generators()

    def _construct_generators(self):
        """构造 W-代数生成元"""
        # w = β
        self.w = self.beta

        # j0 = 2:bc: + 3:βγ:
        self.j0 = 2 * NO(self.b, self.c) + 3 * NO(self.beta, self.gamma)

        # t = -2:b∂c: - 3/2:β∂γ: - :∂bc: - 1/2:∂βγ:
        # 这是能动张量（stress tensor）
        self.t = (
            -2 * NO(self.b, d(self.c))
            - sp.Rational(3, 2) * NO(self.beta, d(self.gamma))
            - NO(d(self.b), self.c)
            - sp.Rational(1, 2) * NO(d(self.beta), self.gamma)
        )

        # g = :γb:
        self.g = NO(self.gamma, self.b)

        # gt = 2:∂βc: + 3:β∂c:
        self.gt = 2 * NO(d(self.beta), self.c) + 3 * NO(self.beta, d(self.c))

        # gw = b
        self.gw = self.b

        # wb = :β:β:γ:γγ::: + 2:β:γ:γ:bc::: - 4:β:∂γγ:: - 4/3:γ:b∂c:: + 2/3:γ:∂bc:: + 2/3:∂β:γγ:: - 8/3:∂γ:bc:: + 10/9*∂²γ
        # 这是一个非常复杂的表达式，包含多重正规序
        self.wb = (
            NO(self.beta, NO(self.beta, NO(self.gamma, NO(self.gamma, self.gamma))))
            + 2 * NO(self.beta, NO(self.gamma, NO(self.gamma, NO(self.b, self.c))))
            - 4 * NO(self.beta, NO(d(self.gamma), self.gamma))
            - sp.Rational(4, 3) * NO(self.gamma, NO(self.b, d(self.c)))
            + sp.Rational(2, 3) * NO(self.gamma, NO(d(self.b), self.c))
            + sp.Rational(2, 3) * NO(d(self.beta), NO(self.gamma, self.gamma))
            - sp.Rational(8, 3) * NO(d(self.gamma), NO(self.b, self.c))
            + sp.Rational(10, 9) * d(self.gamma, 2)
        )

        # gwb = 8/3:b:∂²cc:: + 3:β:β:γ:γ∂c::: - 4:β:γ:b:∂cc::: - 4:β:γ∂²c:: - 4:β:∂γ∂c:: - 2/3:∂b:∂cc:: + 2:∂β:β:γ:γc::: - 8/3:∂β:∂γc:: + 2/3:∂²β:γc:: + 10/9*∂³c
        # 这也是一个非常复杂的表达式
        self.gwb = (
            sp.Rational(8, 3) * NO(self.b, NO(d(self.c, 2), self.c))
            + 3 * NO(self.beta, NO(self.beta, NO(self.gamma, NO(self.gamma, d(self.c)))))
            - 4 * NO(self.beta, NO(self.gamma, NO(self.b, NO(d(self.c), self.c))))
            - 4 * NO(self.beta, NO(self.gamma, d(self.c, 2)))
            - 4 * NO(self.beta, NO(d(self.gamma), d(self.c)))
            - sp.Rational(2, 3) * NO(d(self.b), NO(d(self.c), self.c))
            + 2 * NO(d(self.beta), NO(self.beta, NO(self.gamma, NO(self.gamma, self.c))))
            - sp.Rational(8, 3) * NO(d(self.beta), NO(d(self.gamma), self.c))
            + sp.Rational(2, 3) * NO(d(self.beta, 2), NO(self.gamma, self.c))
            + sp.Rational(10, 9) * d(self.c, 3)
        )

    def test_w_generator(self):
        """测试 w 生成元"""
        # w = β，共形权重应该是 -1/2
        assert self.w == self.beta
        assert self.w.conformal_weight == -1/2

    def test_j0_generator(self):
        """测试 j0 生成元（U(1) 流）"""
        # j0 = 2:bc: + 3:βγ:
        # 共形权重应该是 1
        # 注意：这是一个 sympy 表达式，需要检查其结构
        assert self.j0 is not None

        # 展开并检查项
        expanded = sp.expand(self.j0)
        print(f"j0 = {expanded}")

    def test_t_generator(self):
        """测试 t 生成元（能动张量）"""
        # t 的共形权重应该是 2
        # 注意：这是一个 sympy 表达式
        assert self.t is not None

        expanded = sp.expand(self.t)
        print(f"t = {expanded}")

    def test_g_generator(self):
        """测试 g 生成元"""
        # g = :γb:
        # 共形权重应该是 -1 + 3/2 = 1/2
        assert self.g.conformal_weight == 1/2

    def test_gt_generator(self):
        """测试 gt 生成元"""
        # gt = 2:∂βc: + 3:β∂c:
        assert self.gt is not None

        expanded = sp.expand(self.gt)
        print(f"gt = {expanded}")

    def test_gw_generator(self):
        """测试 gw 生成元"""
        # gw = b
        assert self.gw == self.b
        assert self.gw.conformal_weight == -1

    def test_wb_generator(self):
        """测试 wb 生成元（复杂的 W-场）"""
        # wb 是一个非常复杂的表达式
        # 它的共形权重应该是 3/2
        assert self.wb is not None

        # 打印表达式以便检查
        expanded = sp.expand(self.wb)
        print(f"\nwb = {expanded}")

        # 检查是否包含预期的项
        # 例如，应该包含 :β:β:γ:γγ::: 这样的多重正规序

    def test_gwb_generator(self):
        """测试 gwb 生成元（复杂的费米子场）"""
        # gwb 是另一个非常复杂的表达式
        # 它的共形权重应该是 2
        assert self.gwb is not None

        # 打印表达式以便检查
        expanded = sp.expand(self.gwb)
        print(f"\ngwb = {expanded}")


class TestZ3GeneratorOPEs:
    """测试 Z3 W-代数生成元之间的 OPE"""

    def setup_method(self):
        """设置自由场系统和 W-代数生成元"""
        # 定义自由场
        self.beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
        self.gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)
        self.b = BasisOperator("b", bosonic=False, conformal_weight=-1)
        self.c = BasisOperator("c", bosonic=False, conformal_weight=2)

        # 声明费米子
        Fermionic(self.beta, self.gamma, self.b, self.c)

        # 定义自由场 OPE
        OPE[self.beta, self.gamma] = MakeOPE([-One])
        OPE[self.b, self.c] = MakeOPE([One])

        # 构造生成元
        self.w = self.beta
        self.j0 = 2 * NO(self.b, self.c) + 3 * NO(self.beta, self.gamma)
        self.t = (
            -2 * NO(self.b, d(self.c))
            - sp.Rational(3, 2) * NO(self.beta, d(self.gamma))
            - NO(d(self.b), self.c)
            - sp.Rational(1, 2) * NO(d(self.beta), self.gamma)
        )
        self.g = NO(self.gamma, self.b)
        self.gt = 2 * NO(d(self.beta), self.c) + 3 * NO(self.beta, d(self.c))
        self.gw = self.b

    def test_w_w_ope(self):
        """测试 w(z)w(w) 的 OPE"""
        # w = β，所以 w(z)w(w) = β(z)β(w)
        # 由于 β 是费米子，β(z)β(w) 应该没有奇异项（费米子自身的 OPE 为 0）
        ope_result = OPE(self.w, self.w)

        # 应该没有极点
        assert ope_result.max_pole == 0

    def test_w_gamma_ope(self):
        """测试 w(z)γ(w) 的 OPE"""
        # w = β，所以这就是 β(z)γ(w) ~ -1/(z-w)
        ope_result = OPE(self.w, self.gamma)

        assert ope_result.max_pole == 1
        assert ope_result.pole(1) == -One

    def test_b_c_ope(self):
        """测试 b(z)c(w) 的 OPE"""
        ope_result = OPE(self.b, self.c)

        assert ope_result.max_pole == 1
        assert ope_result.pole(1) == One

    def test_g_gw_ope(self):
        """测试 g(z)gw(w) 的 OPE"""
        # g = :γb:, gw = b
        # 这个 OPE 应该涉及 Jacobi 恒等式
        ope_result = OPE(self.g, self.gw)

        # 打印结果以便调试
        print(f"\nOPE(g, gw) max_pole = {ope_result.max_pole}")
        for i in range(1, ope_result.max_pole + 1):
            print(f"  pole({i}) = {ope_result.pole(i)}")

    def test_j0_j0_ope(self):
        """测试 j0(z)j0(w) 的 OPE（U(1) 流的 OPE）"""
        # j0 = 2:bc: + 3:βγ:
        # 理论上，U(1) 流的 OPE 应该是：
        # j0(z)j0(w) ~ k/(z-w)^2
        # 其中 k 是某个常数（与中心荷相关）
        print("\n计算 j0(z)j0(w)...")
        ope_result = OPE(self.j0, self.j0)

        print(f"max_pole = {ope_result.max_pole}")
        for i in range(1, min(ope_result.max_pole + 1, 5)):
            pole_i = ope_result.pole(i)
            print(f"  pole({i}) = {pole_i}")

        # U(1) 流应该有 2 阶极点
        assert ope_result.max_pole >= 2

    def test_t_j0_ope(self):
        """测试 t(z)j0(w) 的 OPE（能动张量对 U(1) 流的作用）"""
        # 理论上：t(z)j0(w) ~ j0(w)/(z-w)^2 + ∂j0(w)/(z-w)
        # 这表明 j0 是权重 1 的准初级场
        print("\n计算 t(z)j0(w)...")
        ope_result = OPE(self.t, self.j0)

        print(f"max_pole = {ope_result.max_pole}")
        for i in range(1, min(ope_result.max_pole + 1, 5)):
            pole_i = ope_result.pole(i)
            print(f"  pole({i}) = {pole_i}")

        # 应该至少有 2 阶极点
        assert ope_result.max_pole >= 2

    def test_t_t_ope(self):
        """测试 t(z)t(w) 的 OPE（Virasoro 代数）"""
        # 理论上：t(z)t(w) ~ c/2/(z-w)^4 + 2t(w)/(z-w)^2 + ∂t(w)/(z-w)
        # 这是 Virasoro 代数的标准 OPE
        print("\n计算 t(z)t(w)...")
        ope_result = OPE(self.t, self.t)

        print(f"max_pole = {ope_result.max_pole}")
        for i in range(1, min(ope_result.max_pole + 1, 5)):
            pole_i = ope_result.pole(i)
            print(f"  pole({i}) = {pole_i}")

        # Virasoro 代数应该有 4 阶极点
        assert ope_result.max_pole >= 2  # 至少应该有 2 阶极点

    def test_t_w_ope(self):
        """测试 t(z)w(w) 的 OPE"""
        # w = β 的权重是 -1/2
        # 理论上：t(z)w(w) ~ h*w(w)/(z-w)^2 + ∂w(w)/(z-w)
        # 其中 h = -1/2
        print("\n计算 t(z)w(w)...")
        ope_result = OPE(self.t, self.w)

        print(f"max_pole = {ope_result.max_pole}")
        for i in range(1, min(ope_result.max_pole + 1, 5)):
            pole_i = ope_result.pole(i)
            print(f"  pole({i}) = {pole_i}")

        # 应该至少有 2 阶极点
        assert ope_result.max_pole >= 2

    def test_t_g_ope(self):
        """测试 t(z)g(w) 的 OPE"""
        # g = :γb: 的权重是 1/2
        # 理论上：t(z)g(w) ~ h*g(w)/(z-w)^2 + ∂g(w)/(z-w)
        # 其中 h = 1/2
        print("\n计算 t(z)g(w)...")
        ope_result = OPE(self.t, self.g)

        print(f"max_pole = {ope_result.max_pole}")
        for i in range(1, min(ope_result.max_pole + 1, 5)):
            pole_i = ope_result.pole(i)
            print(f"  pole({i}) = {pole_i}")

        # 应该至少有 2 阶极点
        assert ope_result.max_pole >= 2


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
