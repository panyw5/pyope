"""
测试自由场系统（bc 和 βγ）的 OPE

根据 null_states.txt 中的定义：
- OPE[β1,γ1] = MakeOPE[-One (z-w)^-1+Ord[z,w,0]]
- OPE[b1,c1] = MakeOPE[One (z-w)^-1+Ord[z,w,0]]

这意味着：
- β(z)γ(w) ~ -1/(z-w)  （βγ 系统，费米子）
- b(z)c(w) ~ 1/(z-w)   （bc 系统，费米子）
"""

import sys
sys.path.insert(0, '../src')

import pytest
from pyope import (
    BasisOperator,
    OPE,
    MakeOPE,
    One,
    Zero,
    NO,
    bracket,
    Fermionic,
)


class TestBetaGammaSystem:
    """测试 βγ 自由场系统"""

    def setup_method(self):
        """设置 βγ 系统"""
        # 定义 β 和 γ 算符（费米子）
        # β 的共形权重是 -1/2，γ 的共形权重是 3/2
        self.beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
        self.gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)

        # 声明它们是费米子
        Fermionic(self.beta, self.gamma)

        # 定义 OPE: β(z)γ(w) ~ -1/(z-w)
        OPE[self.beta, self.gamma] = MakeOPE([-One])

    def test_beta_gamma_ope(self):
        """测试 β(z)γ(w) 的 OPE"""
        ope_result = OPE(self.beta, self.gamma)

        # 检查只有一个极点（阶数为 1）
        assert ope_result.max_pole == 1

        # 检查极点系数
        assert ope_result.pole(1) == -One
        assert ope_result.pole(2) == 0

    def test_beta_gamma_normal_order(self):
        """测试 βγ 的正规序"""
        # NO(β, γ) = {βγ}_0 = 0（因为 OPE 只有 1 阶极点）
        no_result = NO(self.beta, self.gamma)

        # 正规序应该返回 NormalOrderedOperator
        from pyope.operators import NormalOrderedOperator
        assert isinstance(no_result, NormalOrderedOperator)

    def test_gamma_beta_ope(self):
        """测试 γ(z)β(w) 的 OPE（反向）"""
        # 由于 β 和 γ 都是费米子，应该有反对称性
        # 根据 pyope 的算符交换公式，γ(z)β(w) 应该从 β(z)γ(w) 推导
        # β(z)γ(w) ~ -1/(z-w)
        # 因此 γ(z)β(w) ~ -1/(z-w)（相同符号，因为费米子交换两次）
        ope_result = OPE(self.gamma, self.beta)

        # 应该有一个极点
        assert ope_result.max_pole == 1

        # 系数应该是 -1（与 β(z)γ(w) 相同）
        assert ope_result.pole(1) == -One


class TestBCSystem:
    """测试 bc 自由场系统"""

    def setup_method(self):
        """设置 bc 系统"""
        # 定义 b 和 c 算符（费米子）
        # b 的共形权重是 -1，c 的共形权重是 2
        self.b = BasisOperator("b", bosonic=False, conformal_weight=-1)
        self.c = BasisOperator("c", bosonic=False, conformal_weight=2)

        # 声明它们是费米子
        Fermionic(self.b, self.c)

        # 定义 OPE: b(z)c(w) ~ 1/(z-w)
        OPE[self.b, self.c] = MakeOPE([One])

    def test_bc_ope(self):
        """测试 b(z)c(w) 的 OPE"""
        ope_result = OPE(self.b, self.c)

        # 检查只有一个极点（阶数为 1）
        assert ope_result.max_pole == 1

        # 检查极点系数
        assert ope_result.pole(1) == One
        assert ope_result.pole(2) == 0

    def test_bc_normal_order(self):
        """测试 bc 的正规序"""
        # NO(b, c) = {bc}_0 = 0（因为 OPE 只有 1 阶极点）
        no_result = NO(self.b, self.c)

        # 正规序应该返回 NormalOrderedOperator
        from pyope.operators import NormalOrderedOperator
        assert isinstance(no_result, NormalOrderedOperator)

    def test_cb_ope(self):
        """测试 c(z)b(w) 的 OPE（反向）"""
        # 由于 b 和 c 都是费米子，应该有反对称性
        # b(z)c(w) ~ 1/(z-w)
        # 因此 c(z)b(w) ~ 1/(z-w)（相同符号）
        ope_result = OPE(self.c, self.b)

        # 应该有一个极点
        assert ope_result.max_pole == 1

        # 系数应该是 +1（与 b(z)c(w) 相同）
        assert ope_result.pole(1) == One


class TestFreeFieldNormalOrders:
    """测试自由场的正规序乘积"""

    def setup_method(self):
        """设置自由场系统"""
        # βγ 系统
        self.beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
        self.gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)
        Fermionic(self.beta, self.gamma)
        OPE[self.beta, self.gamma] = MakeOPE([-One])

        # bc 系统
        self.b = BasisOperator("b", bosonic=False, conformal_weight=-1)
        self.c = BasisOperator("c", bosonic=False, conformal_weight=2)
        Fermionic(self.b, self.c)
        OPE[self.b, self.c] = MakeOPE([One])

    def test_normal_order_bc(self):
        """测试 :bc: 的性质"""
        no_bc = NO(self.b, self.c)

        # 共形权重应该是 b 和 c 的权重之和
        # weight(:bc:) = weight(b) + weight(c) = -1 + 2 = 1
        assert no_bc.conformal_weight == 1

    def test_normal_order_beta_gamma(self):
        """测试 :βγ: 的性质"""
        no_bg = NO(self.beta, self.gamma)

        # 共形权重应该是 β 和 γ 的权重之和
        # weight(:βγ:) = weight(β) + weight(γ) = -1/2 + 3/2 = 1
        assert no_bg.conformal_weight == 1

    def test_double_normal_order(self):
        """测试双重正规序 :b:cγ::"""
        # 这在 W-代数构造中很常见
        no_bc = NO(self.b, self.c)
        no_bc_gamma = NO(no_bc, self.gamma)

        # 应该能够构造
        assert no_bc_gamma is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
