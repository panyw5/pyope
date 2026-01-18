"""
bc-βγ 自由场系统完整测试

基于 VOA.wls 中的 bc 和 βγ 系统定义

测试内容：
1. bc 系统（ghost 系统）
2. βγ 系统（twisted ghost 系统）
3. 应力张量的构造
4. 中心荷的计算
5. bc 和 βγ 的对偶关系

参考资料：
- VOA.wls Section 2: bc and βγ Systems
- Polchinski, "String Theory", Volume 1, Chapter 2
"""

import pytest
import sympy as sp
from sympy import Symbol, Rational, simplify, expand

from pyope import (
    BasisOperator,
    d,
    dn,
    One,
    Zero,
    OPE,
    NO,
    bracket,
    MakeOPE,
    Fermionic,
)


@pytest.fixture(autouse=True)
def clear_registry():
    """每个测试前清空注册表"""
    from pyope.registry import ope_registry

    ope_registry.clear()
    yield
    ope_registry.clear()


# ============================================================================
# bc 系统测试 (Fermionic)
# ============================================================================

class TestBCSystem:
    """
    bc ghost 系统测试
    """

    def test_bc_operators_declaration(self):
        lam_val = 2
        b = BasisOperator("b", conformal_weight=lam_val, fermionic=True)
        c = BasisOperator("c", conformal_weight=1 - lam_val, fermionic=True)
        Fermionic(b, c)

        assert b.conformal_weight == lam_val
        assert c.conformal_weight == 1 - lam_val
        assert b.is_bosonic is False
        assert c.is_bosonic is False

    def test_bc_ope(self):
        b = BasisOperator("b", conformal_weight=2, fermionic=True)
        c = BasisOperator("c", conformal_weight=-1, fermionic=True)
        Fermionic(b, c)
        OPE[b, c] = MakeOPE([One])

        result = OPE(b, c)
        assert result.max_pole == 1
        assert result.pole(1) == One

    def test_bc_central_charge(self):
        lam_val = 2
        b = BasisOperator("b", conformal_weight=lam_val, fermionic=True)
        c = BasisOperator("c", conformal_weight=1 - lam_val, fermionic=True)
        Fermionic(b, c)
        OPE[b, c] = MakeOPE([One])

        # 正确的 T 形式: T = -λ NO(b, ∂c) + (1-λ) NO(∂b, c)
        T_bc = -lam_val * NO(b, d(c)) + (1 - lam_val) * NO(d(b), c)

        result = OPE(T_bc, T_bc)
        c_over_2 = result.pole(4)
        
        # λ=2 时, c = -2(6*4 - 6*2 + 1) = -2(13) = -26, 所以 c/2 = -13
        expected = -13 * One
        assert simplify(c_over_2 - expected) == 0

    def test_bc_primary_weights(self):
        lam_val = 2
        b = BasisOperator("b", conformal_weight=lam_val, fermionic=True)
        c = BasisOperator("c", conformal_weight=1 - lam_val, fermionic=True)
        Fermionic(b, c)
        OPE[b, c] = MakeOPE([One])

        T_bc = -lam_val * NO(b, d(c)) + (1 - lam_val) * NO(d(b), c)

        # 验证 b 的权重
        res_b = OPE(T_bc, b)
        assert res_b.pole(2) == lam_val * b
        assert res_b.pole(1) == d(b)

        # 验证 c 的权重
        res_c = OPE(T_bc, c)
        assert res_c.pole(2) == (1 - lam_val) * c
        assert res_c.pole(1) == d(c)


# ============================================================================
# βγ 系统测试 (Bosonic)
# ============================================================================

class TestBetaGammaSystem:
    """
    βγ twisted ghost 系统测试 (玻色子)
    """

    def test_betagamma_ope(self):
        lam_val = Rational(3, 2)
        # βγ 系统是玻色子
        beta = BasisOperator("β", conformal_weight=lam_val, fermionic=False)
        gamma = BasisOperator("γ", conformal_weight=1 - lam_val, fermionic=False)
        
        # OPE 负号约定 (与 VOA.wls 一致)
        OPE[beta, gamma] = MakeOPE([-One])

        result = OPE(beta, gamma)
        assert result.pole(1) == -One

    def test_betagamma_central_charge(self):
        lam_val = Rational(3, 2)
        beta = BasisOperator("β", conformal_weight=lam_val, fermionic=False)
        gamma = BasisOperator("γ", conformal_weight=1 - lam_val, fermionic=False)
        OPE[beta, gamma] = MakeOPE([-One])

        # T = -λ NO(β, ∂γ) + (1-λ) NO(∂β, γ)
        T_bg = -lam_val * NO(beta, d(gamma)) + (1 - lam_val) * NO(d(beta), gamma)

        result = OPE(T_bg, T_bg)
        c_over_2 = result.pole(4)

        # λ=3/2 时, c = 2(6*(9/4) - 6*(3/2) + 1) = 2(13.5 - 9 + 1) = 11, 所以 c/2 = 11/2
        expected = Rational(11, 2) * One
        assert simplify(c_over_2 - expected) == 0

    def test_betagamma_primary_weights(self):
        lam_val = Rational(3, 2)
        beta = BasisOperator("β", conformal_weight=lam_val, fermionic=False)
        gamma = BasisOperator("γ", conformal_weight=1 - lam_val, fermionic=False)
        OPE[beta, gamma] = MakeOPE([-One])

        T_bg = -lam_val * NO(beta, d(gamma)) + (1 - lam_val) * NO(d(beta), gamma)

        # 验证 β 的权重 λ
        res_beta = OPE(T_bg, beta)
        assert res_beta.pole(2) == lam_val * beta
        assert res_beta.pole(1) == d(beta)

        # 验证 γ 的权重 1-λ
        res_gamma = OPE(T_bg, gamma)
        assert res_gamma.pole(2) == (1 - lam_val) * gamma
        assert res_gamma.pole(1) == d(gamma)

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
