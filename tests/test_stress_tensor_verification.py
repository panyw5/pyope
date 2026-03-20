"""
验证 bc-βγ 系统应力张量定义的正确性

测试应力张量：T = -λ NO(b, ∂c) + (1-λ) NO(∂b, c)

验证：
1. OPE(T, b) = λb/(z-w)^2 + ∂b/(z-w)
2. OPE(T, c) = (1-λ)c/(z-w)^2 + ∂c/(z-w)

这确保 b, c 是 Virasoro primary fields，且给出正确的 conformal weights
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


def test_bc_stress_tensor_definition_lambda_2():
    """
    验证 bc 系统应力张量定义（λ=2）

    T = -λ NO(b, ∂c) + (1-λ) NO(∂b, c)
      = -2 NO(b, ∂c) + (-1) NO(∂b, c)
      = -2 NO(b, ∂c) - NO(∂b, c)

    验证：
    1. OPE(T, b) = 2b/(z-w)^2 + ∂b/(z-w)  (b 是 weight 2 的 primary)
    2. OPE(T, c) = -c/(z-w)^2 + ∂c/(z-w)  (c 是 weight -1 的 primary)
    """
    lam_val = 2

    b = BasisOperator("b", conformal_weight=lam_val, fermionic=True)
    c = BasisOperator("c", conformal_weight=1 - lam_val, fermionic=True)
    Fermionic(b, c)

    OPE[b, c] = MakeOPE([One])

    # 正确的应力张量定义
    T = -lam_val * NO(b, d(c)) + (1 - lam_val) * NO(d(b), c)

    print(f"\n应力张量: T = -λ NO(b, ∂c) + (1-λ) NO(∂b, c)")
    print(f"         = -{lam_val} NO(b, ∂c) + {1 - lam_val} NO(∂b, c)")

    # 测试 OPE(T, b)
    result_T_b = OPE(T, b)

    print(f"\n测试 OPE(T, b):")
    print(f"  max_pole = {result_T_b.max_pole}")
    print(f"  pole(2) = {result_T_b.pole(2)}")
    print(f"  pole(1) = {result_T_b.pole(1)}")

    # 验证 pole(2) = λ*b = 2*b
    pole_2_expected = lam_val * b
    pole_2_actual = result_T_b.pole(2)
    print(f"  pole(2) 期望值 = {pole_2_expected}")
    print(f"  pole(2) 实际值 = {pole_2_actual}")
    # 使用字符串比较避免 sympy 简化问题
    assert str(pole_2_actual) == str(pole_2_expected), f"pole(2) 应该是 {lam_val}*b"

    # 验证 pole(1) = ∂b
    pole_1_expected = d(b)
    pole_1_actual = result_T_b.pole(1)
    print(f"  pole(1) 期望值 = {pole_1_expected}")
    print(f"  pole(1) 实际值 = {pole_1_actual}")
    assert str(pole_1_actual) == str(pole_1_expected), "pole(1) 应该是 ∂b"

    print(f"✓ OPE(T, b) = {lam_val}b/(z-w)^2 + ∂b/(z-w) ✓")

    # 测试 OPE(T, c)
    result_T_c = OPE(T, c)

    print(f"\n测试 OPE(T, c):")
    print(f"  max_pole = {result_T_c.max_pole}")
    print(f"  pole(2) = {result_T_c.pole(2)}")
    print(f"  pole(1) = {result_T_c.pole(1)}")

    # 验证 pole(2) = (1-λ)*c = -1*c
    pole_2_expected = (1 - lam_val) * c
    pole_2_actual = result_T_c.pole(2)
    print(f"  pole(2) 期望值 = {pole_2_expected}")
    print(f"  pole(2) 实际值 = {pole_2_actual}")
    assert str(pole_2_actual) == str(pole_2_expected), f"pole(2) 应该是 {1 - lam_val}*c"

    # 验证 pole(1) = ∂c
    pole_1_expected = d(c)
    pole_1_actual = result_T_c.pole(1)
    print(f"  pole(1) 期望值 = {pole_1_expected}")
    print(f"  pole(1) 实际值 = {pole_1_actual}")
    assert str(pole_1_actual) == str(pole_1_expected), "pole(1) 应该是 ∂c"

    print(f"✓ OPE(T, c) = {1 - lam_val}c/(z-w)^2 + ∂c/(z-w) ✓")

    print("\n✅ 应力张量定义正确！b 和 c 都是 Virasoro primary fields")


def test_betagamma_stress_tensor_definition_lambda_3_2():
    """
    验证 βγ 系统应力张量定义（λ=3/2）

    T = -λ NO(β, ∂γ) + (1-λ) NO(∂β, γ)
      = -3/2 NO(β, ∂γ) + (-1/2) NO(∂β, γ)
      = -3/2 NO(β, ∂γ) - 1/2 NO(∂β, γ)

    验证：
    1. OPE(T, β) = (3/2)β/(z-w)^2 + ∂β/(z-w)  (β 是 weight 3/2 的 primary)
    2. OPE(T, γ) = (-1/2)γ/(z-w)^2 + ∂γ/(z-w)  (γ 是 weight -1/2 的 primary)
    """
    lam_val = Rational(3, 2)

    # βγ twisted ghost 系统是玻色子（与 tests/test_bc_betagamma.py 的约定一致）
    beta = BasisOperator("β", conformal_weight=lam_val, fermionic=False)
    gamma = BasisOperator("γ", conformal_weight=1 - lam_val, fermionic=False)

    OPE[beta, gamma] = MakeOPE([-One])

    # 正确的应力张量定义
    T = -lam_val * NO(beta, d(gamma)) + (1 - lam_val) * NO(d(beta), gamma)

    print(f"\n应力张量: T = -λ NO(β, ∂γ) + (1-λ) NO(∂β, γ)")
    print(f"         = -{lam_val} NO(β, ∂γ) + {1 - lam_val} NO(∂β, γ)")

    # 测试 OPE(T, β)
    result_T_beta = OPE(T, beta)

    print(f"\n测试 OPE(T, β):")
    print(f"  max_pole = {result_T_beta.max_pole}")
    print(f"  pole(2) = {result_T_beta.pole(2)}")
    print(f"  pole(1) = {result_T_beta.pole(1)}")

    # 验证 pole(2) = λ*β = (3/2)*β
    pole_2_expected = lam_val * beta
    pole_2_actual = result_T_beta.pole(2)
    print(f"  pole(2) 期望值 = {pole_2_expected}")
    print(f"  pole(2) 实际值 = {pole_2_actual}")
    assert str(pole_2_actual) == str(pole_2_expected), f"pole(2) 应该是 {lam_val}*β"

    # 验证 pole(1) = ∂β
    pole_1_expected = d(beta)
    pole_1_actual = result_T_beta.pole(1)
    print(f"  pole(1) 期望值 = {pole_1_expected}")
    print(f"  pole(1) 实际值 = {pole_1_actual}")
    assert str(pole_1_actual) == str(pole_1_expected), "pole(1) 应该是 ∂β"

    print(f"✓ OPE(T, β) = {lam_val}β/(z-w)^2 + ∂β/(z-w) ✓")

    # 测试 OPE(T, γ)
    result_T_gamma = OPE(T, gamma)

    print(f"\n测试 OPE(T, γ):")
    print(f"  max_pole = {result_T_gamma.max_pole}")
    print(f"  pole(2) = {result_T_gamma.pole(2)}")
    print(f"  pole(1) = {result_T_gamma.pole(1)}")

    # 验证 pole(2) = (1-λ)*γ = (-1/2)*γ
    pole_2_expected = (1 - lam_val) * gamma
    pole_2_actual = result_T_gamma.pole(2)
    print(f"  pole(2) 期望值 = {pole_2_expected}")
    print(f"  pole(2) 实际值 = {pole_2_actual}")
    assert str(pole_2_actual) == str(pole_2_expected), f"pole(2) 应该是 {1 - lam_val}*γ"

    # 验证 pole(1) = ∂γ
    pole_1_expected = d(gamma)
    pole_1_actual = result_T_gamma.pole(1)
    print(f"  pole(1) 期望值 = {pole_1_expected}")
    print(f"  pole(1) 实际值 = {pole_1_actual}")
    assert str(pole_1_actual) == str(pole_1_expected), "pole(1) 应该是 ∂γ"

    print(f"✓ OPE(T, γ) = {1 - lam_val}γ/(z-w)^2 + ∂γ/(z-w) ✓")

    print("\n✅ 应力张量定义正确！β 和 γ 都是 Virasoro primary fields")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
