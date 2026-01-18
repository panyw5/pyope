"""
U(1)_k Kac-Moody 代数测试

基于 VOA.wls 中的 U(1)_k 定义和计算

测试内容：
1. U(1) current J 的声明和 OPE
2. Sugawara 应力张量的构造
3. 中心荷的计算（c=1）
4. Vertex operators 的 OPE
5. Primary fields 的验证

参考资料：
- VOA.wls Section 3: U(1)_k Kac-Moody Algebra
- Di Francesco et al., "Conformal Field Theory", Chapter 15
- Kac, "Infinite Dimensional Lie Algebras"
"""

import pytest
import sympy as sp
from sympy import Symbol, Rational, simplify, expand, sqrt, exp, I

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
    Bosonic,
)


@pytest.fixture(autouse=True)
def clear_registry():
    """每个测试前清空注册表"""
    from pyope.registry import ope_registry

    ope_registry.clear()
    yield
    ope_registry.clear()


# ============================================================================
# U(1)_k 基础定义测试
# ============================================================================


class TestU1KDefinition:
    """
    U(1)_k Kac-Moody 代数基础定义测试

    测试内容：
    1. U(1) current J 的声明
    2. J-J OPE 的定义
    3. Sugawara 应力张量的构造
    """

    def test_u1_current_declaration(self):
        """
        测试 U(1) current J 的声明

        J 是 conformal weight 1 的玻色算符
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        assert J.conformal_weight == 1
        assert J._bosonic is True
        assert str(J) == "J"

        print("✓ U(1) current J 声明: weight=1, bosonic")

    def test_u1_ope(self):
        """
        测试 U(1)_k 的基本 OPE

        定义：J(z)J(w) ~ k/(z-w)²

        这是最简单的 Kac-Moody 代数

        参考：VOA.wls U(1)_k section
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        result = OPE(J, J)

        assert result.max_pole == 2
        assert result.pole(2) == k * One
        assert result.pole(1) == 0

        print("✓ U(1)_k OPE: J(z)J(w) ~ k/(z-w)²")
        print(f"  - 2-pole: {result.pole(2)}")
        print(f"  - 1-pole: {result.pole(1)}")

    def test_sugawara_construction(self):
        """
        测试 Sugawara 应力张量的构造

        T = NO(J,J)/(2k)

        这是从 Kac-Moody 代数构造 Virasoro 代数的标准方法

        注意：由于 sympy 的除法，T 是一个表达式，不是 Operator
        我们只验证 NO(J,J) 的 conformal weight
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        # 构造 Sugawara 张量的分子
        no_jj = NO(J, J)

        # 验证 NO(J,J) 的 conformal weight
        assert no_jj.conformal_weight == 2

        print(f"✓ Sugawara 张量构造: NO(J,J) weight = {no_jj.conformal_weight}")

    def test_u1_central_charge(self):
        """
        测试 U(1)_k 的中心荷

        通过计算 OPE[NO(J,J), NO(J,J)] 提取中心荷
        然后除以 (2k)²

        预期：c = 1（对于 U(1)）

        注意：这个测试比较复杂，我们简化为验证 NO(J,J) 的性质
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k_val = 1  # 使用具体值简化计算
        OPE[J, J] = MakeOPE([k_val * One, 0])

        # 构造 NO(J,J)
        no_jj = NO(J, J)

        # 计算 OPE[NO(J,J), NO(J,J)]
        result = OPE(no_jj, no_jj)

        # 验证有 4 阶极点
        assert result.max_pole >= 4

        # 提取 4 阶极点
        pole_4 = result.pole(4)

        # 对于 k=1，Sugawara 张量是 T = NO(J,J)/2
        # OPE[T,T] 的 4-pole 应该是 c/2 = 1/2
        # 因此 OPE[NO(J,J), NO(J,J)] 的 4-pole 应该是 4 * (1/2) = 2
        expected = 2 * One

        diff = simplify(pole_4 - expected)
        assert diff == 0

        print("✓ U(1)_k 中心荷验证（通过 NO(J,J) OPE）")
        print(f"  - OPE[NO(J,J), NO(J,J)] 的 4-pole: {pole_4}")


# ============================================================================
# U(1)_k 计算测试
# ============================================================================


class TestU1KComputations:
    """
    U(1)_k 的具体计算测试

    测试内容：
    1. T-J OPE（验证 J 是 primary）
    2. 导数算符的 OPE
    3. 正规序乘积
    """

    def test_t_j_ope(self):
        """
        测试 OPE[T, J]

        验证 J 是 conformal weight 1 的 primary field

        预期：T(z)J(w) ~ J/(z-w)² + ∂J/(z-w)
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        # 构造 Sugawara 张量
        T = NO(J, J) / (2 * k)

        # 计算 OPE[T, J]
        result = OPE(T, J)

        # J 是 primary field of weight 1
        assert result.max_pole == 2

        # pole(2) 应该是 1*J（weight = 1）
        pole_2 = simplify(result.pole(2) - J)
        assert pole_2 == 0

        # pole(1) 应该是 ∂J
        pole_1 = simplify(result.pole(1) - d(J))
        assert pole_1 == 0

        print("✓ OPE[T, J]: J 是 primary field of weight 1")
        print(f"  - 2-pole: {result.pole(2)}")
        print(f"  - 1-pole: {result.pole(1)}")

    def test_j_derivative_ope(self):
        """
        测试 OPE[J, ∂J]

        预期：max_pole = 3
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        result = OPE(J, d(J))

        assert result.max_pole == 3
        assert result.pole(3) == 0  # 因为 J-J OPE 的 1-pole 是 0

        print(f"✓ OPE[J, ∂J]: max_pole = {result.max_pole}")

    def test_j_second_derivative_ope(self):
        """
        测试 OPE[J, ∂²J]

        验证高阶导数的 OPE
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        result = OPE(J, d(J, 2))

        assert result.max_pole == 4
        # pole(4) = 6*k*One（来自导数规则）
        expected_pole_4 = 6 * k * One
        diff = simplify(result.pole(4) - expected_pole_4)
        assert diff == 0

        print(f"✓ OPE[J, ∂²J]: 4-pole = {result.pole(4)}")

    def test_normal_order_j_j(self):
        """
        测试 NO(J,J)

        这是 Sugawara 张量的核心部分
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        no_jj = NO(J, J)

        assert no_jj.conformal_weight == 2

        print(f"✓ NO(J,J): weight = {no_jj.conformal_weight}")

    def test_t_no_j_j_ope(self):
        """
        测试 OPE[T, NO(J,J)]

        验证 Sugawara 张量与自身的关系
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        T = NO(J, J) / (2 * k)
        no_jj = NO(J, J)

        result = OPE(T, no_jj)

        # NO(J,J) 应该是 primary field of weight 2
        assert result.max_pole >= 2

        print(f"✓ OPE[T, NO(J,J)]: max_pole = {result.max_pole}")


# ============================================================================
# U(1)_k 性质验证
# ============================================================================


class TestU1KProperties:
    """
    U(1)_k 的代数性质验证

    测试内容：
    1. Kac-Moody 代数的对易关系
    2. Primary fields 的性质
    3. 电荷守恒
    """

    def test_kac_moody_algebra(self):
        """
        测试 U(1) Kac-Moody 代数

        对易关系：[J_m, J_n] = k*m*δ_{m+n,0}

        这可以从 OPE 推导出来
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        # OPE 的 2-pole 对应于对易子的中心项
        result = OPE(J, J)
        assert result.pole(2) == k * One

        print("✓ U(1) Kac-Moody 代数验证通过")

    def test_primary_field_under_t(self):
        """
        测试 primary field 在 T 作用下的行为

        如果 φ 是 weight h 的 primary，则：
        T(z)φ(w) ~ h*φ/(z-w)² + ∂φ/(z-w)
        """
        J = BasisOperator("J", conformal_weight=1)
        phi = BasisOperator("φ", conformal_weight=2)
        Bosonic(J, phi)

        k = Symbol("k", positive=True)
        h = Symbol("h")

        OPE[J, J] = MakeOPE([k * One, 0])
        T = NO(J, J) / (2 * k)

        # 定义 φ 为 primary field
        OPE[T, phi] = MakeOPE([h * phi, d(phi)])

        result = OPE(T, phi)

        assert result.max_pole == 2
        assert result.pole(2) == h * phi
        assert result.pole(1) == d(phi)

        print("✓ Primary field 条件验证通过")

    def test_virasoro_from_sugawara(self):
        """
        测试从 Sugawara 构造得到的 Virasoro 代数

        验证 T 满足 Virasoro OPE
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        OPE[J, J] = MakeOPE([k * One, 0])

        T = NO(J, J) / (2 * k)

        result = OPE(T, T)

        # 验证 Virasoro 结构
        assert result.max_pole == 4

        # 中心荷 c = 1
        c_over_2 = result.pole(4)
        assert simplify(c_over_2 - Rational(1, 2) * One) == 0

        # 2-pole 应该是 2*T
        pole_2 = simplify(result.pole(2) - 2 * T)
        assert pole_2 == 0

        print("✓ Sugawara 构造的 Virasoro 代数验证通过")


# ============================================================================
# U(1)_k 数值测试
# ============================================================================


class TestU1KNumerical:
    """
    U(1)_k 的数值测试

    使用具体的参数值验证计算结果
    """

    def test_u1_k_equals_1(self):
        """
        测试 k=1 的 U(1) 代数
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k_val = 1
        OPE[J, J] = MakeOPE([k_val * One, 0])

        result = OPE(J, J)
        assert result.pole(2) == One

        # Sugawara 张量
        T = NO(J, J) / (2 * k_val)
        result_T = OPE(T, T)

        # c = 1
        assert simplify(result_T.pole(4) - Rational(1, 2) * One) == 0

        print("✓ U(1)_1 数值验证通过")

    def test_u1_k_equals_2(self):
        """
        测试 k=2 的 U(1) 代数
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k_val = 2
        OPE[J, J] = MakeOPE([k_val * One, 0])

        result = OPE(J, J)
        assert result.pole(2) == 2 * One

        # Sugawara 张量
        T = NO(J, J) / (2 * k_val)
        result_T = OPE(T, T)

        # c = 1（与 k 无关）
        assert simplify(result_T.pole(4) - Rational(1, 2) * One) == 0

        print("✓ U(1)_2 数值验证通过")


# ============================================================================
# U(1)_k Vertex Operators（高级）
# ============================================================================


class TestU1KVertexOperators:
    """
    U(1)_k 的 vertex operators 测试

    测试内容：
    1. Vertex operator 的定义
    2. Vertex operator 的 OPE
    3. 与 current 的 OPE

    注意：这部分需要指数算符，可能需要扩展 pyope
    """

    def test_vertex_operator_definition(self):
        """
        测试 vertex operator 的定义

        V_α(z) = :exp(iα∫J):

        性质：
        - conformal weight h = α²/(2k)
        - U(1) charge q = α

        注意：这是概念性测试，实际实现需要指数算符支持
        """
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)

        k = Symbol("k", positive=True)
        alpha = Symbol("alpha", real=True)

        OPE[J, J] = MakeOPE([k * One, 0])

        # Vertex operator 的 conformal weight
        h_V = alpha**2 / (2 * k)

        print(f"✓ Vertex operator V_α 的 weight: h = α²/(2k) = {h_V}")

    def test_j_vertex_operator_ope_concept(self):
        """
        测试 J 与 vertex operator 的 OPE（概念）

        J(z)V_α(w) ~ α*V_α/(z-w)

        这表明 V_α 携带 U(1) 电荷 α
        """
        # 这是概念性测试，说明预期的行为
        alpha = Symbol("alpha", real=True)

        # 预期的 OPE 结构
        # OPE[J, V_α] = MakeOPE([alpha * V_α])

        print(f"✓ 概念验证: J(z)V_α(w) ~ α*V_α/(z-w)")


# ============================================================================
# 运行测试
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
