"""
测试：确保系统中不会出现非法的 operator * operator 表达式

这个测试模块的目的是检测系统中是否会产生以下非法表达式：
1. operator * One（如 b*One）
2. operator1 * operator2（如 b*c，两个算符的普通乘积）
3. scalar * One * operator（如 2*One*b）

合法的表达式：
- scalar * One（如 k*One，纯标量）
- scalar * operator（如 2*b）
- NO(operator1, operator2)（正规序乘积）
"""

import pytest
import sympy as sp
from pyope import BasisOperator, OPE, NO, d, simplify, Bosonic, Fermionic, One, Zero
from pyope.operators import Operator
from pyope.constants import ConstantOperator


def check_for_illegal_operator_products(expr, path="", depth=0, max_depth=20):
    """
    递归检查表达式中是否包含非法的 operator * operator 组合

    Args:
        expr: 要检查的表达式
        path: 当前路径（用于错误报告）
        depth: 当前递归深度
        max_depth: 最大递归深度

    Returns:
        (has_illegal, error_messages): 是否有非法表达式，以及错误信息列表
    """
    if depth > max_depth:
        return False, []

    errors = []

    if isinstance(expr, sp.Mul):
        # 检查 Mul 的参数中是否有 One 和非常数算符同时出现
        has_one = False
        has_non_const_operator = False
        operators = []

        for arg in expr.args:
            if arg is One or arg == One:
                has_one = True
            if isinstance(arg, Operator):
                if not isinstance(arg, ConstantOperator):
                    has_non_const_operator = True
                    operators.append(arg)

        # 检查是否有 One 和非常数算符同时出现
        if has_one and has_non_const_operator:
            errors.append(f"非法表达式 (One * Operator) 在 {path}: {expr}")
            errors.append(f"  args: {expr.args}")

        # 检查是否有多个非常数算符（operator1 * operator2）
        if len(operators) >= 2:
            errors.append(f"非法表达式 (Operator * Operator) 在 {path}: {expr}")
            errors.append(f"  算符: {operators}")

        # 递归检查子表达式
        for i, arg in enumerate(expr.args):
            sub_has_illegal, sub_errors = check_for_illegal_operator_products(
                arg, f"{path}.args[{i}]", depth + 1, max_depth
            )
            errors.extend(sub_errors)

    elif isinstance(expr, sp.Add):
        # 递归检查加法的每一项
        for i, arg in enumerate(expr.args):
            sub_has_illegal, sub_errors = check_for_illegal_operator_products(
                arg, f"{path}.term[{i}]", depth + 1, max_depth
            )
            errors.extend(sub_errors)

    elif isinstance(expr, Operator):
        # 算符本身是合法的，检查其内部结构
        if hasattr(expr, "left") and hasattr(expr, "right"):
            # NormalOrderedOperator
            sub_has_illegal, sub_errors = check_for_illegal_operator_products(
                expr.left, f"{path}.left", depth + 1, max_depth
            )
            errors.extend(sub_errors)
            sub_has_illegal, sub_errors = check_for_illegal_operator_products(
                expr.right, f"{path}.right", depth + 1, max_depth
            )
            errors.extend(sub_errors)
        elif hasattr(expr, "base"):
            # DerivativeOperator
            sub_has_illegal, sub_errors = check_for_illegal_operator_products(
                expr.base, f"{path}.base", depth + 1, max_depth
            )
            errors.extend(sub_errors)

    return len(errors) > 0, errors


class TestNoIllegalOperatorProducts:
    """测试系统中不会产生非法的 operator * operator 表达式"""

    def test_bc_ghost_system(self):
        """测试 bc ghost 系统不产生非法表达式"""
        b = BasisOperator("b", conformal_weight=2, fermionic=True)
        c = BasisOperator("c", conformal_weight=-1, fermionic=True)
        Fermionic(b, c)
        OPE[b, c] = OPE.make([One])

        T_bc = -2 * NO(b, d(c)) + (-1) * NO(d(b), c)
        result = OPE(T_bc, b)

        # 检查所有极点
        for n in range(1, result.max_pole + 1):
            pole_n = result.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

    def test_kac_moody_system(self):
        """测试 Kac-Moody 代数不产生非法表达式"""
        J = BasisOperator("J", conformal_weight=1)
        Bosonic(J)
        k = sp.Symbol("k")

        OPE[J, J] = OPE.make([k * One, 0, 2 * J])

        # 测试 OPE(NO(J,J), J)
        result = OPE(NO(J, J), J)

        for n in range(1, result.max_pole + 1):
            pole_n = result.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

    def test_virasoro_algebra(self):
        """测试 Virasoro 代数不产生非法表达式"""
        T = BasisOperator("T", conformal_weight=2)
        Bosonic(T)
        c = sp.Symbol("c")

        OPE[T, T] = OPE.make([c / 2 * One, 0, 2 * T, d(T)])

        # 测试 OPE(T, T)
        result = OPE(T, T)
        for n in range(1, result.max_pole + 1):
            pole_n = result.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

        # 测试 OPE(NO(T,T), T)
        result2 = OPE(NO(T, T), T)
        for n in range(1, result2.max_pole + 1):
            pole_n = result2.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

    def test_simplify_with_scalar_one(self):
        """测试 simplify 处理包含 k*One 的表达式"""
        J1 = BasisOperator("J1", conformal_weight=1)
        J2 = BasisOperator("J2", conformal_weight=1)
        Bosonic(J1, J2)

        k = sp.Symbol("k")
        OPE[J1, J2] = OPE.make([k * One, 0, 2 * J1])
        OPE[J2, J2] = OPE.make([0, 0, J2])

        # 创建会触发 _expand_nested_no_left 的表达式
        expr = NO(NO(J1, J2), J2)
        result = simplify(expr)

        has_illegal, errors = check_for_illegal_operator_products(
            result, "simplify_result"
        )
        if has_illegal:
            pytest.fail("\n".join(errors))

    def test_derivative_ope(self):
        """测试导数算符的 OPE 不产生非法表达式"""
        b = BasisOperator("b", conformal_weight=2, fermionic=True)
        c = BasisOperator("c", conformal_weight=-1, fermionic=True)
        Fermionic(b, c)
        OPE[b, c] = OPE.make([One])

        # 测试 OPE(d(b), c)
        result = OPE(d(b), c)
        for n in range(1, result.max_pole + 1):
            pole_n = result.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

        # 测试 OPE(b, d(c))
        result2 = OPE(b, d(c))
        for n in range(1, result2.max_pole + 1):
            pole_n = result2.pole(n)
            has_illegal, errors = check_for_illegal_operator_products(
                pole_n, f"pole({n})"
            )
            if has_illegal:
                pytest.fail("\n".join(errors))

    def test_nested_no_expansion(self):
        """测试嵌套 NO 展开不产生非法表达式"""
        J1 = BasisOperator("J1", conformal_weight=1)
        J2 = BasisOperator("J2", conformal_weight=1)
        J3 = BasisOperator("J3", conformal_weight=1)
        Bosonic(J1, J2, J3)

        k = sp.Symbol("k")
        OPE[J1, J2] = OPE.make([k * One, 0, J1])
        OPE[J1, J3] = OPE.make([0, 0, J1])
        OPE[J2, J3] = OPE.make([0, 0, J2])

        # 创建复杂的嵌套 NO
        expr = NO(NO(J1, J2), J3)
        result = simplify(expr)

        has_illegal, errors = check_for_illegal_operator_products(
            result, "nested_no_result"
        )
        if has_illegal:
            pytest.fail("\n".join(errors))
