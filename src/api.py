"""
Core API functions for PyOPE.

This module provides the main user-facing functions:
- OPE: Compute operator product expansions
- NO: Normal ordered products
- bracket: Commutator/anticommutator
"""

from __future__ import annotations
from typing import Union, Optional
import sympy as sp
from .operators import Operator, NormalOrderedOperator, DerivativeOperator
from .local_operator import LocalOperator, OperatorSum, OperatorProduct, extract_scalar_operator
from .ope_data import OPEData
from .registry import ope_registry


def operator_derivative(expr: sp.Expr, order: int = 1) -> sp.Expr:
    """
    对算符表达式求导。

    在共形场论中，算符的导数是对坐标的导数。
    例如：∂T 表示 T 对坐标的导数。

    Args:
        expr: 算符表达式
        order: 导数阶数

    Returns:
        导数后的表达式
    """
    if order == 0:
        return expr

    if order < 0:
        raise ValueError("Derivative order must be non-negative")

    if isinstance(expr, DerivativeOperator):
        # ∂^n A 的导数是 ∂^(n+order) A
        return DerivativeOperator(expr.base, expr.order + order)
    elif isinstance(expr, Operator):
        # A 的导数是 ∂^order A
        return DerivativeOperator(expr, order)
    elif isinstance(expr, sp.Add):
        # 线性性：(A + B)' = A' + B'
        return sum(operator_derivative(term, order) for term in expr.args)
    elif isinstance(expr, sp.Mul):
        # Leibniz 规则：(c*A)' = c*A' (假设 c 是常数)
        scalar_parts = []
        operator_part = None

        for arg in expr.args:
            if isinstance(arg, Operator):
                if operator_part is not None:
                    # 多个算符的乘积需要完整的 Leibniz 规则
                    # 这里暂不实现，因为通常极点是单个算符
                    raise NotImplementedError(
                        "Derivative of product of multiple operators not implemented"
                    )
                operator_part = arg
            else:
                scalar_parts.append(arg)

        if operator_part is None:
            # 纯标量，导数为 0
            return sp.S.Zero

        scalar = sp.Mul(*scalar_parts) if scalar_parts else sp.S.One
        return scalar * operator_derivative(operator_part, order)
    else:
        # 常数或其他，导数为 0
        return sp.S.Zero


def swap_sign(A: Operator, B: Operator) -> int:
    """
    计算交换两个算符的符号因子。

    基于 OPEdefs.m 中的 SwapSign 函数。

    规则:
    - 如果至少一个算符是玻色子,返回 +1
    - 如果两个算符都是费米子,返回 -1

    Args:
        A: 第一个算符
        B: 第二个算符

    Returns:
        +1 或 -1
    """
    # 获取算符的玻色/费米属性
    # 对于导数算符,使用基算符的属性
    if isinstance(A, DerivativeOperator):
        A_bosonic = A.base.is_bosonic
    else:
        A_bosonic = A.is_bosonic

    if isinstance(B, DerivativeOperator):
        B_bosonic = B.base.is_bosonic
    else:
        B_bosonic = B.is_bosonic

    # 如果至少一个是玻色子,返回 +1
    if A_bosonic or B_bosonic:
        return 1
    else:
        # 两个都是费米子,返回 -1
        return -1


def ope_commute_help(B: Operator, A: Operator) -> OPEData:
    """
    计算交换后的 OPE: OPE[B, A] 从 OPE[A, B]。

    基于 OPEdefs.m 第 959-972 行的公式:
    OPE[B,A](q) = SwapSign[A,B] * Sum[(-1)^l / (l-q)! * D^(l-q)[pole_l(OPE[A,B])],
                                      {l, q, maxpole}]

    Args:
        B: 第一个算符 (交换后的顺序)
        A: 第二个算符

    Returns:
        OPEData 包含 OPE[B, A] 的结果
    """
    # 获取 OPE[A, B]
    AB = OPE(A, B)
    max_pole = AB.max_pole

    if max_pole == 0:
        # 如果 OPE[A, B] 为零,则 OPE[B, A] 也为零
        return OPEData()

    # 计算符号因子
    sign = swap_sign(A, B)

    # 计算新的极点
    new_poles = {}
    for q in range(max_pole, 0, -1):
        term = sp.S.Zero

        for l in range(q, max_pole + 1):
            pole_l = AB.pole(l)
            if pole_l != 0:
                # 计算 D^(l-q)[pole_l]
                deriv_order = l - q
                if deriv_order == 0:
                    deriv_pole = pole_l
                else:
                    deriv_pole = operator_derivative(pole_l, deriv_order)

                # 系数: (-1)^l / (l-q)!
                coeff = (-1)**l / sp.factorial(deriv_order)
                term += coeff * deriv_pole

        # 简化结果
        try:
            term = sp.expand(term)
        except:
            pass

        if term != 0:
            new_poles[q] = sign * term

    return OPEData(new_poles)


def ope_composite_help_rq(A: Operator, BC: 'NormalOrderedOperator') -> OPEData:
    """
    计算右复合算符的 OPE: OPE[A, NO[B,C]]。

    基于 OPEdefs.m 第 985-1016 行的公式。

    公式:
    OPE[A, NO[B,C]] = SwapSign[A,B] * Sum[NO[B, pole_j(OPE[A,C])], {j, maxpole_AC}]
                    + Sum[NO[pole_k(OPE[A,B]), C], {k, maxpole_AB}]
                    + Sum[Sum[pole_j(OPE[A,C]) * pole_k(OPE[A,B]), nested terms]]

    Args:
        A: 第一个算符
        BC: 正规序算符 NO[B,C]

    Returns:
        OPEData 包含计算结果
    """
    B = BC.left
    C = BC.right

    # 计算 OPE[A, B] 和 OPE[A, C]
    AB = OPE(A, B)
    AC = OPE(A, C)

    max_pole_AB = AB.max_pole
    max_pole_AC = AC.max_pole

    # 如果两个 OPE 都为零,结果为零
    if max_pole_AB == 0 and max_pole_AC == 0:
        return OPEData()

    # 计算符号因子
    sign = swap_sign(A, B)

    # 收集所有极点
    new_poles = {}

    # 第一项: SwapSign[A,B] * NO[B, pole_j(OPE[A,C])]
    for j in range(max_pole_AC, 0, -1):
        pole_j = AC.pole(j)
        if pole_j != 0:
            # NO[B, pole_j]
            term = sign * NO(B, pole_j)
            if j in new_poles:
                new_poles[j] += term
            else:
                new_poles[j] = term

    # 第二项: NO[pole_k(OPE[A,B]), C]
    for k in range(max_pole_AB, 0, -1):
        pole_k = AB.pole(k)
        if pole_k != 0:
            # NO[pole_k, C]
            term = NO(pole_k, C)
            if k in new_poles:
                new_poles[k] += term
            else:
                new_poles[k] = term

    # 第三项: 嵌套 OPE 项
    # 这部分涉及 OPE[pole_j(OPE[A,C]), pole_k(OPE[A,B])]
    # 对于每对极点,计算它们的 OPE 并加到相应的极点阶数
    for j in range(max_pole_AC, 0, -1):
        pole_j_AC = AC.pole(j)
        if pole_j_AC == 0:
            continue

        for k in range(max_pole_AB, 0, -1):
            pole_k_AB = AB.pole(k)
            if pole_k_AB == 0:
                continue

            # 计算 OPE[pole_j_AC, pole_k_AB]
            # 注意: 这里需要考虑符号
            nested_ope = sign * OPE(pole_j_AC, pole_k_AB)

            # 将嵌套 OPE 的极点加到结果中
            for m in range(nested_ope.max_pole, 0, -1):
                pole_m = nested_ope.pole(m)
                if pole_m != 0:
                    total_order = j + k + m
                    if total_order in new_poles:
                        new_poles[total_order] += pole_m
                    else:
                        new_poles[total_order] = pole_m

    # 简化所有极点
    simplified_poles = {}
    for order, pole in new_poles.items():
        try:
            simplified = sp.expand(pole)
        except:
            simplified = pole
        if simplified != 0:
            simplified_poles[order] = simplified

    return OPEData(simplified_poles)


def ope_composite_help_lq(AB: 'NormalOrderedOperator', C: Operator) -> OPEData:
    """
    计算左复合算符的 OPE: OPE[NO[A,B], C]。

    基于 OPEdefs.m 第 1028-1084 行的公式。

    这个公式更复杂,涉及:
    - OPE[A, OPE[B, C]]
    - 导数项
    - 嵌套求和

    Args:
        AB: 正规序算符 NO[A,B]
        C: 第二个算符

    Returns:
        OPEData 包含计算结果
    """
    A = AB.left
    B = AB.right

    # 计算 OPE[A, C] 和 OPE[B, C]
    AC = OPE(A, C)
    BC = OPE(B, C)

    max_pole_AC = AC.max_pole
    max_pole_BC = BC.max_pole

    # 如果两个 OPE 都为零,结果为零
    if max_pole_AC == 0 and max_pole_BC == 0:
        return OPEData()

    # 收集所有极点
    new_poles = {}

    # 第一项: NO[A, pole_j(OPE[B,C])]
    for j in range(max_pole_BC, 0, -1):
        pole_j = BC.pole(j)
        if pole_j != 0:
            term = NO(A, pole_j)
            if j in new_poles:
                new_poles[j] += term
            else:
                new_poles[j] = term

    # 第二项: SwapSign[A,B] * NO[B, pole_k(OPE[A,C])]
    sign = swap_sign(A, B)
    for k in range(max_pole_AC, 0, -1):
        pole_k = AC.pole(k)
        if pole_k != 0:
            term = sign * NO(B, pole_k)
            if k in new_poles:
                new_poles[k] += term
            else:
                new_poles[k] = term

    # 第三项: 嵌套 OPE 和导数项
    # 对于每个 pole_j(OPE[B,C]),计算 OPE[A, pole_j]
    for j in range(max_pole_BC, 0, -1):
        pole_j_BC = BC.pole(j)
        if pole_j_BC == 0:
            continue

        # 计算 OPE[A, pole_j_BC]
        A_pole_j = OPE(A, pole_j_BC)

        # 将结果的极点加到总结果中
        for m in range(A_pole_j.max_pole, 0, -1):
            pole_m = A_pole_j.pole(m)
            if pole_m != 0:
                total_order = j + m
                if total_order in new_poles:
                    new_poles[total_order] += pole_m
                else:
                    new_poles[total_order] = pole_m

    # 第四项: 对于每个 pole_k(OPE[A,C]),计算 OPE[B, pole_k]
    for k in range(max_pole_AC, 0, -1):
        pole_k_AC = AC.pole(k)
        if pole_k_AC == 0:
            continue

        # 计算 OPE[B, pole_k_AC]
        B_pole_k = OPE(B, pole_k_AC)

        # 将结果的极点加到总结果中,带符号
        for m in range(B_pole_k.max_pole, 0, -1):
            pole_m = B_pole_k.pole(m)
            if pole_m != 0:
                total_order = k + m
                term = sign * pole_m
                if total_order in new_poles:
                    new_poles[total_order] += term
                else:
                    new_poles[total_order] = term

    # 简化所有极点
    simplified_poles = {}
    for order, pole in new_poles.items():
        try:
            simplified = sp.expand(pole)
        except:
            simplified = pole
        if simplified != 0:
            simplified_poles[order] = simplified

    return OPEData(simplified_poles)


def ope_derivative_help_l(A_deriv: DerivativeOperator, B: Operator) -> OPEData:
    """
    计算左导数算符的 OPE: OPE[∂^i A, B]。

    基于 OPEdefs.m 第 910-920 行的公式：
    OPE[∂^i A, B] = (-1)^i * Sum[Pochhammer[j,i] * pole_j(OPE[A,B]), {j, maxpole, 1, -1}]

    其中 Pochhammer[j,i] = j*(j+1)*...*(j+i-1) 是上升阶乘 (rising factorial)。

    Args:
        A_deriv: 导数算符 ∂^i A
        B: 第二个算符

    Returns:
        OPEData 包含计算结果
    """
    base = A_deriv.base
    order = A_deriv.order

    # 获取 OPE[A, B]
    AB = OPE(base, B)
    max_pole = AB.max_pole

    if max_pole == 0:
        # 如果 OPE[A, B] 为零，则导数也为零
        return OPEData()

    # 计算新的极点
    new_poles = {}
    for j in range(max_pole, 0, -1):
        pole_j = AB.pole(j)
        if pole_j != 0:
            # Pochhammer[j, order] = j * (j+1) * ... * (j+order-1)
            # 在 SymPy 中是 rf(j, order) - rising factorial
            coeff = sp.rf(j, order)  # rising factorial
            new_pole_order = j + order
            term = (-1)**order * coeff * pole_j

            # 简化表达式
            try:
                term = sp.expand(term)
            except:
                pass

            if term != 0:
                new_poles[new_pole_order] = term

    return OPEData(new_poles)


def ope_derivative_help_r(A: Operator, B_deriv: DerivativeOperator) -> OPEData:
    """
    计算右导数算符的 OPE: OPE[A, ∂^i B]。

    基于 OPEdefs.m 第 937-948 行的公式：
    OPE[A, ∂^i B] = Sum[Binomial[i,k] * Pochhammer[j-k,k] * ∂^(i-k)[pole_j(OPE[A,B])],
                        {j, maxpole+i, 1, -1}, {k, 0, min(i, j-1)}]

    Args:
        A: 第一个算符
        B_deriv: 导数算符 ∂^i B

    Returns:
        OPEData 包含计算结果
    """
    base = B_deriv.base
    order = B_deriv.order  # i

    # 获取 OPE[A, B]
    AB = OPE(A, base)
    max_pole = AB.max_pole

    if max_pole == 0:
        return OPEData()

    # 预计算所有需要的导数
    # derivatives[deriv_order][j-1] 存储 ∂^deriv_order[pole_j]
    derivatives = {}
    derivatives[0] = {}
    for j in range(1, max_pole + 1):
        derivatives[0][j] = AB.pole(j)

    for deriv_order in range(1, order + 1):
        derivatives[deriv_order] = {}
        for j in range(1, max_pole + 1):
            pole_val = derivatives[0][j]
            if pole_val != 0:
                derivatives[deriv_order][j] = operator_derivative(pole_val, deriv_order)
            else:
                derivatives[deriv_order][j] = sp.S.Zero

    # 计算新的极点
    new_poles = {}
    for j in range(max_pole + order, 0, -1):
        result = sp.S.Zero

        # k 的范围：max(0, j - max_pole) <= k <= min(order, j-1)
        k_min = max(0, j - max_pole)
        k_max = min(order, j - 1)

        for k in range(k_min, k_max + 1):
            j_orig = j - k  # 原始极点的阶数
            if j_orig >= 1 and j_orig <= max_pole:
                deriv_order = order - k
                pole_val = derivatives[deriv_order].get(j_orig, sp.S.Zero)

                if pole_val != 0:
                    # Binomial[order, k] * Pochhammer[j-k, k]
                    coeff = sp.binomial(order, k) * sp.rf(j - k, k)
                    result += coeff * pole_val

        # 简化结果
        try:
            result = sp.expand(result)
        except:
            pass

        if result != 0:
            new_poles[j] = result

    return OPEData(new_poles)


class OPEFunction:
    """
    Main OPE computation function.

    This class implements the OPE computation with automatic linearity expansion.
    """

    def __call__(
        self, A: Union[Operator, sp.Expr], B: Union[Operator, sp.Expr]
    ) -> OPEData:
        """
        Compute OPE of A(z) with B(w).

        Args:
            A: First operator or expression
            B: Second operator or expression

        Returns:
            OPEData containing the poles of the OPE
        """
        # Handle OperatorSum: OPE(A, B + C) = OPE(A, B) + OPE(A, C)
        if isinstance(B, OperatorSum):
            return sum(self(A, term) for term in B.terms)

        if isinstance(A, OperatorSum):
            return sum(self(term, B) for term in A.terms)

        # Handle OperatorProduct: OPE(c*A, B) = c * OPE(A, B)
        if isinstance(B, OperatorProduct):
            return B.coeff * self(A, B.operator)

        if isinstance(A, OperatorProduct):
            return A.coeff * self(A.operator, B)

        # Handle legacy sp.Add (for backward compatibility)
        if isinstance(B, sp.Add):
            return sum(self(A, term) for term in B.args)

        if isinstance(A, sp.Add):
            return sum(self(term, B) for term in A.args)

        # Handle legacy sp.Mul (for backward compatibility)
        if isinstance(B, sp.Mul):
            scalar, op = self._extract_scalar(B)
            if scalar is not None:
                return scalar * self(A, op)

        if isinstance(A, sp.Mul):
            scalar, op = self._extract_scalar(A)
            if scalar is not None:
                return scalar * self(op, B)

        # Handle zero
        if A == 0 or B == 0:
            return OPEData()

        # Now A and B should be operators
        if not isinstance(A, Operator):
            raise TypeError(f"Expected Operator, got {type(A)}")
        if not isinstance(B, Operator):
            raise TypeError(f"Expected Operator, got {type(B)}")

        # Dispatch to appropriate computation method
        return self._compute(A, B)

    def _extract_scalar(
        self, expr: sp.Expr
    ) -> tuple[Optional[sp.Expr], Optional[Operator]]:
        """
        Extract scalar coefficient from expression.

        Args:
            expr: Expression to analyze

        Returns:
            Tuple of (scalar, operator) or (None, None)
        """
        if not isinstance(expr, sp.Mul):
            return None, None

        # Separate scalar and operator parts
        scalar_parts = []
        operator_part = None

        for arg in expr.args:
            if isinstance(arg, Operator):
                if operator_part is not None:
                    # Multiple operators, can't extract scalar
                    return None, None
                operator_part = arg
            else:
                scalar_parts.append(arg)

        if operator_part is None:
            return None, None

        scalar = sp.Mul(*scalar_parts) if scalar_parts else sp.S.One
        return scalar, operator_part

    def _compute(self, A: Operator, B: Operator) -> OPEData:
        """
        Compute OPE between two operators.

        Args:
            A: First operator
            B: Second operator

        Returns:
            OPEData result
        """
        # Check if A is a normal ordered operator
        if isinstance(A, NormalOrderedOperator):
            return ope_composite_help_lq(A, B)

        # Check if B is a normal ordered operator
        if isinstance(B, NormalOrderedOperator):
            return ope_composite_help_rq(A, B)

        # Check if A is a derivative operator
        if isinstance(A, DerivativeOperator):
            return ope_derivative_help_l(A, B)

        # Check if B is a derivative operator
        if isinstance(B, DerivativeOperator):
            return ope_derivative_help_r(A, B)

        # Check for defined OPE
        result = ope_registry.lookup(A, B)
        if result is not None:
            return result

        # Try commutation: if OPE[B, A] is defined, compute OPE[A, B] from it
        result_ba = ope_registry.lookup(B, A)
        if result_ba is not None:
            return ope_commute_help(A, B)

        # Default: OPE is zero (regular)
        return OPEData()

    def __setitem__(
        self, key: tuple[Operator, Operator], value: Union[OPEData, list]
    ) -> None:
        """
        Define an OPE: OPE[A, B] = ...

        Args:
            key: Tuple of (A, B)
            value: OPEData or list of poles
        """
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("Key must be a tuple of two operators")

        A, B = key

        if isinstance(value, list):
            value = OPEData.from_list(value)

        if not isinstance(value, OPEData):
            raise TypeError(f"Value must be OPEData or list, got {type(value)}")

        ope_registry.define(A, B, value)

    @staticmethod
    def make(poles: Union[list, OPEData]) -> OPEData:
        """
        Create OPEData from poles list (Mathematica-style).

        Args:
            poles: List of poles [pole_n, ..., pole_1] or OPEData

        Returns:
            OPEData object
        """
        if isinstance(poles, OPEData):
            return poles
        return OPEData.from_list(poles)


class NOFunction:
    """
    Normal ordered product function.

    Computes NO[A, B] = (AB)(z) using point-splitting convention.
    """

    def __call__(
        self, A: Union[Operator, sp.Expr], B: Union[Operator, sp.Expr]
    ) -> Union[NormalOrderedOperator, sp.Expr]:
        """
        Compute normal ordered product of A and B.

        Args:
            A: First operator or expression
            B: Second operator or expression

        Returns:
            NormalOrderedOperator or simplified expression
        """
        # Handle OperatorSum
        if isinstance(B, OperatorSum):
            return sum(self(A, term) for term in B.terms)

        if isinstance(A, OperatorSum):
            return sum(self(term, B) for term in A.terms)

        # Handle OperatorProduct
        if isinstance(B, OperatorProduct):
            return B.coeff * self(A, B.operator)

        if isinstance(A, OperatorProduct):
            return A.coeff * self(A.operator, B)

        # Handle legacy sp.Add (for backward compatibility)
        if isinstance(B, sp.Add):
            return sum(self(A, term) for term in B.args)

        if isinstance(A, sp.Add):
            return sum(self(term, B) for term in A.args)

        # Handle legacy sp.Mul (for backward compatibility)
        if isinstance(B, sp.Mul):
            scalar, op = self._extract_scalar(B)
            if scalar is not None:
                return scalar * self(A, op)

        if isinstance(A, sp.Mul):
            scalar, op = self._extract_scalar(A)
            if scalar is not None:
                return scalar * self(op, B)

        # Handle zero
        if A == 0 or B == 0:
            return sp.S.Zero

        # Now A and B should be operators
        if not isinstance(A, Operator):
            raise TypeError(f"Expected Operator, got {type(A)}")
        if not isinstance(B, Operator):
            raise TypeError(f"Expected Operator, got {type(B)}")

        # Create normal ordered operator
        return NormalOrderedOperator(A, B)

    def _extract_scalar(
        self, expr: sp.Expr
    ) -> tuple[Optional[sp.Expr], Optional[Operator]]:
        """Extract scalar coefficient from expression."""
        if not isinstance(expr, sp.Mul):
            return None, None

        scalar_parts = []
        operator_part = None

        for arg in expr.args:
            if isinstance(arg, Operator):
                if operator_part is not None:
                    return None, None
                operator_part = arg
            else:
                scalar_parts.append(arg)

        if operator_part is None:
            return None, None

        scalar = sp.Mul(*scalar_parts) if scalar_parts else sp.S.One
        return scalar, operator_part


def bracket(
    A: Operator, B: Operator, anticommutator: bool = False
) -> sp.Expr:
    """
    Compute commutator or anticommutator.

    [A, B] = AB - BA (commutator)
    {A, B} = AB + BA (anticommutator)

    Args:
        A: First operator
        B: Second operator
        anticommutator: If True, compute {A,B}, otherwise [A,B]

    Returns:
        Expression for the (anti)commutator
    """
    # This is a placeholder - full implementation requires OPE computation
    # For now, just return the symbolic expression
    AB = NO(A, B)
    BA = NO(B, A)

    if anticommutator:
        return AB + BA
    else:
        return AB - BA


# Create global instances
OPE = OPEFunction()
NO = NOFunction()


__all__ = ["OPE", "NO", "bracket"]
