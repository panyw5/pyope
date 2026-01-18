"""
表达式化简模块

本模块提供将 OPE 计算结果化简为规范形式的功能。

主要函数：
- simplify(expr): 将表达式化简为排序的 NO product 线性组合
"""

from typing import Any, Dict, List, Tuple
import sympy as sp
from sympy import Add, Mul

from .operators import (
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
)
from .local_operator import (
    is_local_operator,
    extract_scalar_operator,
    collect_operator_terms,
)
from .constants import Zero, One
from .registry import ope_registry


def simplify(expr: Any, expand_derivatives: bool = True) -> Any:
    """
    化简 OPE 表达式为规范形式

    将表达式化简为排序的 NO product 的线性组合。执行以下操作：
    1. 展开嵌套的 NO 乘积
    2. 在 NO 内部按算符顺序排列
    3. 合并同类项
    4. 标准化导数表示
    5. （可选）应用莱布尼茨法则展开正规序的导数

    Args:
        expr: 要化简的表达式（可以是 Operator、OPEData 或符号表达式）
        expand_derivatives: 是否自动展开正规序算符的导数（默认 True）
                           当为 True 时，应用莱布尼茨法则：
                           d^n(NO(A,B)) = Σ_{k=0}^{n} C(n,k) * NO(d^k(A), d^{n-k}(B))

    Returns:
        化简后的表达式

    Examples:
        >>> T = BasisOperator("T")
        >>> J = BasisOperator("J")
        >>> expr = NO(NO(T,J), J)
        >>> simplified = simplify(expr)

        >>> # 自动展开导数（默认行为）
        >>> from pyope import d
        >>> expr = d(NO(T, J))
        >>> simplify(expr)  # 返回 NO(d(T), J) + NO(T, d(J))

        >>> # 禁用导数展开
        >>> simplify(expr, expand_derivatives=False)  # 返回 d(NO(T, J))
    """
    # 处理零
    if expr == 0 or expr == Zero:
        return Zero

    # 处理单位
    if expr == One:
        return One

    # 处理 OPEData
    from .ope_data import OPEData

    if isinstance(expr, OPEData):
        return _simplify_ope_data(expr, expand_derivatives)

    # 处理加法：分别化简每一项
    if isinstance(expr, Add):
        simplified_terms = [simplify(term, expand_derivatives) for term in expr.args]
        return sp.Add(*simplified_terms)

    # 处理标量乘法
    if isinstance(expr, Mul):
        coeff, op = extract_scalar_operator(expr)
        if coeff != 1:
            return coeff * simplify(op, expand_derivatives)

    # 处理算符
    if isinstance(expr, Operator):
        return _simplify_operator(expr, expand_derivatives)

    # 其他情况直接返回
    return expr


def _simplify_ope_data(
    ope_data: "OPEData", expand_derivatives: bool = True
) -> "OPEData":
    """
    化简 OPEData 对象

    Args:
        ope_data: OPEData 实例
        expand_derivatives: 是否展开导数

    Returns:
        化简后的 OPEData
    """
    from .ope_data import OPEData

    new_poles = {}
    for q, coeff in ope_data.poles.items():
        simplified_coeff = simplify(coeff, expand_derivatives)
        if simplified_coeff != 0:
            new_poles[q] = simplified_coeff

    return OPEData(new_poles)


def _simplify_operator(op: Operator, expand_derivatives: bool = True) -> Any:
    """
    化简单个算符

    Args:
        op: 算符实例
        expand_derivatives: 是否展开导数

    Returns:
        化简后的表达式
    """
    # 处理 DerivativeOperator：应用莱布尼茨法则展开 d(NO(...))
    if isinstance(op, DerivativeOperator) and expand_derivatives:
        base = op.base
        n = op.order

        # 检查 base 是否为 NormalOrderedOperator
        if isinstance(base, NormalOrderedOperator):
            # 应用莱布尼茨法则: d^n(NO(A,B)) = Σ_{k=0}^{n} C(n,k) * NO(d^k(A), d^{n-k}(B))
            from .api import NO
            from .operators import d as derivative_operator
            from .cache import cached_binomial

            # 先化简左右算符
            left = simplify(base.left, expand_derivatives)
            right = simplify(base.right, expand_derivatives)

            # 生成莱布尼茨展开的各项
            terms = []
            for k in range(n + 1):
                coeff = cached_binomial(n, k)

                # 计算 d^k(left) 和 d^(n-k)(right)
                left_deriv = derivative_operator(left, k) if k > 0 else left
                right_deriv = derivative_operator(right, n - k) if n - k > 0 else right

                # 构造 NO 项
                terms.append(coeff * NO(left_deriv, right_deriv))

            # 递归化简展开后的表达式
            result = sp.Add(*terms) if len(terms) > 1 else terms[0]
            return simplify(result, expand_derivatives)

        # 对于非 NO 的 base，保持 DerivativeOperator 结构
        # 但可以递归化简 base（可选）
        simplified_base = simplify(base, expand_derivatives)
        if simplified_base != base:
            return DerivativeOperator(simplified_base, n)
        return op

    # BasisOperator 和 DerivativeOperator（未展开）已经是最简形式
    if isinstance(op, (BasisOperator, DerivativeOperator)):
        return op

    # 处理 NormalOrderedOperator
    if isinstance(op, NormalOrderedOperator):
        return _simplify_normal_ordered(op, expand_derivatives)

    return op


def _simplify_normal_ordered(
    no_op: NormalOrderedOperator, expand_derivatives: bool = True
) -> Any:
    """
    化简正规序算符

    处理：
    1. 递归化简左右算符
    2. 展开嵌套的 NO
    3. 排序内部算符

    Args:
        no_op: NormalOrderedOperator 实例
        expand_derivatives: 是否展开导数

    Returns:
        化简后的表达式
    """
    from .api import NO

    # 递归化简左右算符
    left = simplify(no_op.left, expand_derivatives)
    right = simplify(no_op.right, expand_derivatives)

    # 如果左侧或右侧是加法，分配
    if isinstance(left, Add):
        return sp.Add(*[NO(term, right) for term in left.args])
    if isinstance(right, Add):
        return sp.Add(*[NO(left, term) for term in right.args])

    # 处理标量乘法
    if isinstance(left, Mul):
        coeff, op = extract_scalar_operator(left)
        if coeff != 1:
            return coeff * NO(op, right)
    if isinstance(right, Mul):
        coeff, op = extract_scalar_operator(right)
        if coeff != 1:
            return coeff * NO(left, op)

    # 处理嵌套的 NO
    # 如果左侧是 NO，展开为 NO(NO(A,B), C)
    if isinstance(left, NormalOrderedOperator):
        # 这需要使用 OPE 来展开
        # 暂时保持原样，因为完整展开需要 OPE 信息
        pass

    # 如果右侧是 NO，展开为 NO(A, NO(B,C))
    if isinstance(right, NormalOrderedOperator):
        # 同样暂时保持原样
        pass

    # 检查算符顺序
    # 只有当左右都是 BasisOperator 或 DerivativeOperator 时才进行交换
    # 嵌套的 NO 已经在上面处理过了（虽然目前是 pass）
    if isinstance(left, (BasisOperator, DerivativeOperator)) and isinstance(
        right, (BasisOperator, DerivativeOperator)
    ):
        order = ope_registry.compare_operators(left, right)
        if order < 0:
            # 需要交换顺序: NO(B, A) -> 使用 eq (2.3.16) 计算
            # 公式 (2.3.16): [BA]_q = (-1)^{|A||B|} \sum_{l >= q} \frac{(-1)^l}{(l-q)!} \partial^{(l-q)} [AB]_l
            # 对于 q=0 (正规序): [BA]_0 = (-1)^{|A||B|} \sum_{l >= 0} \frac{(-1)^l}{l!} \partial^l [AB]_l

            # 1. 计算符号因子 (-1)^{|A||B|}
            from .local_operator import get_operator_parity

            parity_A = get_operator_parity(right)  # right 是 A
            parity_B = get_operator_parity(left)  # left 是 B
            swap_sign = (-1) ** (parity_A * parity_B)

            # 2. 计算 OPE(A, B)（注意顺序：right 是 A，left 是 B）
            from .api import OPE
            from .operators import d as derivative_operator

            ope_AB = OPE(right, left)  # OPE(A, B)

            # 3. 应用公式 (2.3.16) 对 q=0
            # [BA]_0 = (-1)^{|A||B|} \sum_{l >= 0} \frac{(-1)^l}{l!} \partial^l [AB]_l

            # 检查是否有非零极点
            max_pole = ope_AB.max_pole
            has_nonzero_poles = any(ope_AB.pole(l) != 0 for l in range(max_pole + 1))

            if not has_nonzero_poles:
                # OPE 未定义或全为 0，直接返回正则顺序的 NO
                # 根据 eq (2.3.16)，[BA]_0 = swap_sign * [AB]_0
                # 当 [AB]_l = 0 对所有 l 时，[BA]_0 = swap_sign * NO(A,B)
                return swap_sign * NO(right, left)

            # 有非零极点，应用完整公式
            from .constants import Zero as ZeroConst

            result = ZeroConst

            # l=0 项: [AB]_0
            pole_0 = ope_AB.pole(0)
            if pole_0 != 0 and pole_0 != ZeroConst:
                result = swap_sign * pole_0

            # l >= 1 项: \frac{(-1)^l}{l!} \partial^l [AB]_l
            for l in range(1, max_pole + 1):
                pole_l = ope_AB.pole(l)
                if pole_l != 0 and pole_l != ZeroConst:
                    # 计算 \partial^l [AB]_l
                    deriv_pole = derivative_operator(pole_l, l)
                    # 系数: swap_sign * (-1)^l / l!
                    coeff = swap_sign * ((-1) ** l) / sp.factorial(l)
                    result = result + coeff * deriv_pole

            # 递归化简结果
            return simplify(result, expand_derivatives)

    # 创建简化的 NO
    return NO(left, right)


def collect_normal_ordered_terms(expr: Any) -> Dict[Tuple, Any]:
    """
    收集并合并正规序项

    将表达式中的所有 NO(...) 项按结构分组并合并系数

    Args:
        expr: 表达式

    Returns:
        字典，键为 NO 的规范形式，值为系数

    Examples:
        >>> expr = 2*NO(T,J) + 3*NO(T,J) + NO(J,J)
        >>> collect_normal_ordered_terms(expr)
        {('T', 'J'): 5, ('J', 'J'): 1}
    """
    terms = {}

    if isinstance(expr, Add):
        for term in expr.args:
            sub_terms = collect_normal_ordered_terms(term)
            for key, coeff in sub_terms.items():
                terms[key] = terms.get(key, 0) + coeff
    elif isinstance(expr, Mul):
        coeff, op = extract_scalar_operator(expr)
        if isinstance(op, NormalOrderedOperator):
            key = _make_no_key(op)
            terms[key] = terms.get(key, 0) + coeff
        else:
            # 非 NO 项
            key = ("other", str(op))
            terms[key] = terms.get(key, 0) + coeff
    elif isinstance(expr, NormalOrderedOperator):
        key = _make_no_key(expr)
        terms[key] = terms.get(key, 0) + 1
    else:
        # 其他项
        key = ("other", str(expr))
        terms[key] = terms.get(key, 0) + 1

    return terms


def _make_no_key(no_op: NormalOrderedOperator) -> Tuple:
    """
    为 NO 算符创建唯一键

    Args:
        no_op: NormalOrderedOperator 实例

    Returns:
        元组键
    """
    left_key = _operator_to_key(no_op.left)
    right_key = _operator_to_key(no_op.right)
    return ("NO", left_key, right_key)


def _operator_to_key(op: Any) -> Tuple:
    """
    将算符转换为可哈希的键

    Args:
        op: 算符

    Returns:
        元组键
    """
    if isinstance(op, BasisOperator):
        return ("Basis", op.name, op.is_bosonic)
    elif isinstance(op, DerivativeOperator):
        return ("Deriv", _operator_to_key(op.base), op.order)
    elif isinstance(op, NormalOrderedOperator):
        return _make_no_key(op)
    else:
        return ("Other", str(op))
