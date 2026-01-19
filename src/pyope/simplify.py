"""
表达式化简模块

本模块提供将 OPE 计算结果化简为规范形式的功能。

主要函数：
- simplify(expr): 将表达式化简为排序的 NO product 线性组合
"""

from typing import Any, Dict, List, Tuple
import sympy as sp
from sympy import Add, Mul, Pow

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


def normalize_identity(expr: Any) -> Any:
    """
    正规化表达式中的 One 常数

    规则：
    - b*One, One*b -> b (如果 b 是算符)
    - NO(X, One), NO(One, X) -> X (如果 X 是算符)
    - c*One -> c*One (如果 c 是纯标量，保留 One)
    - One - 1 -> 0
    - One + (-1) -> 0

    Args:
        expr: 要正规化的表达式

    Returns:
        正规化后的表达式
    """
    if expr == 0 or expr == Zero:
        return Zero

    if expr == One:
        return One

    # 处理 One - 1, One + (-1) 等
    if isinstance(expr, Add):
        # 检查是否有 One 和 -1 或 1 和 -One
        has_one = One in expr.args or 1 in expr.args
        has_neg_one = -One in expr.args or -1 in expr.args

        if has_one and has_neg_one:
            # 移除 One 和 -1 (或 1 和 -One)
            new_args = []
            removed_one = False
            removed_neg_one = False

            for arg in expr.args:
                if not removed_one and (arg == One or arg == 1):
                    removed_one = True
                    continue
                if not removed_neg_one and (arg == -One or arg == -1):
                    removed_neg_one = True
                    continue
                new_args.append(normalize_identity(arg))

            if not new_args:
                return Zero
            if len(new_args) == 1:
                return new_args[0]
            return sp.Add(*new_args)

        # 递归处理每一项
        return sp.Add(*[normalize_identity(arg) for arg in expr.args])

    # 处理 Mul: b*One -> b
    if isinstance(expr, Mul):
        # 先递归处理每个因子（特别是 Add 类型的因子）
        normalized_args = [normalize_identity(arg) for arg in expr.args]

        # 检查是否有 Zero
        if Zero in normalized_args or 0 in normalized_args:
            return Zero

        # 分离标量和算符
        scalar_part = []
        operator_part = []
        has_one = False

        for arg in normalized_args:
            if arg == One:
                has_one = True
            elif isinstance(arg, Operator):
                operator_part.append(arg)
            else:
                scalar_part.append(arg)

        # 如果有算符，移除 One
        if operator_part and has_one:
            if scalar_part:
                return sp.Mul(*scalar_part, *operator_part)
            else:
                if len(operator_part) == 1:
                    return operator_part[0]
                return sp.Mul(*operator_part)

        # 如果没有算符，保留 One（纯标量情况）
        if normalized_args != list(expr.args):
            return sp.Mul(*normalized_args)
        return expr

    # 处理 NormalOrderedOperator: NO(X, One) -> X
    if isinstance(expr, NormalOrderedOperator):
        left = normalize_identity(expr.left)
        right = normalize_identity(expr.right)

        if left == One:
            return right
        if right == One:
            return left

        if left != expr.left or right != expr.right:
            from .api import NO

            return NO(left, right)

    # 处理 Pow: Operator**n -> 嵌套的 NO
    if isinstance(expr, Pow):
        base = expr.base
        exp = expr.exp

        # 只处理算符的整数次幂
        if isinstance(base, Operator) and exp.is_integer and exp >= 2:
            from .api import NO

            result = base
            for _ in range(int(exp) - 1):
                result = NO(result, base)
            return result

    return expr


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
        # 过滤掉 Zero 项
        non_zero_terms = [t for t in simplified_terms if t != Zero and t != 0]
        if not non_zero_terms:
            return Zero
        if len(non_zero_terms) == 1:
            return non_zero_terms[0]
        return sp.Add(*non_zero_terms)
    # 处理标量乘法
    if isinstance(expr, Mul):
        coeff, op = extract_scalar_operator(expr)

        # 关键：如果 op 本身是多个算符的普通乘积 Mul(op1, op2, ...)，
        # 则 extract_scalar_operator 无法继续降低复杂度，会导致 simplify 无限递归。
        # 在 VOA 语义下，这类乘积应当被解释为嵌套的正规序算符。
        if (
            isinstance(op, Mul)
            and all(isinstance(a, Operator) for a in op.args)
            and len(op.args) >= 2
        ):
            from .api import NO

            nested = op.args[0]
            for factor in op.args[1:]:
                nested = NO(nested, factor)
            op = nested

        simplified_op = simplify(op, expand_derivatives)
        # 关键：把 (-1)*Zero, k*Zero 等形式规约到 Zero
        if simplified_op == Zero or simplified_op == 0:
            return Zero
        if coeff != 1:
            return coeff * simplified_op

    # 处理算符
    if isinstance(expr, Operator):
        result = _simplify_operator(expr, expand_derivatives)
        return normalize_identity(result)

    # 其他情况：正规化后返回
    return normalize_identity(expr)


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
        if simplified_coeff != 0 and simplified_coeff != Zero:
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

    # 如果左侧或右侧是加法，分配（保持标量系数，避免 NO 接收到 Mul）
    if isinstance(left, Add):
        terms = []
        for term in left.args:
            if isinstance(term, Mul):
                coeff, op = extract_scalar_operator(term)
                if coeff != 1:
                    terms.append(coeff * NO(op, right))
                    continue
            terms.append(NO(term, right))
        return sp.Add(*terms)

    if isinstance(right, Add):
        terms = []
        for term in right.args:
            if isinstance(term, Mul):
                coeff, op = extract_scalar_operator(term)
                # coeff!=1: 直接外提
                if coeff != 1:
                    terms.append(coeff * NO(left, op))
                    continue
                # coeff==1 但 op 仍可能是 Mul(Operator, Operator)
                if isinstance(op, Operator):
                    terms.append(NO(left, op))
                    continue
                # op 不是单个 Operator（例如 Mul(op1, op2, ...)），
                # 把它解释为多个算符的乘积，并左结合构造嵌套 NO。
                if isinstance(op, Mul):
                    nested = left
                    for factor in op.args:
                        if not isinstance(factor, Operator):
                            # 遇到纯标量，外提
                            nested = sp.Mul(factor, nested)
                            continue
                        nested = NO(nested, factor)
                    terms.append(nested)
                    continue

                terms.append(NO(left, term))
                continue

            terms.append(NO(left, term))
        return sp.Add(*terms)

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
    # 参考 OPEdefs.m line 1467-1471：只在 NOOrder[A,C]>0 且有非零 OPE 时展开
    if isinstance(left, NormalOrderedOperator):
        A, B = left.left, left.right
        C = right
        # 检查是否需要展开：只有当 A 应该排在 C 之后时才展开
        order_AC = ope_registry.compare_operators(A, C)
        if order_AC > 0:
            # 检查是否有非零的 OPE，如果没有就不展开
            from .api import OPE

            ope_AC = OPE(A, C)
            ope_BC = OPE(B, C) if not _operators_equal(A, B) else ope_AC
            has_nonzero_ope = ope_AC.max_pole > 0 or ope_BC.max_pole > 0

            if has_nonzero_ope:
                # A 应该排在 C 之后，且有非零 OPE，需要展开成 A(BC) 形式
                return _expand_nested_no_left(A, B, C, expand_derivatives)
            # 否则保持原结构，不展开

        # A 应该排在 C 之前或相等，或者没有非零 OPE，顺序已经正确，不需要展开
        # 但仍需递归化简内层的 NO(A,B)
        simplified_left = simplify(left, expand_derivatives)
        if simplified_left != left:
            return NO(simplified_left, right)
        return NO(left, right)

    # 如果右侧是 NO，展开为 NO(A, NO(B,C))
    # 参考 OPEdefs.m line 1473-1475：只在需要重新排序时展开
    if isinstance(right, NormalOrderedOperator):
        A, C = right.left, right.right
        # NOOrder[A,B]>0 意味着需要展开
        # 在我们的实现中，compare_operators 返回相反的符号
        order_AB = ope_registry.compare_operators(A, left)
        if order_AB < 0:
            # A 的 position 更小，应该在 left 之前，需要展开
            return _expand_nested_no_right(left, A, C, expand_derivatives)
        # 否则顺序已经正确，不需要展开
        return NO(left, simplify(right, expand_derivatives))

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


def _expand_nested_no_left(
    A: Any, B: Any, C: Any, expand_derivatives: bool = True
) -> Any:
    """
    展开 NO(NO(A,B), C) 使用 Jacobi 恒等式

    对应 OPEdefs.m 中的 NOCompositeHelpLQ

    公式 (Thielemans eq 3.3.4):
    (AB)C = A(BC) + Σ_{l=1}^∞ (1/l!) (∂^l A {BC}_l) + (-1)^{|A||B|} Σ_{l=1}^∞ (1/l!) (∂^l B {AC}_l)

    Args:
        A: NO 左侧的左算符
        B: NO 左侧的右算符
        C: NO 右侧的算符
        expand_derivatives: 是否展开导数

    Returns:
        展开后的表达式
    """
    from .api import OPE, NO
    from .operators import d as derivative_operator
    from .local_operator import get_operator_parity
    from .constants import Zero as ZeroConst

    parity_A = get_operator_parity(A)
    parity_B = get_operator_parity(B)
    sign = (-1) ** (parity_A * parity_B)

    ope_AC = OPE(A, C)
    ope_BC = OPE(B, C) if not _operators_equal(A, B) else ope_AC

    result = NO(A, NO(B, C))

    for l in range(1, ope_BC.max_pole + 1):
        pole_l = ope_BC.pole(l)
        if pole_l != 0 and pole_l != ZeroConst:
            deriv_A = derivative_operator(A, l)
            coeff = sp.Rational(1, sp.factorial(l))
            # 如果 pole_l 是算符，构造 NO；否则直接乘以标量
            if isinstance(pole_l, Operator):
                result = result + coeff * NO(deriv_A, pole_l)
            else:
                result = result + coeff * pole_l * deriv_A

    for l in range(1, ope_AC.max_pole + 1):
        pole_l = ope_AC.pole(l)
        if pole_l != 0 and pole_l != ZeroConst:
            deriv_B = derivative_operator(B, l)
            coeff = sign * sp.Rational(1, sp.factorial(l))
            # 如果 pole_l 是算符，构造 NO；否则直接乘以标量
            if isinstance(pole_l, Operator):
                result = result + coeff * NO(deriv_B, pole_l)
            else:
                result = result + coeff * pole_l * deriv_B

    # 不在这里递归调用 simplify，避免无限递归
    # 返回展开的结果，让外层的 simplify 处理
    return result


def _expand_nested_no_right(
    B: Any, A: Any, C: Any, expand_derivatives: bool = True
) -> Any:
    """
    展开 NO(B, NO(A,C)) 使用 Thielemans eq (3.3.9-3.3.10)

    对应 OPEdefs.m 中的 NOCompositeHelpRQ

    本函数实现的是将 B(AC) 改写成 A(BC) 的方向（对应 OPEdefs.m 的
    NOCompositeHelpRQ[B,A,C]），它可由 Thielemans eq (3.3.9-3.3.10) 交换 A/B 得到。

    公式:
        B(AC) = (-1)^{|A||B|} A(BC) + (NOCommuteHelp[B,A])C

    其中 NOCommuteHelp[B,A] = NO[B,A] - (-1)^{|A||B|} NO[A,B]

    这个公式来自 Thielemans 论文 eq (3.3.9-3.3.10)：
        [A[BC]_0]_0 = (-1)^{|A||B|} [B[AC]_0]_0 +
                      [[AB]_0 - (-1)^{|A||B|} [BA]_0] C]_0

    交换 A 和 B 的角色得到：
        [B[AC]_0]_0 = (-1)^{|A||B|} [A[BC]_0]_0 +
                      [[BA]_0 - (-1)^{|A||B|} [AB]_0] C]_0

    Args:
        B: NO 外层左侧的算符 (outer_left)
        A: NO 内层左侧的算符 (inner_left)
        C: NO 内层右侧的算符 (inner_right)
        expand_derivatives: 是否展开导数

    Returns:
        展开后的表达式
    """
    from .api import NO
    from .local_operator import get_operator_parity

    parity_A = get_operator_parity(A)
    parity_B = get_operator_parity(B)
    sign = (-1) ** (parity_A * parity_B)

    # B(AC) = sign * A(BC) + (NOCommuteHelp[B,A])C
    commute_term = _compute_no_commute_help(B, A, expand_derivatives=expand_derivatives)
    result = sign * NO(A, NO(B, C)) + NO(commute_term, C)

    # 不在这里递归调用 simplify，避免无限递归
    return result


def _compute_no_commute_help(A: Any, B: Any, expand_derivatives: bool = True) -> Any:
    """
    计算 NOCommuteHelp[A,B] = NO[A,B] - sign * NO[B,A]

    对应 OPEdefs.m 中的 NOCommuteHelpQ

    使用公式: -Σ_{m=1}^∞ (-1)^m / m! ∂^m {AB}_m
    注意：前面有负号！

    Args:
        A: 第一个算符
        B: 第二个算符
        expand_derivatives: 是否展开导数

    Returns:
        对易子项
    """
    from .api import OPE
    from .operators import d as derivative_operator
    from .constants import Zero as ZeroConst

    ope_AB = OPE(A, B)
    result = ZeroConst

    for m in range(1, ope_AB.max_pole + 1):
        pole_m = ope_AB.pole(m)
        if pole_m != 0 and pole_m != ZeroConst:
            deriv_pole = derivative_operator(pole_m, m)
            # 注意：OPEdefs.m 的公式是 -(-1)^m / m!，即整体有负号
            coeff = -((-1) ** m) * sp.Rational(1, sp.factorial(m))
            result = result + coeff * deriv_pole

    # 不在这里递归调用 simplify，避免无限递归
    return result


def _operators_equal(op1: Any, op2: Any) -> bool:
    """
    检查两个算符是否相等

    Args:
        op1: 第一个算符
        op2: 第二个算符

    Returns:
        是否相等
    """
    if type(op1) != type(op2):
        return False

    if isinstance(op1, BasisOperator):
        return op1.name == op2.name and op1.is_bosonic == op2.is_bosonic
    elif isinstance(op1, DerivativeOperator):
        return op1.order == op2.order and _operators_equal(op1.base, op2.base)
    elif isinstance(op1, NormalOrderedOperator):
        return _operators_equal(op1.left, op2.left) and _operators_equal(
            op1.right, op2.right
        )
    else:
        return op1 == op2
