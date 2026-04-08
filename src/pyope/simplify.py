"""
表达式化简模块

本模块提供将 OPE 计算结果化简为规范形式的功能。

主要函数：
- simplify(expr): 将表达式化简为排序的 NO product 线性组合
"""

from functools import lru_cache
from typing import Any, Dict, List, Tuple
import sympy as sp
from sympy import Add, Mul, Pow

from .backend import get_compute_backend
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
from .utils import _multiply_bracket_operator


def simplify(
    expr: Any, expand_derivatives: bool = True, preserve_nested_structure: bool = False
) -> Any:
    """
    化简 OPE 表达式为规范形式

    将表达式化简为排序的 NO product 的线性组合。执行以下操作：
    1. 展开嵌套的 NO 乘积（默认将左嵌套转换为右嵌套）
    2. 在 NO 内部按算符顺序排列（实现两种正规排序）
       - 算符之间按标准顺序排列：NO(C, NO(A, B)) -> NO(A, NO(B, C)) + additional terms
       - NO 化成右嵌套：NO(NO(...), ...) -> NO(..., NO(..., NO(...)))
    3. 合并同类项
    4. 标准化导数表示
    5. （可选）应用莱布尼茨法则展开正规序的导数

    Args:
        expr: 要化简的表达式（可以是 Operator、OPEData 或符号表达式）
        expand_derivatives: 是否自动展开正规序算符的导数（默认 True）
                           当为 True 时，应用莱布尼茨法则：
                           d^n(NO(A,B)) = Σ_{k=0}^{n} C(n,k) * NO(d^k(A), d^{n-k}(B))
        preserve_nested_structure: 是否保持嵌套 NO 的结构不变（默认 False）
                                  当为 False 时，自动将 NO(NO(A,B),C) 展开为右嵌套形式 NO(A,NO(B,C)) + 校正项
                                  并确保所有算符按标准顺序排列
                                  当为 True 时，仅在需要重排序时才展开嵌套 NO

    Returns:
        化简后的表达式

    Examples:
        >>> T = BasisOperator("T")
        >>> J = BasisOperator("J")
        >>> expr = NO(NO(T,J), J)
        >>> simplified = simplify(expr)  # 默认展开为右嵌套并排序

        >>> # 保持嵌套结构
        >>> simplify(expr, preserve_nested_structure=True)  # 仅在需要排序时展开

        >>> # 自动展开导数（默认行为）
        >>> from pyope import d
        >>> expr = d(NO(T, J))
        >>> simplify(expr)  # 返回 NO(d(T), J) + NO(T, d(J))

        >>> # 禁用导数展开
        >>> simplify(expr, expand_derivatives=False)  # 返回 d(NO(T, J))
    """
    backend = get_compute_backend()

    if backend == "wolfram" and _should_delegate_wolfram_simplify(expr):
        from .backend import compute_backend
        from .wolfram_backend import expand_with_wolfram

        with compute_backend("sympy"):
            evaluated = expand_with_wolfram(expr)
            if evaluated != expr:
                return simplify(
                    evaluated, expand_derivatives, preserve_nested_structure
                )

    # 处理零
    if expr == 0 or expr == Zero:
        return Zero

    # 处理单位
    if expr == One:
        return One

    # 处理 OPEData
    from .ope_data import OPEData

    if isinstance(expr, OPEData):
        return _simplify_ope_data(expr, expand_derivatives, preserve_nested_structure)

    # 处理加法：分别化简每一项
    if isinstance(expr, Add):
        simplified_terms = [
            simplify(term, expand_derivatives, preserve_nested_structure)
            for term in expr.args
        ]
        # 过滤掉 Zero 项
        non_zero_terms = [t for t in simplified_terms if t != Zero and t != 0]
        if not non_zero_terms:
            return Zero
        if len(non_zero_terms) == 1:
            return non_zero_terms[0]
        return sp.Add(*non_zero_terms)

    # 处理幂运算：Operator**n -> 嵌套的 NO（必须在 Mul 之前处理）
    if isinstance(expr, Pow):
        base = expr.base
        exp = expr.exp

        # 只处理算符的整数次幂
        if isinstance(base, Operator) and exp.is_integer and exp >= 2:
            from .api import NO

            result = base
            for _ in range(int(exp) - 1):
                result = NO(result, base)
            return simplify(result, expand_derivatives, preserve_nested_structure)

        # 对于非算符的幂，递归简化 base
        simplified_base = simplify(base, expand_derivatives, preserve_nested_structure)
        if simplified_base != base:
            return simplify(
                simplified_base**exp, expand_derivatives, preserve_nested_structure
            )
        return expr

    # 处理标量乘法
    if isinstance(expr, Mul):
        # 先递归简化所有参数（特别是 Pow）
        simplified_args = [
            simplify(arg, expand_derivatives, preserve_nested_structure)
            for arg in expr.args
        ]

        # 如果有参数被简化了，重新构造表达式并递归简化
        if simplified_args != list(expr.args):
            new_expr = sp.Mul(*simplified_args)
            return simplify(new_expr, expand_derivatives, preserve_nested_structure)

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

        simplified_op = simplify(op, expand_derivatives, preserve_nested_structure)
        # 关键：把 (-1)*Zero, k*Zero 等形式规约到 Zero
        if simplified_op == Zero or simplified_op == 0:
            return Zero
        if coeff != 1:
            return coeff * simplified_op

    # 处理算符
    if isinstance(expr, Operator):
        result = _simplify_operator(expr, expand_derivatives, preserve_nested_structure)
        return result

    # 其他情况：直接返回
    return expr


def _should_delegate_wolfram_simplify(expr: Any) -> bool:
    if isinstance(expr, (Add, NormalOrderedOperator)):
        return True
    if isinstance(expr, Mul) and expr.has(Operator):
        return True
    return False


def _simplify_ope_data(
    ope_data: Any,
    expand_derivatives: bool = True,
    preserve_nested_structure: bool = False,
) -> Any:
    """
    化简 OPEData 对象

    Args:
        ope_data: OPEData 实例
        expand_derivatives: 是否展开导数
        preserve_nested_structure: 是否保持嵌套结构

    Returns:
        化简后的 OPEData
    """
    from .ope_data import OPEData

    new_poles = {}
    for q, coeff in ope_data.poles.items():
        simplified_coeff = simplify(
            coeff, expand_derivatives, preserve_nested_structure
        )
        if simplified_coeff != 0 and simplified_coeff != Zero:
            new_poles[q] = simplified_coeff

    return OPEData(new_poles)


def _simplify_operator(
    op: Operator,
    expand_derivatives: bool = True,
    preserve_nested_structure: bool = False,
) -> Any:
    """
    化简单个算符

    Args:
        op: 算符实例
        expand_derivatives: 是否展开导数
        preserve_nested_structure: 是否保持嵌套结构

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
            left = simplify(base.left, expand_derivatives, preserve_nested_structure)
            right = simplify(base.right, expand_derivatives, preserve_nested_structure)

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
            return simplify(result, expand_derivatives, preserve_nested_structure)

        # 对于非 NO 的 base，保持 DerivativeOperator 结构
        # 但可以递归化简 base（可选）
        simplified_base = simplify(base, expand_derivatives, preserve_nested_structure)
        if simplified_base != base:
            return DerivativeOperator(simplified_base, n)
        return op

    # BasisOperator 和 DerivativeOperator（未展开）已经是最简形式
    if isinstance(op, (BasisOperator, DerivativeOperator)):
        return op

    # 处理 NormalOrderedOperator
    if isinstance(op, NormalOrderedOperator):
        return _simplify_normal_ordered(
            op, expand_derivatives, preserve_nested_structure
        )

    return op


def _simplify_normal_ordered(
    no_op: NormalOrderedOperator,
    expand_derivatives: bool = True,
    preserve_nested_structure: bool = False,
) -> Any:
    """
    化简正规序 NO(..., ...) 算符

    其中 ... 本身可能也是 NO 算符

    处理：
    1. 递归化简左右算符
    2. 对 q = 0 的正规序重排应用 Thielemans / OPEdefs 规则
    3. 展开嵌套的 NO（根据 preserve_nested_structure 参数）
    4. 排序内部算符

    说明：
    - 本函数只负责 q = 0 的正规序重排与规范化。
    - q >= 1 的奇异 OPE 计算在 api.py 中完成，对应 manual eqs. (9), (10), (18)。
    - 这里使用的 q = 0 规则对应 manual eqs. (14), (15), (16), (19)，
      以及 OPEdefs.m 中的 NOCommuteHelpQ / NOCompositeHelpLQ / NOCompositeHelpRQ。

    Args:
        no_op: NormalOrderedOperator 实例
        expand_derivatives: 是否展开导数
        preserve_nested_structure: 是否保持嵌套结构（默认 False，即默认展开左嵌套）

    Returns:
        化简后的表达式
    """
    from .api import NO

    # 递归化简左右算符
    left = simplify(no_op.left, expand_derivatives, preserve_nested_structure)
    right = simplify(no_op.right, expand_derivatives, preserve_nested_structure)

    # <=================>
    # ====== 分配律 ======
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
        # 关键：分配后生成的新 NO 项还需要进一步规范化（重排/展开嵌套等）
        # 否则在 Jacobi 检查等场景会残留“形式不一样但数学上为零”的表达式。
        return simplify(sp.Add(*terms), expand_derivatives, preserve_nested_structure)

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
        # 同上：分配后需要再跑一遍 simplify 以完成规范化。
        return simplify(sp.Add(*terms), expand_derivatives, preserve_nested_structure)

    # 处理标量乘法
    if isinstance(left, Mul):
        coeff, op = extract_scalar_operator(left)
        if coeff != 1:
            return coeff * NO(op, right)
    if isinstance(right, Mul):
        coeff, op = extract_scalar_operator(right)
        if coeff != 1:
            return coeff * NO(left, op)

    # ===== 分配律结束 ======
    # <====================>

    # <========================>
    # ====== 嵌套 NO 展开 =======

    # 处理嵌套的 NO
    # 如果左侧是 NO，则利用式子
    # NO(NO(A, B), C) = NO(A, NO(B, C))
    #                   + sum_{l=1}^maxBC 1/l NO(∂^l A, {BC}_l)
    #                   + sum_{l=1}^maxAC sign * 1/l NO(∂^l B, {AC}_l)
    if isinstance(left, NormalOrderedOperator):
        A, B = left.left, left.right
        C = right

        if not preserve_nested_structure:
            # 默认规范化路径：对左嵌套应用 q = 0 的左复合正规序公式。
            # 对应 manual eq. (19) / OPEdefs.m 的 NOCompositeHelpLQ。
            from .api import OPE

            # 对应 OPEdefs.m 的 maxAC, maxBC
            ope_AC = OPE(A, C)
            ope_BC = OPE(B, C) if not _operators_equal(A, B) else ope_AC
            has_nonzero_ope = ope_AC.max_pole > 0 or ope_BC.max_pole > 0

            if has_nonzero_ope:
                # 只有当 OPE(A,C) 或 OPE(B,C) 含奇异部分时，eq. (19) 的修正项才非零。
                expanded = _expand_nested_no_left(A, B, C, expand_derivatives)
                # 递归简化展开后的结果
                return simplify(expanded, expand_derivatives, preserve_nested_structure)
            else:
                # 当 OPE(A,C)=OPE(B,C)=0 时，eq. (19) 只剩首项：
                # [[AB]_0 C]_0 = [A [BC]_0]_0。
                # 对应 OPEdefs.m 的 NOCompositeHelpLQ 首项，不能额外乘整体 SwapSign。
                # 递归简化转换后的结果，因为 B 或 C 可能仍然是嵌套的 NO
                result = NO(A, NO(B, C))
                return simplify(result, expand_derivatives, preserve_nested_structure)
        else:
            order_AC = ope_registry.compare_operators(A, C)
            if order_AC < 0:
                # A 应该排在 C 之后，需要展开把 C 移到前面
                from .api import OPE

                ope_AC = OPE(A, C)
                ope_BC = OPE(B, C) if not _operators_equal(A, B) else ope_AC
                has_nonzero_ope = ope_AC.max_pole > 0 or ope_BC.max_pole > 0

                if has_nonzero_ope:
                    # 有非零 OPE，需要用 Jacobi 恒等式展开
                    expanded = _expand_nested_no_left(A, B, C, expand_derivatives)
                    # 递归简化展开后的结果
                    return simplify(
                        expanded, expand_derivatives, preserve_nested_structure
                    )
                # 如果 A 与 C、B 与 C 都没有奇异 OPE，但又必须交换次序，
                # 这里只剩 super-commutativity 的整体 parity 符号。
                from .local_operator import get_operator_parity

                parity_AB = (get_operator_parity(A) + get_operator_parity(B)) % 2
                parity_C = get_operator_parity(C)
                sign = (-1) ** (parity_AB * parity_C)
                return sign * NO(C, NO(A, B))

            # left 在入口处已经按同配置完成 simplify。
            return NO(left, right)

    # 如果右侧是 NO，展开为 NO(B, NO(A,C))
    # 需要确保所有算符都按标准顺序排列
    if isinstance(right, NormalOrderedOperator):
        A, C = right.left, right.right
        B = left

        # 检查整体是否满足标准序
        # NO(B, NO(A, C)) 的扁平化形式是 [B, A, C]
        # 需要检查所有相邻对的顺序

        # 1. 检查 B 和 A 的顺序
        order_BA = ope_registry.compare_operators(B, A)
        # 2. 检查 A 和 C 的顺序（内层 NO）
        order_AC = ope_registry.compare_operators(A, C)

        # 根据 compare_operators 的语义：
        # > 0: 顺序正确（第一个算符应该在前）
        # = 0: 两个算符相等
        # < 0: 需要交换（第一个算符应该在后）

        # 情况 1: B 和 A 需要交换（order_BA < 0），即 A 应该排在 B 前面
        if order_BA < 0:
            # A 应该排在 B 前面，使用 _expand_nested_no_right 展开
            expanded = _expand_nested_no_right(B, A, C, expand_derivatives)
            # 递归简化展开后的结果
            return simplify(expanded, expand_derivatives, preserve_nested_structure)

        # 情况 2: B 和 A 的顺序正确，但 A 和 C 需要交换（order_AC < 0）
        # 需要先简化内层 NO(A, C)，这会触发内层的重排
        if order_AC < 0:
            # right 在入口处已经按同配置完成 simplify；若仍不满足次序，
            # 直接把当前组合重新走一遍整体规范化路径即可。
            return simplify(
                NO(B, right),
                expand_derivatives,
                preserve_nested_structure,
            )

        # 顺序已经正确，不需要再次递归化简内层
        return NO(left, right)

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

            # FIX: l=0 项是 NO(A,B)，不是 ope_AB.pole(0)
            # 因为 OPE(A,B) 只包含奇异部分（l>=1），0 阶由 NO 定义
            result = swap_sign * NO(right, left)

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
            return simplify(result, expand_derivatives, preserve_nested_structure)

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
        if isinstance(op, Add):
            distributed = sp.Add(*[coeff * term for term in op.args])
            return collect_normal_ordered_terms(distributed)
        if isinstance(op, NormalOrderedOperator):
            key = _make_no_key(op)
            terms[key] = terms.get(key, 0) + coeff
        else:
            # 非 NO 项
            key = ("other", op)
            terms[key] = terms.get(key, 0) + coeff
    elif isinstance(expr, NormalOrderedOperator):
        key = _make_no_key(expr)
        terms[key] = terms.get(key, 0) + 1
    else:
        # 其他项
        key = ("other", expr)
        terms[key] = terms.get(key, 0) + 1

    return terms


def combine_normal_ordered_terms(expr: Any) -> Any:
    """
    合并表达式中多个 NO 项的线性组合。

    这个函数会把结构完全相同的正规序项视为同类项，并把它们的系数相加；
    非 NO 项也会按表达式本身合并。

    Args:
        expr: 输入表达式

    Returns:
        合并同类项后的表达式

    Examples:
        >>> expr = 2*NO(T, J) + 3*NO(T, J) - NO(J, J)
        >>> combine_normal_ordered_terms(expr)
        5*NO(T, J) - NO(J, J)
    """
    terms = collect_normal_ordered_terms(expr)

    rebuilt_terms = []
    for key, coeff in terms.items():
        coeff = sp.sympify(coeff)
        if coeff == 0:
            continue

        operator = _key_to_expr(key)
        if operator == One:
            rebuilt_terms.append(coeff)
        elif coeff == 1:
            rebuilt_terms.append(operator)
        else:
            rebuilt_terms.append(coeff * operator)

    if not rebuilt_terms:
        return Zero
    if len(rebuilt_terms) == 1:
        return rebuilt_terms[0]
    return sp.Add(*rebuilt_terms)


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
        return ("Other", op)


def _key_to_expr(key: Tuple) -> Any:
    """将 collect_normal_ordered_terms 生成的键还原为表达式。"""
    from .api import NO

    key_type = key[0]

    if key_type == "NO":
        return NO(_key_to_expr(key[1]), _key_to_expr(key[2]))
    if key_type == "Basis":
        name, is_bosonic = key[1], key[2]
        return BasisOperator(name, fermionic=not is_bosonic)
    if key_type == "Deriv":
        return DerivativeOperator(_key_to_expr(key[1]), key[2])
    if key_type in {"Other", "other"}:
        return key[1]

    raise ValueError(f"Unknown normal-ordered term key type: {key_type}")


def flatten_right_no(expr: Any) -> list:
    """
    将右结合的 NO 扁平化为算符列表

    例如：NO(a, NO(b, NO(c, d))) -> [a, b, c, d]

    Args:
        expr: 表达式（可能是 NormalOrderedOperator 或其他算符）

    Returns:
        算符列表
    """
    if not isinstance(expr, NormalOrderedOperator):
        return [expr]

    left = expr.left
    right = expr.right
    return [left] + flatten_right_no(right)


def is_standard_order_no(expr: Any) -> bool:
    """
    检查 NormalOrderedOperator 是否满足标准序（canonical order）

    标准序要求：
    1. 所有相邻算符对都满足 compare(f_i, f_{i+1}) >= 0
    2. 即 f_i 应该排在 f_{i+1} 前面或相等

    Args:
        expr: NormalOrderedOperator 实例

    Returns:
        True 如果满足标准序，False 否则
    """
    if not isinstance(expr, NormalOrderedOperator):
        return True  # 非 NO 算符总是满足标准序

    factors = flatten_right_no(expr)

    # 检查所有相邻对
    for i in range(len(factors) - 1):
        order = ope_registry.compare_operators(factors[i], factors[i + 1])
        if order < 0:
            # factors[i] 应该排在 factors[i+1] 后面，不满足标准序
            return False

    return True


def _expand_nested_no_left(
    A: Any, B: Any, C: Any, expand_derivatives: bool = True
) -> Any:
    """
    展开 嵌套正规序算符 NO(NO(A,B), C)

    根据 Thielemans 文献,将 NO(NO(A,B), C) 重构为规范形式:

        NO(NO(A,B), C) → NO(A, NO(B,C)) + 校正项

    展开公式 (对应 OPEdefs.m 中的 NOCompositeHelpLQ):

        NO(A, NO(B,C)) +
        Σ_{l=1}^{maxBC} 1/l! · NO(∂^l A, {BC}_l) +
        Σ_{l=1}^{maxAC} sign · 1/l! · NO(∂^l B, {AC}_l)

    其中:
        sign = (-1)^{|A||B|} (由 A 和 B 的奇偶性决定的符号因子)
        {BC}_l = pole(OPE(B,C), l) 是 OPE(B,C) 的 l 阶极点系数
        {AC}_l = pole(OPE(A,C), l) 是 OPE(A,C) 的 l 阶极点系数
        maxBC = OPE(B,C).max_pole 是 OPE(B,C) 的最大极点数
        maxAC = OPE(A,C).max_pole 是 OPE(A,C) 的最大极点数
        ∂^l A 表示 A 的 l 阶导数

    与 Thielemans 方程的关系:
    - 该公式对应 Thielemans eq (3.3.19)
    - eq (3.3.19) 是 eq (2.3.24) 在 q=0 时的特殊情况，结合 eq (3.3.6) 得到
    - 这是正规序乘积 (normal ordered product) 情形下的 Jacobi 恒等式
    - 注意：eq (3.3.4) 是针对 [A[BC]_0]_q (q ≥ 1) 的情况，不适用于此函数

    Args:
        A: NO 左侧的左算符 (NO(A,B) 中的 A)
        B: NO 左侧的右算符 (NO(A,B) 中的 B)
        C: NO 右侧的算符
        expand_derivatives: 是否展开导数 (当前未使用,保留以保持接口一致性)

    Returns:
        Any: 展开后的表达式,为 NO(A, NO(B,C)) 及其校正项的和

    Note:
        该函数不会递归调用 simplify,以避免无限递归。
        外层调用 simplify 时会继续处理展开后的结果。

    See also:
        _expand_nested_no_right: 处理 NO(A, NO(B,C)) 的展开
        _simplify_normal_ordered: 调用此函数的化简函数
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

            # 根据 Thielemans eq (3.3.19)，所有项都应该是 NO(∂^l A, {BC}_l)
            # 即使 pole_l 不是单个算符，也应该用 NO 包裹
            result = result + coeff * NO(deriv_A, pole_l)

    for l in range(1, ope_AC.max_pole + 1):
        pole_l = ope_AC.pole(l)
        if pole_l != 0 and pole_l != ZeroConst:
            deriv_B = derivative_operator(B, l)
            coeff = sign * sp.Rational(1, sp.factorial(l))

            # 根据 Thielemans eq (3.3.19)，所有项都应该是 NO(∂^l B, {AC}_l)
            # 即使 pole_l 不是单个算符，也应该用 NO 包裹
            result = result + coeff * NO(deriv_B, pole_l)

    # 不在这里递归调用 simplify，避免无限递归
    # 返回展开的结果，让外层的 simplify 处理
    return result


def _expand_nested_no_right(
    B: Any, A: Any, C: Any, expand_derivatives: bool = True
) -> Any:
    """
    展开 NO(B, NO(A,C)) 使用 q = 0 的右复合正规序公式。

    对应 OPEdefs.m 中的 NOCompositeHelpRQ

    本函数实现的是将 B(AC) 改写成 A(BC) 的方向（对应 OPEdefs.m 的
    NOCompositeHelpRQ[B,A,C]），它可由 Thielemans eq (3.3.9-3.3.10) 交换 A/B 得到。

    公式:
        NO(B, NO(A,C)) = (-1)^{|A||B|} NO(A, NO(B,C)) + NO(NOCommuteHelp[B,A], C)

    其中 NOCommuteHelp[B,A] = NO[B,A] - (-1)^{|A||B|} NO[A,B]

    这个公式对应 manual eqs. (15)-(16)（论文中对应 q = 0 Jacobi 重排）：
        eq (3.3.9): [A[BC]_0]_0 = (-1)^{|A||B|} [B[AC]_0]_0 + ...
        eq (3.3.10): ... + [([AB]_0 - (-1)^{|A||B|} [BA]_0) C]_0

    交换 A 和 B 的角色得到：
        [B[AC]_0]_0 = (-1)^{|A||B|} [A[BC]_0]_0 + [([BA]_0 - (-1)^{|A||B|} [AB]_0) C]_0

    与 Thielemans 方程的关系:
    - eq (3.3.9) 和 (3.3.10) 合起来给出正规序乘积的重排公式
    - eq (3.3.10) 来自 eq (2.3.21)
    - 这是用于处理 NO(B, NO(A,C)) 形式的嵌套正规序

    Args:
        B: NO 外层左侧的算符 (outer_left)
        A: NO 内层左侧的算符 (inner_left)
        C: NO 内层右侧的算符 (inner_right)
        expand_derivatives: 是否展开导数 (当前未使用,保留以保持接口一致性)

    Returns:
        Any: 展开后的表达式

    Note:
        该函数不会递归调用 simplify,以避免无限递归。
        外层调用 simplify 时会继续处理展开后的结果。

    See also:
        _expand_nested_no_left: 处理 NO(NO(A,B), C) 的展开
        _compute_no_commute_help: 计算 NOCommuteHelp[A,B]
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
    return _compute_no_commute_help_cached(
        A, B, expand_derivatives, ope_registry.version
    )


@lru_cache(maxsize=None)
def _compute_no_commute_help_cached(
    A: Any, B: Any, expand_derivatives: bool, registry_version: int
) -> Any:
    del registry_version
    """
    计算 NOCommuteHelp[A,B] = NO[A,B] - sign * NO[B,A]

    对应 OPEdefs.m 中的 NOCommuteHelpQ，也即 manual eq. (14)
    里除去首项 super-commutative 部分后的“导数修正项”。
    它只适用于 q = 0 的正规序交换，不用于 q >= 1 的 OPE。

    使用公式: -Σ_{m=1}^∞ (-1)^m / m! ∂^m {AB}_m
    注意：前面有负号。

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


def expand_nested_no(
    expr: Any, max_depth: Any = None, expand_derivatives: bool = False
) -> Any:
    """
    展开嵌套正规序算符，但不进行进一步的规范化

    这个函数实现 OPEdefs.m 风格的展开：应用 Thielemans 重排公式
    将嵌套的 NO 展开为和式，但保留三重嵌套项，不继续规约到二重项。

    与 simplify() 的区别：
    - simplify(): 完全规范化，将所有嵌套 NO 规约到"标准基底"（通常是二重 NO + 导数项）
    - expand_nested_no(): 只展开指定深度，保留中间形式（如三重嵌套 NO）

    Args:
        expr: 要展开的表达式
        max_depth: 最大展开深度（默认 None，表示完全展开但不规约三重项）
        expand_derivatives: 是否展开导数（默认 False，保持 d(NO(...)) 形式）

    Returns:
        展开后的表达式

    Examples:
        >>> J_plus = BasisOperator("J⁺", conformal_weight=1)
        >>> J_zero = BasisOperator("J⁰", conformal_weight=1)
        >>> J_minus = BasisOperator("J⁻", conformal_weight=1)
        >>> # 设置 OPE...
        >>> expr = NO(J_minus, NO(J_plus, J_zero))
        >>> result = expand_nested_no(expr)
        >>> # 返回: NO(J⁺, NO(J⁰, J⁻)) + 2*NO(J⁺, ∂J⁻) - NO(∂J⁰, J⁰)
        >>> # 而不是 simplify 的结果: 2*NO(J⁺, ∂J⁻)
    """
    # 处理零和单位
    if expr == 0 or expr == Zero:
        return Zero
    if expr == One:
        return One

    # 处理 OPEData
    from .ope_data import OPEData

    if isinstance(expr, OPEData):
        new_poles = {}
        for q, coeff in expr.poles.items():
            expanded_coeff = expand_nested_no(coeff, max_depth, expand_derivatives)
            if expanded_coeff != 0 and expanded_coeff != Zero:
                new_poles[q] = expanded_coeff
        return OPEData(new_poles)

    # 处理加法：分别展开每一项
    if isinstance(expr, Add):
        expanded_terms = [
            expand_nested_no(term, max_depth, expand_derivatives) for term in expr.args
        ]
        non_zero_terms = [t for t in expanded_terms if t != Zero and t != 0]
        if not non_zero_terms:
            return Zero
        if len(non_zero_terms) == 1:
            return non_zero_terms[0]
        return sp.Add(*non_zero_terms)

    # 处理乘法：提取标量系数
    if isinstance(expr, Mul):
        coeff, op = extract_scalar_operator(expr)
        expanded_op = expand_nested_no(op, max_depth, expand_derivatives)
        if expanded_op == Zero or expanded_op == 0:
            return Zero
        if coeff != 1:
            return coeff * expanded_op
        return expanded_op

    # 处理算符
    if isinstance(expr, Operator):
        return _expand_nested_no_operator(expr, max_depth, expand_derivatives)

    # 其他情况：直接返回
    return expr


def _expand_nested_no_operator(
    op: Operator, max_depth: Any = None, expand_derivatives: bool = False
) -> Any:
    """
    展开单个算符中的嵌套正规序

    Args:
        op: 算符实例
        max_depth: 最大展开深度（None 表示完全展开）
        expand_derivatives: 是否展开导数

    Returns:
        展开后的表达式
    """
    from .api import NO

    # 处理 DerivativeOperator
    if isinstance(op, DerivativeOperator):
        if expand_derivatives and isinstance(op.base, NormalOrderedOperator):
            from .operators import d as derivative_operator
            from .cache import cached_binomial

            base = op.base
            n = op.order
            left = expand_nested_no(base.left, max_depth, expand_derivatives)
            right = expand_nested_no(base.right, max_depth, expand_derivatives)

            terms = []
            for k in range(n + 1):
                coeff = cached_binomial(n, k)
                left_deriv = derivative_operator(left, k) if k > 0 else left
                right_deriv = derivative_operator(right, n - k) if n - k > 0 else right
                terms.append(coeff * NO(left_deriv, right_deriv))

            result = sp.Add(*terms) if len(terms) > 1 else terms[0]
            return expand_nested_no(result, max_depth, expand_derivatives)
        else:
            expanded_base = expand_nested_no(op.base, max_depth, expand_derivatives)
            if expanded_base != op.base:
                return DerivativeOperator(expanded_base, op.order)
            return op

    if isinstance(op, BasisOperator):
        return op

    # 处理 NormalOrderedOperator
    if isinstance(op, NormalOrderedOperator):
        if max_depth is not None and max_depth <= 0:
            return op

        left = op.left
        right = op.right

        # 计算下一层的深度
        next_depth = None if max_depth is None else max_depth - 1

        # 递归展开左右算符
        left = expand_nested_no(left, next_depth, expand_derivatives)
        right = expand_nested_no(right, next_depth, expand_derivatives)

        # 处理左嵌套：NO(NO(A,B), C)
        if isinstance(left, NormalOrderedOperator):
            A, B = left.left, left.right
            C = right

            expanded = _expand_nested_no_left(A, B, C, expand_derivatives=False)
            return expand_nested_no(expanded, max_depth, expand_derivatives)

        # 处理右嵌套：NO(B, NO(A,C))
        if isinstance(right, NormalOrderedOperator):
            B = left
            A, C = right.left, right.right

            order_AB = ope_registry.compare_operators(A, B)
            if order_AB > 0:
                expanded = _expand_nested_no_right(B, A, C, expand_derivatives=False)
                return expand_nested_no(expanded, max_depth, expand_derivatives)

            # 即使不需要重排 B 和 A，也要检查内层 NO(A,C) 是否需要重排
            order_AC = ope_registry.compare_operators(A, C)
            if order_AC < 0:
                # A 应该排在 C 后面，需要重排内层
                # 使用完整的 OPE 公式：NO(A,C) -> NO(C,A) + 收缩项
                # 然后将结果放入外层：NO(B, NO(A,C)) -> NO(B, NO(C,A)) + NO(B, 收缩项)
                from .api import OPE
                from .operators import d as derivative_operator
                from .local_operator import get_operator_parity
                from .constants import Zero as ZeroConst

                parity_A = get_operator_parity(A)
                parity_C = get_operator_parity(C)
                swap_sign = (-1) ** (parity_A * parity_C)

                ope_CA = OPE(C, A)
                max_pole = ope_CA.max_pole

                # 构建重排后的表达式
                terms = []

                # 第一项：NO(B, NO(C,A))（三重嵌套项，总是存在）
                terms.append(swap_sign * NO(B, NO(C, A)))

                # 收缩项：NO(B, ∂^l [CA]_l)
                for l in range(1, max_pole + 1):
                    pole_l = ope_CA.pole(l)
                    if pole_l != 0 and pole_l != ZeroConst:
                        deriv_pole = derivative_operator(pole_l, l)
                        coeff = swap_sign * ((-1) ** l) / sp.factorial(l)
                        terms.append(coeff * NO(B, deriv_pole))

                result = sp.Add(*terms) if len(terms) > 1 else terms[0]
                return result

        return NO(left, right)

    return op
