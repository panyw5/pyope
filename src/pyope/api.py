"""
OPE 计算 API 模块

本模块实现了 VOA 中 OPE 计算的核心函数：
- OPE(A, B): 计算算符积展开 A(z)B(w)
- NO(A, B): 计算正规序乘积 (AB)(z)
- bracket(A, B, n): 计算 bracket {AB}_n(z)
- MakeOPE: 创建 OPEData 的便捷函数
"""

from typing import Union, Any, List
import sympy as sp
from sympy import Add, Mul, Number, Integer, binomial, factorial

from .operators import (
    Operator,
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    d as derivative,
)
from .ope_data import OPEData
from .constants import Zero, One
from .registry import ope_registry, OPEDefiner
from .local_operator import (
    is_local_operator,
    extract_scalar_operator,
)


class OPEComputer(OPEDefiner):
    """
    OPE 计算器类

    继承 OPEDefiner 并添加 OPE 计算功能。
    """

    def __call__(self, left: Any, right: Any) -> OPEData:
        """
        计算 OPE: OPE(A, B)

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            OPEData 实例
        """
        return _compute_ope(left, right)

    @staticmethod
    def make(data: Union[List, OPEData]) -> OPEData:
        """
        创建 OPEData: OPE.make([...])

        Args:
            data: 极点列表或 OPEData 实例

        Returns:
            OPEData 实例
        """
        if isinstance(data, OPEData):
            return data
        elif isinstance(data, list):
            return OPEData.from_list(data)
        else:
            raise TypeError(f"MakeOPE expects list or OPEData, got {type(data)}")


# 创建全局 OPE 实例
OPE = OPEComputer(ope_registry)
"""
全局 OPE 计算器

支持三种使用方式：
1. 定义 OPE: OPE[A, B] = OPEData(...)
2. 计算 OPE: OPE(A, B)
3. 创建 OPEData: OPE.make([...])

Examples:
    >>> T = BasisOperator("T", bosonic=True)
    >>> OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])
    >>> ope_result = OPE(T, T)
"""


def MakeOPE(data: Union[List, OPEData]) -> OPEData:
    """
    创建 OPEData 的便捷函数

    支持从列表创建 OPEData（Mathematica 风格）。

    Args:
        data: 极点列表或 OPEData 实例
            - 如果是列表：按从高阶到低阶排列 [pole_n, ..., pole_2, pole_1]
            - 如果是 OPEData：直接返回

    Returns:
        OPEData 实例

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> c = sp.Symbol("c")
        >>> ope = MakeOPE([c/2*One, 0, 2*T, d(T)])
        >>> # 等价于 OPEData({4: c/2*One, 2: 2*T, 1: d(T)})
    """
    if isinstance(data, OPEData):
        return data
    elif isinstance(data, list):
        return OPEData.from_list(data)
    else:
        raise TypeError(f"MakeOPE expects list or OPEData, got {type(data)}")


def _compute_ope(left: Any, right: Any) -> OPEData:
    """
    计算两个算符的 OPE（内部实现）

    实现 OPE 计算的核心算法，按照以下顺序应用规则：
    1. 处理零算符
    2. 处理线性性（加法和标量乘法）
    3. 处理导数算符
    4. 查询注册表中的基本 OPE
    5. 处理正规序算符
    6. 返回零 OPE（未定义的情况）

    Args:
        left: 左侧算符
        right: 右侧算符

    Returns:
        OPEData 实例
    """
    # 规则 1: 处理零算符
    # OPE(0, B) = 0, OPE(A, 0) = 0
    if left == 0 or left == Zero or right == 0 or right == Zero:
        return OPEData({})

    # 规则 2: 线性性 - 右侧加法
    # OPE(A, B+C) = OPE(A,B) + OPE(A,C)
    if isinstance(right, Add):
        result = OPEData({})
        for term in right.args:
            result = result + _compute_ope(left, term)
        return result

    # 规则 3: 线性性 - 左侧加法
    # OPE(A+B, C) = OPE(A,C) + OPE(B,C)
    if isinstance(left, Add):
        result = OPEData({})
        for term in left.args:
            result = result + _compute_ope(term, right)
        return result

    # 规则 4: 标量乘法 - 右侧
    # OPE(A, c*B) = c*OPE(A,B)
    if isinstance(right, Mul):
        coeff, op = extract_scalar_operator(right)
        if coeff != 1:
            return coeff * _compute_ope(left, op)

    # 规则 5: 标量乘法 - 左侧
    # OPE(c*A, B) = c*OPE(A,B)
    if isinstance(left, Mul):
        coeff, op = extract_scalar_operator(left)
        if coeff != 1:
            return coeff * _compute_ope(op, right)

    # 规则 6: 左侧导数算符
    # [∂A, B]_q = -(q-1)[A,B]_{q-1}
    if isinstance(left, DerivativeOperator):
        return _ope_derivative_left(left, right)

    # 规则 7: 右侧导数算符
    # [A, ∂B]_q = (q-1)[A,B]_{q-1} + ∂[A,B]_q
    if isinstance(right, DerivativeOperator):
        return _ope_derivative_right(left, right)

    # 规则 8: 查询注册表
    # 对于基本算符，从注册表查询 OPE
    if isinstance(left, BasisOperator) and isinstance(right, BasisOperator):
        # 首先检查算符顺序
        order = ope_registry.compare_operators(left, right)

        if order < 0:
            # 顺序错误，需要使用对称性公式
            return _ope_commute_help(left, right)

        # 顺序正确或相同，查询注册表
        ope_data = ope_registry.get_ope(left, right)
        if ope_data is not None:
            return ope_data
        # 未定义的 OPE 返回零
        return OPEData({})

    # 规则 9: 右侧正规序算符
    # OPE(A, NO(B,C)) 使用 Jacobi 恒等式
    if isinstance(right, NormalOrderedOperator):
        return _ope_composite_right(left, right)

    # 规则 10: 左侧正规序算符
    # OPE(NO(A,B), C) 使用 Jacobi 恒等式
    if isinstance(left, NormalOrderedOperator):
        return _ope_composite_left(left, right)

    # 默认：未定义的 OPE 返回零
    return OPEData({})


def _ope_derivative_left(left: DerivativeOperator, right: Any) -> OPEData:
    """
    计算左侧导数算符的 OPE

    公式：[∂^n A, B]_q = (-1)^n * (q-1)_n * [A,B]_{q-n}
    其中 (q-1)_n 是 Pochhammer 符号

    对于 base_ope 中的 pole(p)，它对 derivative_ope 中的 pole(q) 有贡献，
    其中 q = p + n（极点阶数增加）

    Args:
        left: 导数算符 ∂^n A
        right: 右侧算符 B

    Returns:
        OPEData 实例
    """
    base = left.base
    order = left.order

    # 计算基础算符的 OPE
    base_ope = _compute_ope(base, right)

    # 应用导数规则
    new_poles = {}
    for p, coeff in base_ope.poles.items():
        # 对于 base_ope 的 pole(p)，它贡献到 derivative_ope 的 pole(q)
        # 其中 q = p + order
        q = p + order

        # Pochhammer 符号: (q-1)_n = (q-1)(q-2)...(q-n)
        pochhammer = 1
        for i in range(order):
            pochhammer *= (q - 1 - i)

        new_coeff = ((-1) ** order) * pochhammer * coeff
        new_poles[q] = new_coeff

    return OPEData(new_poles)


def _ope_derivative_right(left: Any, right: DerivativeOperator) -> OPEData:
    """
    计算右侧导数算符的 OPE

    公式：[A, ∂^n B]_q = Σ_{k=0}^{n} C(n,k) * (q-1)_k * ∂^{n-k} [A,B]_{q-k}

    简化版本（n=1）：[A, ∂B]_q = (q-1)[A,B]_{q-1} + ∂[A,B]_q

    Args:
        left: 左侧算符 A
        right: 导数算符 ∂^n B

    Returns:
        OPEData 实例
    """
    base = right.base
    order = right.order

    # 计算基础算符的 OPE
    base_ope = _compute_ope(left, base)

    # 应用导数规则
    new_poles = {}

    # 对于 base_ope 中的每个极点 p，它对 derivative_ope 中的极点 p+k 有贡献
    for p, coeff in base_ope.poles.items():
        # 对于每个 k，计算 [A,B]_p 对 [A,∂^n B]_{p+k} 的贡献
        for k in range(order + 1):
            # 新的极点阶数：q = p + k
            new_q = p + k

            # 计算 Pochhammer 符号 (q-1)_k = (p+k-1)_k
            pochhammer = 1
            for i in range(k):
                pochhammer *= (new_q - 1 - i)

            # 计算二项式系数 C(n, k)
            binom_coeff = binomial(order, k)

            # 新的系数：对原系数求 (n-k) 阶导数
            if order - k > 0:
                new_coeff = binom_coeff * pochhammer * derivative(coeff, order - k)
            else:
                new_coeff = binom_coeff * pochhammer * coeff

            # 累加到结果中
            if new_q in new_poles:
                new_poles[new_q] = new_poles[new_q] + new_coeff
            else:
                new_poles[new_q] = new_coeff

    return OPEData(new_poles)


def _ope_composite_right(left: Any, right: NormalOrderedOperator) -> OPEData:
    """
    计算 OPE(A, NO(B,C))

    使用 Jacobi 恒等式：
    [A, [BC]_0]_q = (-1)^{|A||B|} [B, [AC]_q]_0 + [[AB]_q, C]_0
                    + Σ_{l=1}^{q-1} C(q-1, l) [[AB]_{q-l}, C]_l

    Args:
        left: 左侧算符 A
        right: 正规序算符 NO(B,C)

    Returns:
        OPEData 实例
    """
    A = left
    B = right.left
    C = right.right

    # 获取 parity
    parity_A = _get_parity(A)
    parity_B = _get_parity(B)
    sign = (-1) ** (parity_A * parity_B)

    result = OPEData({})

    # 计算 OPE(A, B)
    ope_AB = _compute_ope(A, B)

    # 对于每个极点 q
    for q, bracket_AB_q in ope_AB.poles.items():
        if q == 0:
            # q=0: [[AB]_0, C]_0 = NO(NO(A,B), C)
            no_AB = NO(A, B)
            no_result = NO(no_AB, C)
            result = result + OPEData({0: no_result})
        else:
            # 主项: [[AB]_q, C]_0
            # 使用 NO 函数来处理标量乘法
            no_term = NO(bracket_AB_q, C)
            result = result + OPEData({q: no_term})

            # 求和项: Σ_{l=1}^{q-1} C(q-1, l) [[AB]_{q-l}, C]_l
            for l in range(1, q):
                binom_coeff = binomial(q - 1, l)
                bracket_AB_ql = ope_AB.pole(q - l)
                if bracket_AB_ql != 0:
                    # 计算 [[AB]_{q-l}, C]_l
                    ope_bracket_C = _compute_ope(bracket_AB_ql, C)
                    pole_l = ope_bracket_C.pole(l)
                    if pole_l != 0:
                        if q in result.poles:
                            result._poles[q] = result._poles[q] + binom_coeff * pole_l
                        else:
                            result._poles[q] = binom_coeff * pole_l

    # 交换项: (-1)^{|A||B|} [B, [AC]_q]_0
    # 这需要计算所有 q 的贡献
    for q in range(1, ope_AB.max_pole + 1):
        # 计算 [AC]_q
        ope_AC = _compute_ope(A, C)
        bracket_AC_q = ope_AC.pole(q)
        if bracket_AC_q != 0:
            # 计算 [B, [AC]_q]_0 = NO(B, [AC]_q)
            # 使用 NO 函数来处理标量乘法
            no_B_AC = NO(B, bracket_AC_q)
            if q in result.poles:
                result._poles[q] = result._poles[q] + sign * no_B_AC
            else:
                result._poles[q] = sign * no_B_AC

    return result


def _ope_composite_left(left: NormalOrderedOperator, right: Any) -> OPEData:
    """
    计算 OPE(NO(A,B), C)

    使用完整的 Jacobi 恒等式公式（基于 OPEdefs.m 的实现）：

    OPE[NO[A,B], C] =
      Σ_{q=1}^{maxBC} Σ_{l=0}^{maxBC-q} NO[∂^l A, {BC}_{l+q}] / l!
      + sign * Σ_{q=1}^{maxAC} Σ_{l=0}^{maxAC-q} NO[∂^l B, {AC}_{l+q}] / l!

    其中 sign = (-1)^(|A||B|)

    Args:
        left: 正规序算符 NO(A,B)
        right: 右侧算符 C

    Returns:
        OPEData 实例
    """
    A = left.left
    B = left.right
    C = right

    # 获取 parity 和交换符号
    parity_A = _get_parity(A)
    parity_B = _get_parity(B)
    sign = (-1) ** (parity_A * parity_B)

    # 计算 OPE(A, C) 和 OPE(B, C)
    ope_AC = _compute_ope(A, C)
    ope_BC = _compute_ope(B, C)

    max_AC = ope_AC.max_pole
    max_BC = ope_BC.max_pole

    # 如果两个 OPE 都是零，直接返回零
    if max_AC == 0 and max_BC == 0:
        return OPEData({0: NormalOrderedOperator(left, right)})

    result = OPEData({})

    # 第一项: Σ_{q=1}^{maxBC} Σ_{l=0}^{maxBC-q} NO[∂^l A, {BC}_{l+q}] / l!
    for q in range(1, max_BC + 1):
        pole_sum = 0
        for l in range(0, max_BC - q + 1):
            # 计算 ∂^l A
            if l == 0:
                deriv_A = A
            else:
                deriv_A = derivative(A, l)

            # 获取 {BC}_{l+q}
            bracket_BC = ope_BC.pole(l + q)
            if bracket_BC != 0:
                # 计算 NO[∂^l A, {BC}_{l+q}] / l!
                no_term = NO(deriv_A, bracket_BC)
                pole_sum = pole_sum + no_term / factorial(l)

        if pole_sum != 0:
            result._poles[q] = pole_sum

    # 第二项: sign * Σ_{q=1}^{maxAC} Σ_{l=0}^{maxAC-q} NO[∂^l B, {AC}_{l+q}] / l!
    for q in range(1, max_AC + 1):
        pole_sum = 0
        for l in range(0, max_AC - q + 1):
            # 计算 ∂^l B
            if l == 0:
                deriv_B = B
            else:
                deriv_B = derivative(B, l)

            # 获取 {AC}_{l+q}
            bracket_AC = ope_AC.pole(l + q)
            if bracket_AC != 0:
                # 计算 NO[∂^l B, {AC}_{l+q}] / l!
                no_term = NO(deriv_B, bracket_AC)
                pole_sum = pole_sum + no_term / factorial(l)

        if pole_sum != 0:
            # 累加到结果中
            if q in result._poles:
                result._poles[q] = result._poles[q] + sign * pole_sum
            else:
                result._poles[q] = sign * pole_sum

    return result


def _get_parity(operator: Any) -> int:
    """
    获取算符的 parity

    Args:
        operator: 算符

    Returns:
        parity 值（0 或 1）
    """
    if isinstance(operator, Operator):
        return operator.parity

    # 对于复合表达式，尝试从注册表获取
    if isinstance(operator, BasisOperator):
        parity = ope_registry.get_parity(operator)
        if parity is not None:
            return parity

    # 默认返回 0（玻色子）
    return 0


def bracket(left: Any, right: Any, n: int = None, anticommutator: bool = None) -> Any:
    """
    计算 bracket {AB}_n 或对易子/反对易子

    有两种使用方式：
    1. bracket(A, B, n): 从 OPE(A, B) 中提取第 n 阶极点的系数
    2. bracket(A, B, anticommutator=False/True): 计算对易子或反对易子

    Args:
        left: 左侧算符 A
        right: 右侧算符 B
        n: 极点阶数（可选）
        anticommutator: 如果为 False 计算对易子，如果为 True 计算反对易子（可选）

    Returns:
        第 n 阶极点的系数（LocalOperator）或对易子/反对易子

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> OPE[T, T] = MakeOPE([c/2*One, 0, 2*T, d(T)])
        >>> bracket(T, T, 2)  # 返回 2*T
        >>> bracket(T, T, 0)  # 返回 NO(T, T)
        >>> bracket(T, J, anticommutator=False)  # 返回 [T, J] = NO(T,J) - NO(J,T)
    """
    if anticommutator is not None:
        # 计算对易子或反对易子
        if anticommutator:
            # {A, B} = NO(A, B) + NO(B, A)
            return NO(left, right) + NO(right, left)
        else:
            # [A, B] = NO(A, B) - NO(B, A)
            return NO(left, right) - NO(right, left)
    elif n is not None:
        # 提取第 n 阶极点
        ope_result = _compute_ope(left, right)
        return ope_result.pole(n)
    else:
        raise ValueError("Either 'n' or 'anticommutator' must be specified")


def NO(left: Any, right: Any) -> Any:
    """
    计算正规序乘积 (AB)

    正规序乘积定义为 OPE 的 0 阶极点：NO(A, B) = {AB}_0

    Args:
        left: 左侧算符 A
        right: 右侧算符 B

    Returns:
        NormalOrderedOperator 或简化后的表达式

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> J = BasisOperator("J", bosonic=True)
        >>> NO(T, J)  # 返回 NormalOrderedOperator(T, J)
    """
    # 处理零算符
    if left == 0 or left == Zero:
        return 0
    if right == 0 or right == Zero:
        return 0

    # 处理单位算符
    if left == One:
        return right
    if right == One:
        return left

    # 处理线性性
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

    # 确保 left 和 right 都是 Operator 实例
    if not isinstance(left, Operator):
        raise TypeError(f"NO requires Operator instances for left operand, got {type(left)}")
    if not isinstance(right, Operator):
        raise TypeError(f"NO requires Operator instances for right operand, got {type(right)}")

    # 创建正规序算符
    return NormalOrderedOperator(left, right)


def _ope_commute_help(left: Any, right: Any) -> OPEData:
    """
    计算 OPE(B,A) from OPE(A,B) 使用对称性公式

    类似于 OPEdefs.m 中的 OPECommuteHelp

    公式：
    [BA](q) = SwapSign[A,B] * Sum[(-1)^l / (l-q)! * D^(l-q) [[AB](l)], {l,q,MaxPole[AB]}]

    其中：
    - SwapSign[A,B] = (-1)^(|A||B|)
    - D^n 表示对 w 求 n 阶导数

    Args:
        left: 对应 OPE[B, A] 中的 B（顺序错误的左侧）
        right: 对应 OPE[B, A] 中的 A（顺序错误的右侧）

    Returns:
        OPEData 实例表示 OPE[B, A]
    """
    # 计算正确顺序的 OPE: OPE[A, B]
    ope_AB = _compute_ope(right, left)  # 注意：right 是 A，left 是 B

    if ope_AB.max_pole == 0:
        return OPEData({})

    # 计算 SwapSign
    from .local_operator import get_operator_parity
    parity_A = get_operator_parity(right)
    parity_B = get_operator_parity(left)
    swap_sign = ((-1) ** (parity_A * parity_B))

    # 应用公式
    max_pole = ope_AB.max_pole
    new_poles = {}

    for q in range(max_pole, 0, -1):
        # term[q] 累积 pole q 的系数
        term_q = swap_sign * ((-1) ** q) * ope_AB.pole(q)

        # 添加导数项的贡献
        for l in range(q + 1, max_pole + 1):
            # pole(l) 对 pole(q) 的贡献通过 (l-q) 阶导数
            pole_l = ope_AB.pole(l)
            if pole_l != 0:
                # 计算 D^(l-q) [pole_l] / (l-q)!
                deriv_order = l - q
                deriv_pole = derivative(pole_l, deriv_order)

                # 系数: swap_sign * (-1)^l / (l-q)!
                coeff = swap_sign * ((-1) ** l) / sp.factorial(deriv_order)
                term_q = term_q + coeff * deriv_pole

        if term_q != 0:
            new_poles[q] = term_q

    return OPEData(new_poles)
