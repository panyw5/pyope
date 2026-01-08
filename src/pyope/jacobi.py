"""
Jacobi 恒等式检查函数

实现 OPEdefs.m 中的 OPEJacobi 函数，用于验证 Jacobi 恒等式。

Jacobi 恒等式：
    [A, [B, C]_q]_m - (-1)^(|A||B|) [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p} = 0

其中：
- [A, B]_n 表示 OPE(A, B) 的第 n 阶极点
- |A| 表示算符 A 的 parity
- C(n, k) 表示二项式系数

参考：
- OPEdefs.m 第 1601-1637 行：OPEJacobi 实现
"""

from typing import Any, List
import sympy as sp
from sympy import binomial

from .api import OPE, bracket
from .operators import Operator
from .ope_data import OPEData
from .local_operator import get_operator_parity


def check_jacobi_identity(A: Any, B: Any, C: Any, simplify_func=None) -> List[List[Any]]:
    """
    检查三个算符的 Jacobi 恒等式

    计算 Jacobi 恒等式的左侧，如果恒等式成立，结果应该全为 0。

    Jacobi 恒等式：
        [A, [B, C]_q]_m - sign * [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p} = 0

    其中 sign = (-1)^(|A||B|)

    Args:
        A: 第一个算符
        B: 第二个算符
        C: 第三个算符
        simplify_func: 简化函数（可选），默认为 sp.expand

    Returns:
        二维列表，每个元素是 Jacobi 恒等式在特定 (m, n) 处的值
        如果恒等式成立，所有元素应该为 0

    Examples:
        >>> from pyope import BasisOperator, check_jacobi_identity
        >>> import sympy as sp
        >>> c = sp.Symbol('c')
        >>> T = BasisOperator("T", bosonic=True, conformal_weight=2)
        >>> OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])
        >>> result = check_jacobi_identity(T, T, T)
        >>> # result 应该是全零矩阵

    Reference:
        OPEdefs.m lines 1601-1637: OPEJacobi implementation
    """
    if simplify_func is None:
        simplify_func = sp.expand

    # 计算 parity 和符号
    parity_A = get_operator_parity(A)
    parity_B = get_operator_parity(B)
    sign = (-1) ** (parity_A * parity_B)

    # 计算基本 OPE
    ope_AB = OPE(A, B)
    ope_BC = OPE(B, C)

    max_AB = ope_AB.max_pole
    max_BC = ope_BC.max_pole

    # 计算 AnBC[n] = OPE[A, {BC}_n]
    AnBC = {}
    for n in range(1, max_BC + 1):
        bracket_BC_n = bracket(B, C, n)
        if bracket_BC_n != 0:
            AnBC[n] = OPE(A, bracket_BC_n)
        else:
            AnBC[n] = OPEData({})

    # 计算 ABnC[n] = OPE[{AB}_n, C]
    ABnC = {}
    for n in range(1, max_AB + 1):
        bracket_AB_n = bracket(A, B, n)
        if bracket_AB_n != 0:
            ABnC[n] = OPE(bracket_AB_n, C)
        else:
            ABnC[n] = OPEData({})

    # 计算 BnAC[n] = OPE[B, {AC}_n]
    # 如果 A == B，则 BnAC = AnBC, AC = BC
    if A == B:
        BnAC = AnBC
        ope_AC = ope_BC
        max_AC = max_BC
    else:
        ope_AC = OPE(A, C)
        max_AC = ope_AC.max_pole
        BnAC = {}
        for n in range(1, max_AC + 1):
            bracket_AC_n = bracket(A, C, n)
            if bracket_AC_n != 0:
                BnAC[n] = OPE(B, bracket_AC_n)
            else:
                BnAC[n] = OPEData({})

    # 计算最大极点
    max_AnBC = max([AnBC[n].max_pole for n in range(1, max_BC + 1)] + [0])
    max_BnAC = max([BnAC[n].max_pole for n in range(1, max_AC + 1)] + [0])
    max_ABnC = max([ABnC[n].max_pole for n in range(1, max_AB + 1)] + [0])

    max_n = max(max_AnBC, max_AC, max_AB)
    max_m = max(max_BC, max_BnAC, max_ABnC)

    # 构建结果矩阵
    result = []
    for m in range(1, max_m + 1):
        row = []
        for n in range(1, max_n + 1):
            # 第一项: {A, {BC}_m}_n
            term1 = 0
            if m in AnBC:
                term1 = AnBC[m].pole(n)

            # 第二项: -sign * {B, {AC}_n}_m
            term2 = 0
            if n in BnAC:
                term2 = sign * BnAC[n].pole(m)

            # 第三项: Σ_p C(n-1, p-1) {{AB}_p, C}_{m+n-p}
            term3 = 0
            for p in range(1, n + 1):
                if p in ABnC:
                    pole_index = m + n - p
                    pole_value = ABnC[p].pole(pole_index)
                    if pole_value != 0:
                        binom_coeff = binomial(n - 1, p - 1)
                        term3 = term3 + binom_coeff * pole_value

            # Jacobi 恒等式
            jacobi_value = term1 - term2 - term3

            # 简化
            if simplify_func is not None and jacobi_value != 0:
                jacobi_value = simplify_func(jacobi_value)

            row.append(jacobi_value)
        result.append(row)

    return result


def verify_jacobi_identity(A: Any, B: Any, C: Any, simplify_func=None) -> bool:
    """
    验证 Jacobi 恒等式是否成立

    Args:
        A: 第一个算符
        B: 第二个算符
        C: 第三个算符
        simplify_func: 简化函数（可选）

    Returns:
        True 如果 Jacobi 恒等式成立（所有项为 0），否则 False

    Examples:
        >>> from pyope import BasisOperator, verify_jacobi_identity
        >>> import sympy as sp
        >>> c = sp.Symbol('c')
        >>> T = BasisOperator("T", bosonic=True, conformal_weight=2)
        >>> OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])
        >>> verify_jacobi_identity(T, T, T)
        True
    """
    result = check_jacobi_identity(A, B, C, simplify_func)

    # 检查所有元素是否为 0
    for row in result:
        for value in row:
            if value != 0:
                return False

    return True
