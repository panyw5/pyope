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


def simplify(expr: Any) -> Any:
    """
    化简 OPE 表达式为规范形式

    将表达式化简为排序的 NO product 的线性组合。执行以下操作：
    1. 展开嵌套的 NO 乘积
    2. 在 NO 内部按算符顺序排列
    3. 合并同类项
    4. 标准化导数表示

    Args:
        expr: 要化简的表达式（可以是 Operator、OPEData 或符号表达式）

    Returns:
        化简后的表达式

    Examples:
        >>> T = BasisOperator("T", bosonic=True)
        >>> J = BasisOperator("J", bosonic=True)
        >>> expr = NO(NO(T,J), J)
        >>> simplified = simplify(expr)
    """
    # 处理零
    if expr == 0 or expr == Zero:
        return 0

    # 处理单位
    if expr == One:
        return One

    # 处理 OPEData
    from .ope_data import OPEData
    if isinstance(expr, OPEData):
        return _simplify_ope_data(expr)

    # 处理加法：分别化简每一项
    if isinstance(expr, Add):
        simplified_terms = [simplify(term) for term in expr.args]
        return sp.Add(*simplified_terms)

    # 处理标量乘法
    if isinstance(expr, Mul):
        coeff, op = extract_scalar_operator(expr)
        if coeff != 1:
            return coeff * simplify(op)

    # 处理算符
    if isinstance(expr, Operator):
        return _simplify_operator(expr)

    # 其他情况直接返回
    return expr


def _simplify_ope_data(ope_data: 'OPEData') -> 'OPEData':
    """
    化简 OPEData 对象

    Args:
        ope_data: OPEData 实例

    Returns:
        化简后的 OPEData
    """
    from .ope_data import OPEData

    new_poles = {}
    for q, coeff in ope_data.poles.items():
        simplified_coeff = simplify(coeff)
        if simplified_coeff != 0:
            new_poles[q] = simplified_coeff

    return OPEData(new_poles)


def _simplify_operator(op: Operator) -> Any:
    """
    化简单个算符

    Args:
        op: 算符实例

    Returns:
        化简后的表达式
    """
    # BasisOperator 和 DerivativeOperator 已经是最简形式
    if isinstance(op, (BasisOperator, DerivativeOperator)):
        return op

    # 处理 NormalOrderedOperator
    if isinstance(op, NormalOrderedOperator):
        return _simplify_normal_ordered(op)

    return op


def _simplify_normal_ordered(no_op: NormalOrderedOperator) -> Any:
    """
    化简正规序算符

    处理：
    1. 递归化简左右算符
    2. 展开嵌套的 NO
    3. 排序内部算符

    Args:
        no_op: NormalOrderedOperator 实例

    Returns:
        化简后的表达式
    """
    from .api import NO

    # 递归化简左右算符
    left = simplify(no_op.left)
    right = simplify(no_op.right)

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
    if isinstance(left, (BasisOperator, DerivativeOperator)) and \
       isinstance(right, (BasisOperator, DerivativeOperator)):
        
        order = ope_registry.compare_operators(left, right)
        if order < 0:
            # 需要交换顺序: NO(B, A) -> NO(A, B) + 修正项
            # 修正项来自 OPE(B, A) 的极点部分
            # 公式: NO(B, A) = (-1)^{|A||B|} NO(A, B) + \sum_{n >= 1} rac{1}{n!} \partial^n \{BA\}_n
            
            # 1. 计算符号因子 (-1)^{|A||B|}
            parity_sign = 1
            if left.parity == 1 and right.parity == 1:
                parity_sign = -1
                
            # 2. 计算修正项
            # 需要计算 OPE(B, A)
            from .api import OPE
            from .operators import d as derivative_operator
            
            # 注意：这里 left 是 B，right 是 A
            # 我们要计算 OPE(left, right)
            try:
                ope_data = OPE(left, right)
                correction = 0
                
                # 遍历极点，计算 rac{1}{n!} \partial^n \{BA\}_n
                # 注意：OPEData 中的极点系数就是 \{BA\}_n
                for n, coeff in ope_data.poles.items():
                    if n >= 1:
                        term = derivative_operator(coeff, n) / sp.factorial(n)
                        correction += term
                        
                # 递归化简修正项
                correction = simplify(correction)
                
                # 返回交换后的结果
                return parity_sign * NO(right, left) + correction
                
            except Exception as e:
                # 如果无法计算 OPE（例如未定义），则不交换，或者发出警告
                # 目前保持原样
                pass

    # 创建简化的 NO
    return NO(left, right)


def canonicalize(expr: Any) -> Any:
    """
    将表达式规范化

    比 simplify 更激进的化简，会：
    1. 完全展开所有嵌套的 NO
    2. 强制按字母顺序排列算符
    3. 合并所有同类项

    Args:
        expr: 要规范化的表达式

    Returns:
        规范化后的表达式

    Note:
        此函数需要完整的 OPE 信息才能正确展开
    """
    # TODO: 实现完整的规范化
    # 需要：
    # 1. NOCommuteHelp 来处理算符交换
    # 2. NO 展开规则
    # 3. 项的收集和合并
    return simplify(expr)


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
            key = ('other', str(op))
            terms[key] = terms.get(key, 0) + coeff
    elif isinstance(expr, NormalOrderedOperator):
        key = _make_no_key(expr)
        terms[key] = terms.get(key, 0) + 1
    else:
        # 其他项
        key = ('other', str(expr))
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
    return ('NO', left_key, right_key)


def _operator_to_key(op: Any) -> Tuple:
    """
    将算符转换为可哈希的键

    Args:
        op: 算符

    Returns:
        元组键
    """
    if isinstance(op, BasisOperator):
        return ('Basis', op.name, op.is_bosonic)
    elif isinstance(op, DerivativeOperator):
        return ('Deriv', _operator_to_key(op.base), op.order)
    elif isinstance(op, NormalOrderedOperator):
        return _make_no_key(op)
    else:
        return ('Other', str(op))
