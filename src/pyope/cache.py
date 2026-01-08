"""
性能优化：缓存模块

本模块提供 OPE 计算的缓存机制，以避免重复计算。
"""

from functools import lru_cache
from typing import Any, Tuple, Hashable
import sympy as sp
from sympy import Add, Mul


def make_operator_key(expr: Any) -> Hashable:
    """
    为算符表达式创建可哈希的键

    将算符表达式转换为可用于缓存键的元组形式。

    Args:
        expr: 算符或表达式

    Returns:
        可哈希的元组键

    Examples:
        >>> from pyope import BasisOperator
        >>> T = BasisOperator("T", bosonic=True)
        >>> key = make_operator_key(T)
        >>> isinstance(key, tuple)
        True
    """
    from .operators import (
        BasisOperator,
        DerivativeOperator,
        NormalOrderedOperator
    )
    from .constants import ConstantOperator

    # 处理 None 和零
    if expr is None or expr == 0:
        return ('zero',)

    # BasisOperator
    if isinstance(expr, BasisOperator):
        return ('basis', expr.name, expr.is_bosonic, expr._conformal_weight)

    # DerivativeOperator
    if isinstance(expr, DerivativeOperator):
        base_key = make_operator_key(expr.base)
        return ('deriv', base_key, expr.order)

    # NormalOrderedOperator
    if isinstance(expr, NormalOrderedOperator):
        left_key = make_operator_key(expr.left)
        right_key = make_operator_key(expr.right)
        return ('no', left_key, right_key)

    # ConstantOperator
    if isinstance(expr, ConstantOperator):
        return ('const', expr.name)

    # Sympy Add (加法)
    if isinstance(expr, Add):
        # 对项进行排序以确保相同表达式有相同的键
        term_keys = tuple(sorted(
            [make_operator_key(arg) for arg in expr.args],
            key=str
        ))
        return ('add', term_keys)

    # Sympy Mul (乘法)
    if isinstance(expr, Mul):
        # 提取标量和算符
        from .local_operator import extract_scalar_operator
        coeff, op = extract_scalar_operator(expr)

        # 如果是纯标量
        if coeff == expr:
            # 尝试将 sympy 表达式转换为字符串
            return ('scalar', str(expr))

        # 标量 * 算符
        op_key = make_operator_key(op)
        return ('mul', str(coeff), op_key)

    # Sympy Symbol 或其他表达式
    if isinstance(expr, sp.Expr):
        return ('expr', str(expr))

    # 其他类型：尝试使用字符串表示
    try:
        return ('other', str(expr), type(expr).__name__)
    except:
        # 最后的兜底：使用 id
        return ('id', id(expr))


def make_ope_cache_key(left: Any, right: Any) -> Tuple[Hashable, Hashable]:
    """
    为 OPE 计算创建缓存键

    Args:
        left: 左侧算符
        right: 右侧算符

    Returns:
        (left_key, right_key) 元组
    """
    left_key = make_operator_key(left)
    right_key = make_operator_key(right)
    return (left_key, right_key)


class OPECache:
    """
    OPE 计算结果的缓存

    使用字典存储已计算的 OPE 结果，避免重复计算。
    """

    def __init__(self, maxsize: int = 1024):
        """
        初始化缓存

        Args:
            maxsize: 最大缓存条目数
        """
        self.maxsize = maxsize
        self._cache = {}
        self._access_count = {}
        self.hits = 0
        self.misses = 0

    def get(self, left: Any, right: Any):
        """
        从缓存获取 OPE 结果

        Args:
            left: 左侧算符
            right: 右侧算符

        Returns:
            缓存的 OPEData 或 None（未命中）
        """
        try:
            key = make_ope_cache_key(left, right)
            if key in self._cache:
                self.hits += 1
                self._access_count[key] = self._access_count.get(key, 0) + 1
                return self._cache[key]
            else:
                self.misses += 1
                return None
        except Exception:
            # 如果无法创建键，跳过缓存
            self.misses += 1
            return None

    def put(self, left: Any, right: Any, result):
        """
        将 OPE 结果放入缓存

        Args:
            left: 左侧算符
            right: 右侧算符
            result: OPEData 结果
        """
        try:
            # 检查缓存大小
            if len(self._cache) >= self.maxsize:
                self._evict_lru()

            key = make_ope_cache_key(left, right)
            self._cache[key] = result
            self._access_count[key] = 1
        except Exception:
            # 如果无法创建键，跳过缓存
            pass

    def _evict_lru(self):
        """移除最少使用的缓存项"""
        if not self._cache:
            return

        # 找到访问次数最少的键
        lru_key = min(self._access_count.items(), key=lambda x: x[1])[0]
        del self._cache[lru_key]
        del self._access_count[lru_key]

    def clear(self):
        """清空缓存"""
        self._cache.clear()
        self._access_count.clear()
        self.hits = 0
        self.misses = 0

    def stats(self):
        """
        返回缓存统计信息

        Returns:
            包含命中率等信息的字典
        """
        total = self.hits + self.misses
        hit_rate = self.hits / total if total > 0 else 0

        return {
            'hits': self.hits,
            'misses': self.misses,
            'total_requests': total,
            'hit_rate': hit_rate,
            'cache_size': len(self._cache),
            'max_size': self.maxsize
        }


# 全局缓存实例
_global_ope_cache = OPECache(maxsize=2048)


def get_ope_cache() -> OPECache:
    """获取全局 OPE 缓存实例"""
    return _global_ope_cache


# 数值计算缓存

@lru_cache(maxsize=512)
def cached_pochhammer(q: int, n: int) -> int:
    """
    缓存的 Pochhammer 符号计算

    计算 (q-1)_n = (q-1)(q-2)...(q-n)

    Args:
        q: 参数 q
        n: 阶数 n

    Returns:
        Pochhammer 符号的值
    """
    result = 1
    for i in range(n):
        result *= (q - 1 - i)
    return result


@lru_cache(maxsize=512)
def cached_binomial(n: int, k: int):
    """
    缓存的二项式系数计算

    使用 sympy 的 binomial 但添加缓存。

    Args:
        n: 上标
        k: 下标

    Returns:
        C(n, k)
    """
    from sympy import binomial
    return binomial(n, k)


@lru_cache(maxsize=512)
def cached_factorial(n: int):
    """
    缓存的阶乘计算

    Args:
        n: 整数

    Returns:
        n!
    """
    from sympy import factorial
    return factorial(n)
