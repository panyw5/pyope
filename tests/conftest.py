"""
Pytest configuration for pyope tests

This file contains shared fixtures and configuration for all tests.
"""

from pathlib import Path

import pytest
import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    Fermionic,
    OPE,
    NO,
    MakeOPE,
    d,
    dn,
    One,
    Zero,
    simplify,
)
from pyope.cache import get_ope_cache
from pyope.registry import ope_registry


def pytest_sessionstart(session):
    """Fail fast if tests import pyope from outside this worktree."""
    del session

    import pyope

    repo_root = Path(__file__).resolve().parents[1]
    expected_package_dir = (repo_root / "src" / "pyope").resolve()
    actual_package_dir = Path(pyope.__file__).resolve().parent

    if actual_package_dir != expected_package_dir:
        raise RuntimeError(
            "pytest imported pyope from the wrong source tree: "
            f"expected {expected_package_dir}, got {actual_package_dir}. "
            'Run tests with PYTHONPATH="$PWD/src" python -m pytest ... '
            "or use a worktree-local editable install."
        )


@pytest.fixture(autouse=True, scope="function")
def disable_cache_for_tests():
    """
    Disable OPE cache for all tests to ensure test isolation.

    This fixture automatically runs before each test function and:
    1. Disables the global OPE cache
    2. Clears the registry
    3. Re-enables the cache after the test

    This ensures that tests don't interfere with each other through
    cached results.
    """
    cache = get_ope_cache()

    # Disable cache and clear registry before test
    cache.disable()
    cache.clear()
    ope_registry.clear()

    yield

    # Re-enable cache after test (for performance tests if needed)
    cache.enable()
    cache.clear()
    ope_registry.clear()


@pytest.fixture
def sl2_algebra():
    """
    SL(2) Kac-Moody 代数 fixture

    定义 sl(2) 的三个生成元 J+, J0, J- 及其 OPE 规则。

    Returns:
        dict: 包含以下键的字典
            - Jplus: J+ 算符
            - Jzero: J0 算符
            - Jminus: J- 算符
            - k: level 参数（sympy Symbol）
            - clear: 清理函数（测试后调用）

    OPE 规则（level k）：
        J0(z) J0(w) = k/2 / (z-w)^2
        J0(z) J+(w) = J+(w) / (z-w)
        J0(z) J-(w) = -J-(w) / (z-w)
        J+(z) J-(w) = k / (z-w)^2 + 2*J0(w) / (z-w)

    Examples:
        >>> def test_example(sl2_algebra):
        ...     Jp = sl2_algebra['Jplus']
        ...     J0 = sl2_algebra['Jzero']
        ...     result = simplify(NO(Jp, J0))
        ...     sl2_algebra['clear']()  # 清理 OPE 定义
    """
    # 定义算符
    Jplus = BasisOperator("Jplus")
    Jzero = BasisOperator("Jzero")
    Jminus = BasisOperator("Jminus")

    # 声明为玻色算符
    Bosonic(Jplus, Jzero, Jminus)

    # 定义 level 参数
    k = sp.Symbol("k")

    # 定义 OPE 规则
    # J0 * J0 = k/2 / (z-w)^2
    OPE[Jzero, Jzero] = MakeOPE(
        [
            sp.Rational(1, 2) * k * One,  # (z-w)^{-2}
            Zero,  # (z-w)^{-1}
        ]
    )

    # J0 * J+ = J+ / (z-w)
    OPE[Jzero, Jplus] = MakeOPE(
        [
            Zero,  # (z-w)^{-2}
            Jplus,  # (z-w)^{-1}
        ]
    )

    # J0 * J- = -J- / (z-w)
    OPE[Jzero, Jminus] = MakeOPE(
        [
            Zero,  # (z-w)^{-2}
            -Jminus,  # (z-w)^{-1}
        ]
    )

    # J+ * J- = k / (z-w)^2 + 2*J0 / (z-w)
    OPE[Jplus, Jminus] = MakeOPE(
        [
            k * One,  # (z-w)^{-2}
            2 * Jzero,  # (z-w)^{-1}
        ]
    )

    def clear():
        """清理 OPE 定义"""
        # 清理这些算符的 OPE 定义
        keys_to_remove = []
        for key in ope_registry._opes.keys():
            if any(op in [Jplus, Jzero, Jminus] for op in key):
                keys_to_remove.append(key)
        for key in keys_to_remove:
            del ope_registry._opes[key]

    return {
        "Jplus": Jplus,
        "Jzero": Jzero,
        "Jminus": Jminus,
        "k": k,
        "clear": clear,
    }


@pytest.fixture
def w3_algebra():
    """
    W3 代数 fixture

    定义 W3 代数的两个生成元 T, W 及其 OPE 规则。

    Returns:
        dict: 包含以下键的字典
            - T: 能量-动量张量
            - W: W 算符（spin-3）
            - c: 中心荷参数（sympy Symbol）
            - beta: β 参数（sympy Symbol）
            - Lambda: 辅助算符 Λ = NO(T,T) - (3/10)*∂²T
            - clear: 清理函数（测试后调用）

    OPE 规则：
        T(z) T(w) = c/2 / (z-w)^4 + 2*T(w) / (z-w)^2 + ∂T(w) / (z-w)
        T(z) W(w) = 3*W(w) / (z-w)^2 + ∂W(w) / (z-w)
        W(z) W(w) = c / (z-w)^6 + 2*T(w) / (z-w)^4 + ∂T(w) / (z-w)^3
                    + (2*β*Λ(w) + (3/10)*∂²T(w)) / (z-w)^2
                    + (β*∂Λ(w) + (1/15)*∂³T(w)) / (z-w)

    其中 Λ(w) = NO(T,T)(w) - (3/10)*∂²T(w)

    Examples:
        >>> def test_example(w3_algebra):
        ...     T = w3_algebra['T']
        ...     W = w3_algebra['W']
        ...     result = OPE(T, T)
        ...     w3_algebra['clear']()  # 清理 OPE 定义
    """
    # 定义算符
    T = BasisOperator("T")
    W = BasisOperator("W")

    # 声明为玻色算符
    Bosonic(T, W)

    # 定义参数
    c = sp.Symbol("c")
    beta = sp.Symbol("beta")

    # 定义辅助算符 Λ = NO(T,T) - (3/10)*∂²T
    Lambda = NO(T, T) - sp.Rational(3, 10) * dn(2, T)

    # 定义 OPE 规则
    # T * T = c/2 / (z-w)^4 + 2*T / (z-w)^2 + ∂T / (z-w)
    OPE[T, T] = MakeOPE(
        [
            sp.Rational(1, 2) * c * One,  # (z-w)^{-4}
            Zero,  # (z-w)^{-3}
            2 * T,  # (z-w)^{-2}
            d(T),  # (z-w)^{-1}
        ]
    )

    # T * W = 3*W / (z-w)^2 + ∂W / (z-w)
    OPE[T, W] = MakeOPE(
        [
            Zero,  # (z-w)^{-3}
            Zero,  # (z-w)^{-2} 实际上应该是 3*W，但索引从最高极点开始
            3 * W,  # (z-w)^{-2}
            d(W),  # (z-w)^{-1}
        ]
    )

    # W * W = c / (z-w)^6 + 2*T / (z-w)^4 + ∂T / (z-w)^3
    #         + (2*β*Λ + (3/10)*∂²T) / (z-w)^2
    #         + (β*∂Λ + (1/15)*∂³T) / (z-w)
    OPE[W, W] = MakeOPE(
        [
            c * One,  # (z-w)^{-6}
            Zero,  # (z-w)^{-5}
            2 * T,  # (z-w)^{-4}
            d(T),  # (z-w)^{-3}
            2 * beta * Lambda + sp.Rational(3, 10) * dn(2, T),  # (z-w)^{-2}
            beta * d(Lambda) + sp.Rational(1, 15) * dn(3, T),  # (z-w)^{-1}
        ]
    )

    def clear():
        """清理 OPE 定义"""
        # 清理这些算符的 OPE 定义
        keys_to_remove = []
        for key in ope_registry._opes.keys():
            if any(op in [T, W] for op in key):
                keys_to_remove.append(key)
        for key in keys_to_remove:
            del ope_registry._opes[key]

    return {
        "T": T,
        "W": W,
        "c": c,
        "beta": beta,
        "Lambda": Lambda,
        "clear": clear,
    }


# 辅助函数：检查 API 能力
def supports_derivative():
    """检查是否支持导数操作"""
    try:
        from pyope import d, dn

        return True
    except ImportError:
        return False


def supports_no():
    """检查是否支持正规序操作"""
    try:
        from pyope import NO

        return True
    except ImportError:
        return False


def supports_simplify():
    """检查是否支持 simplify 操作"""
    try:
        from pyope import simplify

        return True
    except ImportError:
        return False


# pytest markers
def pytest_configure(config):
    """配置 pytest markers"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "mathematica_ref: marks tests that reference Mathematica results"
    )
    config.addinivalue_line(
        "markers", "requires_derivative: marks tests that require derivative support"
    )
    config.addinivalue_line(
        "markers", "requires_no: marks tests that require normal ordering support"
    )
