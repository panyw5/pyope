"""
性能基准测试脚本

测试 pyope 库在不同场景下的性能，特别是缓存优化的效果。

测试场景：
1. 基础 OPE 计算（Virasoro 代数）
2. 复合 OPE 计算（正规序算符）
3. Jacobi 恒等式验证
4. 缓存命中率和加速比分析
"""

import time
import sys
from typing import Callable, Dict, Any, List, Tuple
import sympy as sp

# 添加 src 到路径
sys.path.insert(0, '/Users/lelouch/pyope/src')

from pyope import (
    BasisOperator,
    OPE,
    NO,
    d,
    One,
    get_ope_cache,
    check_jacobi_identity,
    Bosonic,
    Fermionic,
)
from pyope.ope_data import OPEData


class PerformanceBenchmark:
    """性能基准测试类"""

    def __init__(self):
        self.results = {}

    def time_function(self, func: Callable, name: str, iterations: int = 1) -> Dict[str, Any]:
        """
        测量函数执行时间

        Args:
            func: 要测试的函数
            name: 测试名称
            iterations: 迭代次数

        Returns:
            包含执行时间和统计信息的字典
        """
        # 预热
        try:
            func()
        except Exception as e:
            print(f"预热失败: {e}")

        # 清除缓存以获得公平比较
        cache = get_ope_cache()
        cache.clear()

        # 开始计时
        start_time = time.time()

        for _ in range(iterations):
            try:
                func()
            except Exception as e:
                print(f"测试 {name} 失败: {e}")
                return {
                    'name': name,
                    'error': str(e),
                    'time': 0,
                    'iterations': iterations,
                }

        end_time = time.time()
        total_time = end_time - start_time
        avg_time = total_time / iterations

        # 获取缓存统计
        cache_stats = cache.stats()

        result = {
            'name': name,
            'total_time': total_time,
            'avg_time': avg_time,
            'iterations': iterations,
            'cache_stats': cache_stats,
        }

        self.results[name] = result
        return result

    def compare_with_without_cache(
        self,
        func: Callable,
        name: str,
        iterations: int = 10
    ) -> Dict[str, Any]:
        """
        比较带缓存和不带缓存的性能

        Args:
            func: 要测试的函数（会多次调用以触发缓存）
            name: 测试名称
            iterations: 迭代次数

        Returns:
            比较结果字典
        """
        # 测试不带缓存（每次都清除缓存）
        cache = get_ope_cache()

        start_time = time.time()
        for _ in range(iterations):
            cache.clear()
            try:
                func()
            except Exception as e:
                print(f"无缓存测试 {name} 失败: {e}")
        end_time = time.time()

        time_without_cache = (end_time - start_time) / iterations

        # 测试带缓存
        cache.clear()

        start_time = time.time()
        for _ in range(iterations):
            try:
                func()
            except Exception as e:
                print(f"带缓存测试 {name} 失败: {e}")
        end_time = time.time()

        time_with_cache = (end_time - start_time) / iterations
        cache_stats = cache.stats()

        # 计算加速比
        speedup = time_without_cache / time_with_cache if time_with_cache > 0 else 0

        result = {
            'name': name,
            'time_without_cache': time_without_cache,
            'time_with_cache': time_with_cache,
            'speedup': speedup,
            'cache_stats': cache_stats,
            'iterations': iterations,
        }

        self.results[f"{name}_comparison"] = result
        return result

    def print_results(self):
        """打印所有测试结果"""
        print("\n" + "="*80)
        print("性能基准测试结果")
        print("="*80)

        for name, result in self.results.items():
            print(f"\n【{name}】")

            if 'error' in result:
                print(f"  ❌ 错误: {result['error']}")
                continue

            if 'speedup' in result:
                # 比较测试结果
                print(f"  迭代次数: {result['iterations']}")
                print(f"  无缓存平均时间: {result['time_without_cache']*1000:.2f} ms")
                print(f"  有缓存平均时间: {result['time_with_cache']*1000:.2f} ms")
                print(f"  ⚡ 加速比: {result['speedup']:.2f}x")

                if 'cache_stats' in result:
                    stats = result['cache_stats']
                    print(f"  缓存命中率: {stats['hit_rate']*100:.1f}%")
                    print(f"  缓存命中: {stats['hits']}, 未命中: {stats['misses']}")
                    print(f"  缓存大小: {stats['cache_size']}/{stats['max_size']}")
            else:
                # 单次测试结果
                print(f"  迭代次数: {result['iterations']}")
                print(f"  总时间: {result['total_time']*1000:.2f} ms")
                print(f"  平均时间: {result['avg_time']*1000:.2f} ms")

                if 'cache_stats' in result:
                    stats = result['cache_stats']
                    print(f"  缓存命中率: {stats['hit_rate']*100:.1f}%")
                    print(f"  缓存命中: {stats['hits']}, 未命中: {stats['misses']}")


def setup_virasoro():
    """设置 Virasoro 代数"""
    c = sp.Symbol('c')
    T = BasisOperator("T", bosonic=True, conformal_weight=2)

    # 定义 Virasoro OPE: T(z)T(w) ~ c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
    OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])

    return T, c


def setup_superconformal():
    """设置超共形代数 (N=1)"""
    c = sp.Symbol('c')
    T = BasisOperator("T", bosonic=True, conformal_weight=2)
    G = BasisOperator("G", bosonic=False, conformal_weight=sp.Rational(3, 2))

    # Virasoro OPE
    OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])

    # TG OPE: T(z)G(w) ~ (3/2)G/(z-w)^2 + ∂G/(z-w)
    OPE[T, G] = OPE.make([sp.Rational(3,2)*G, d(G)])

    # GG OPE: G(z)G(w) ~ (2/3)c/(z-w)^3 + 2T/(z-w)
    OPE[G, G] = OPE.make([sp.Rational(2,3)*c*One, 0, 2*T])

    return T, G, c


# 测试场景定义

def benchmark_basic_ope():
    """基准测试 1: 基础 OPE 计算"""
    print("\n" + "="*80)
    print("基准测试 1: 基础 OPE 计算 (Virasoro)")
    print("="*80)

    T, c = setup_virasoro()
    benchmark = PerformanceBenchmark()

    # 测试 1: 简单 OPE 计算
    def test_simple_ope():
        return OPE(T, T)

    benchmark.time_function(test_simple_ope, "T×T OPE", iterations=100)

    # 测试 2: 导数 OPE
    def test_derivative_ope():
        dT = d(T)
        return OPE(dT, T)

    benchmark.time_function(test_derivative_ope, "∂T×T OPE", iterations=100)

    # 测试 3: 二阶导数 OPE
    def test_second_derivative_ope():
        d2T = d(T, 2)
        return OPE(d2T, T)

    benchmark.time_function(test_second_derivative_ope, "∂²T×T OPE", iterations=100)

    benchmark.print_results()
    return benchmark


def benchmark_composite_ope():
    """基准测试 2: 复合 OPE 计算"""
    print("\n" + "="*80)
    print("基准测试 2: 复合 OPE 计算")
    print("="*80)

    T, c = setup_virasoro()
    benchmark = PerformanceBenchmark()

    # 测试 1: NO(T, T) × T
    def test_no_tt_t():
        NO_TT = NO(T, T)
        return OPE(NO_TT, T)

    benchmark.compare_with_without_cache(test_no_tt_t, "NO(T,T)×T", iterations=20)

    # 测试 2: T × NO(T, T)
    def test_t_no_tt():
        NO_TT = NO(T, T)
        return OPE(T, NO_TT)

    benchmark.compare_with_without_cache(test_t_no_tt, "T×NO(T,T)", iterations=20)

    # 测试 3: NO(T, ∂T) × T
    def test_no_t_dt_t():
        dT = d(T)
        NO_TdT = NO(T, dT)
        return OPE(NO_TdT, T)

    benchmark.compare_with_without_cache(test_no_t_dt_t, "NO(T,∂T)×T", iterations=20)

    benchmark.print_results()
    return benchmark


def benchmark_jacobi_identity():
    """基准测试 3: Jacobi 恒等式验证"""
    print("\n" + "="*80)
    print("基准测试 3: Jacobi 恒等式验证")
    print("="*80)

    T, c = setup_virasoro()
    benchmark = PerformanceBenchmark()

    # 测试 Virasoro Jacobi 恒等式
    def test_virasoro_jacobi():
        return check_jacobi_identity(T, T, T, simplify_func=sp.expand)

    benchmark.compare_with_without_cache(
        test_virasoro_jacobi,
        "Virasoro Jacobi(T,T,T)",
        iterations=5
    )

    benchmark.print_results()
    return benchmark


def benchmark_superconformal():
    """基准测试 4: 超共形代数"""
    print("\n" + "="*80)
    print("基准测试 4: 超共形代数 (N=1)")
    print("="*80)

    T, G, c = setup_superconformal()
    benchmark = PerformanceBenchmark()

    # 测试 1: T×G OPE
    def test_tg_ope():
        return OPE(T, G)

    benchmark.time_function(test_tg_ope, "T×G OPE", iterations=50)

    # 测试 2: G×G OPE
    def test_gg_ope():
        return OPE(G, G)

    benchmark.time_function(test_gg_ope, "G×G OPE", iterations=50)

    # 测试 3: NO(T,G) × G
    def test_no_tg_g():
        NO_TG = NO(T, G)
        return OPE(NO_TG, G)

    benchmark.compare_with_without_cache(test_no_tg_g, "NO(T,G)×G", iterations=10)

    # 测试 4: Jacobi 恒等式 (T, G, G)
    def test_superconformal_jacobi():
        return check_jacobi_identity(T, G, G, simplify_func=sp.expand)

    benchmark.compare_with_without_cache(
        test_superconformal_jacobi,
        "Superconformal Jacobi(T,G,G)",
        iterations=3
    )

    benchmark.print_results()
    return benchmark


def benchmark_cache_effectiveness():
    """基准测试 5: 缓存效果分析"""
    print("\n" + "="*80)
    print("基准测试 5: 缓存效果分析")
    print("="*80)

    T, c = setup_virasoro()
    cache = get_ope_cache()

    print("\n【测试场景】: 重复计算相同的 OPE")
    print("  预期: 第一次计算后，后续计算应该从缓存获取，速度显著提升")

    # 清除缓存
    cache.clear()

    # 第一次计算（缓存未命中）
    start = time.time()
    for _ in range(100):
        OPE(T, T)
    first_run = time.time() - start

    stats_after_first = cache.stats()
    print(f"\n  第一轮 (100次): {first_run*1000:.2f} ms")
    print(f"    缓存命中率: {stats_after_first['hit_rate']*100:.1f}%")

    # 第二次计算（应该全部命中缓存）
    start = time.time()
    for _ in range(100):
        OPE(T, T)
    second_run = time.time() - start

    stats_after_second = cache.stats()
    print(f"\n  第二轮 (100次): {second_run*1000:.2f} ms")
    print(f"    缓存命中率: {stats_after_second['hit_rate']*100:.1f}%")
    print(f"    加速比: {first_run/second_run:.2f}x")

    # 测试不同复杂度的计算
    print("\n【测试场景】: 不同复杂度计算的缓存效果")

    test_cases = [
        ("T×T", lambda: OPE(T, T)),
        ("∂T×T", lambda: OPE(d(T), T)),
        ("∂²T×T", lambda: OPE(d(T, 2), T)),
        ("NO(T,T)×T", lambda: OPE(NO(T, T), T)),
        ("T×NO(T,T)", lambda: OPE(T, NO(T, T))),
    ]

    cache.clear()

    for name, func in test_cases:
        # 第一次（无缓存）
        start = time.time()
        func()
        time_first = time.time() - start

        # 第二次（有缓存）
        start = time.time()
        func()
        time_second = time.time() - start

        speedup = time_first / time_second if time_second > 0 else float('inf')
        print(f"\n  {name}:")
        print(f"    首次: {time_first*1000:.3f} ms")
        print(f"    缓存: {time_second*1000:.3f} ms")
        print(f"    加速: {speedup:.1f}x")


def main():
    """运行所有基准测试"""
    print("="*80)
    print("PyOPE 性能基准测试套件")
    print("="*80)
    print("\n本测试套件将评估:")
    print("  1. 基础 OPE 计算性能")
    print("  2. 复合 OPE 计算性能")
    print("  3. Jacobi 恒等式验证性能")
    print("  4. 超共形代数计算性能")
    print("  5. 缓存机制有效性")
    print("\n开始测试...\n")

    # 运行所有基准测试
    try:
        benchmark_basic_ope()
    except Exception as e:
        print(f"\n❌ 基础 OPE 测试失败: {e}")

    try:
        benchmark_composite_ope()
    except Exception as e:
        print(f"\n❌ 复合 OPE 测试失败: {e}")

    try:
        benchmark_jacobi_identity()
    except Exception as e:
        print(f"\n❌ Jacobi 恒等式测试失败: {e}")

    try:
        benchmark_superconformal()
    except Exception as e:
        print(f"\n❌ 超共形代数测试失败: {e}")

    try:
        benchmark_cache_effectiveness()
    except Exception as e:
        print(f"\n❌ 缓存效果测试失败: {e}")

    # 最终总结
    print("\n" + "="*80)
    print("测试完成")
    print("="*80)

    # 打印最终缓存统计
    cache = get_ope_cache()
    final_stats = cache.stats()

    print("\n【最终缓存统计】")
    print(f"  总请求: {final_stats['total_requests']}")
    print(f"  缓存命中: {final_stats['hits']}")
    print(f"  缓存未命中: {final_stats['misses']}")
    print(f"  命中率: {final_stats['hit_rate']*100:.1f}%")
    print(f"  缓存使用: {final_stats['cache_size']}/{final_stats['max_size']}")
    print()


if __name__ == "__main__":
    main()
