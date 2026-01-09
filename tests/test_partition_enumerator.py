"""
测试基于整数分拆的算符枚举器

验证新的枚举器能否正确生成任意数量生成元的乘积
"""

from fractions import Fraction
from pyope.null_states import OperatorEnumerator, integer_partitions
from pyope.operators import BasisOperator


def test_integer_partitions():
    """测试整数分拆生成器"""
    # 测试整数
    partitions_3 = integer_partitions(3)
    print(f"\n整数 3 的分拆: {partitions_3}")
    assert [3] in partitions_3
    assert [2, 1] in partitions_3
    assert [1, 1, 1] in partitions_3

    # 测试半整数
    partitions_2_5 = integer_partitions(Fraction(5, 2))
    print(f"\n半整数 5/2 的分拆: {partitions_2_5}")
    assert [Fraction(5, 2)] in partitions_2_5
    assert [2, Fraction(1, 2)] in partitions_2_5
    assert [Fraction(3, 2), 1] in partitions_2_5

    print("\n✓ 整数分拆测试通过")


def test_simple_enumerator():
    """测试简单的算符枚举"""
    # 创建简单的生成元
    J = BasisOperator('J', 0, 1)  # 权重为 1

    generators = {
        'J': {'op': J, 'weight': Fraction(1)}
    }

    enumerator = OperatorEnumerator(generators)

    # 测试 level=2
    ops_level2 = enumerator.enumerate_operators(Fraction(2))
    print(f"\nLevel 2 算符数量: {len(ops_level2)}")
    print("算符列表:")
    for op in ops_level2:
        print(f"  {op}")

    # 应该包含：
    # - ∂J (单个生成元的导数)
    # - NO(J, J) (两个生成元的乘积)
    assert len(ops_level2) >= 2

    # 测试 level=3
    ops_level3 = enumerator.enumerate_operators(Fraction(3))
    print(f"\nLevel 3 算符数量: {len(ops_level3)}")
    print("算符列表:")
    for op in ops_level3:
        print(f"  {op}")

    # 应该包含：
    # - ∂²J
    # - NO(∂J, J)
    # - NO(J, ∂J)
    # - NO(J, J, J) (三个生成元的乘积！)
    assert len(ops_level3) >= 4

    # 测试 level=4
    ops_level4 = enumerator.enumerate_operators(Fraction(4))
    print(f"\nLevel 4 算符数量: {len(ops_level4)}")
    print("算符列表:")
    for op in ops_level4:
        print(f"  {op}")

    # 应该包含：
    # - ∂³J
    # - NO(∂²J, J), NO(∂J, ∂J), NO(J, ∂²J)
    # - NO(∂J, J, J), NO(J, ∂J, J), NO(J, J, ∂J)
    # - NO(J, J, J, J) (四个生成元的乘积！)
    assert len(ops_level4) >= 8

    print("\n✓ 简单枚举器测试通过")


def test_z3_w_algebra_level4():
    """测试 Z₃ W-algebra 在 level=4 的枚举"""
    from pyope.operators import BasisOperator

    # 定义 Z₃ W-algebra 的生成元
    # 根据 null_states.txt 的定义
    W = BasisOperator('W', 0, Fraction(3, 2))      # 权重 3/2
    J = BasisOperator('J', 0, 1)                    # 权重 1
    Wb = BasisOperator('Wb', 0, Fraction(3, 2))    # 权重 3/2
    T = BasisOperator('T', 0, 2)                    # 权重 2
    G = BasisOperator('G', 0, Fraction(3, 2))      # 权重 3/2
    Gt = BasisOperator('Gt', 0, Fraction(3, 2))    # 权重 3/2
    Gw = BasisOperator('Gw', 0, 2)                  # 权重 2
    Gwb = BasisOperator('Gwb', 0, 2)                # 权重 2

    generators = {
        'W': {'op': W, 'weight': Fraction(3, 2)},
        'J': {'op': J, 'weight': Fraction(1)},
        'Wb': {'op': Wb, 'weight': Fraction(3, 2)},
        'T': {'op': T, 'weight': Fraction(2)},
        'G': {'op': G, 'weight': Fraction(3, 2)},
        'Gt': {'op': Gt, 'weight': Fraction(3, 2)},
        'Gw': {'op': Gw, 'weight': Fraction(2)},
        'Gwb': {'op': Gwb, 'weight': Fraction(2)},
    }

    enumerator = OperatorEnumerator(generators)

    # 枚举 level=4 的算符
    ops = enumerator.enumerate_operators(Fraction(4))

    print(f"\n=== Z₃ W-algebra Level 4 枚举结果 ===")
    print(f"总算符数量: {len(ops)}")
    print(f"\n算符列表:")
    for i, op in enumerate(ops, 1):
        print(f"{i:3d}. {op}")

    # 根据 gemini 的分析，Mathematica 在 level=4 生成了 61 个算符
    # 但经过去重后应该是更少的数量
    # 让我们检查是否包含关键的多生成元乘积

    ops_str = [str(op) for op in ops]

    # 检查是否包含四个生成元的乘积（例如 J^4）
    has_four_generators = any('NO' in s and s.count(',') >= 3 for s in ops_str)
    print(f"\n包含四生成元乘积: {has_four_generators}")

    # 检查是否包含三个生成元的乘积
    has_three_generators = any('NO' in s and s.count(',') >= 2 for s in ops_str)
    print(f"包含三生成元乘积: {has_three_generators}")

    assert has_three_generators, "应该包含三生成元乘积"
    assert has_four_generators, "应该包含四生成元乘积"

    print("\n✓ Z₃ W-algebra Level 4 枚举测试通过")
    print(f"\n注意：Mathematica 参考结果是 61 个算符（去重前）")
    print(f"当前实现生成了 {len(ops)} 个算符")


if __name__ == '__main__':
    test_integer_partitions()
    test_simple_enumerator()
    test_z3_w_algebra_level4()
    print("\n" + "="*50)
    print("所有测试通过！")
    print("="*50)
