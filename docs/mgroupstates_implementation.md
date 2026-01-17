# mgroupstates 实现说明

## 概述

本文档说明 pyope 库中 `QuantumNumberGrouper` 类如何对应 Mathematica `null_states.wls` 中的分组逻辑。

## Mathematica 分组流程

在 `null_states.wls` 中，分组流程为：

```mathematica
(* 第 1 步：计算所有态的 m 量子数 *)
mtset = mtotal[#] & /@ fulllist;

(* 第 2 步：按 m 分组 - mgroupstates *)
mgroupstates = Table[
  Append[{}, finalgen[[#]]]& /@ Position[mtset, i] // Flatten,
  {i, Min[mtset], Max[mtset], 1/2}
] /. {} -> Nothing

(* 第 3 步：在 m 扇区内计算 r 量子数 *)
rtset = Table[
  Table[rtotal[mgfulllist[[i,j]]], {j,Length[mgfulllist[[i]]]}],
  {i,Length[mgfulllist]}
]

(* 第 4 步：在 m 扇区内按 r 分组 - fullgrouping *)
fullgrouping = Table[
  Table[Append[{},mgroupstates[[j,#]]]& /@ Position[rtset[[j]],i] // Flatten,
    {i,Min[rtset[[j]]],Max[rtset[[j]]],1/2}],
  {j,Length[rtset]}
] /. {} -> Nothing
```

### 数据结构

- `mgroupstates`: `List[ List[state] ]` - 二维嵌套列表
  - 外层索引：m 扇区（按 m 值隐式排序）
  - 内层列表：该 m 扇区内的态

- `fullgrouping`: `List[ List[ List[state] ] ]` - 三维嵌套列表
  - 第 1 层：m 扇区
  - 第 2 层：该 m 扇区内的 r 子扇区
  - 第 3 层：该 (m, r) 扇区的态列表

## Python 实现

### 方法对应关系

| Mathematica 变量/操作 | Python 方法 | 返回类型 |
|---------------------|------------|---------|
| `mgroupstates` | `group_by_m_only()` | `Dict[Fraction, List[Operator]]` |
| `fullgrouping` | `group_by_r_within_m()` | `Dict[Tuple[Fraction, Fraction], List[Operator]]` |
| `mparam = Union[mtset]` | `get_m_values()` | `List[Fraction]` |
| `sfullgrouping` | `group_by_m_only(..., only_non_negative_m=True)` | 同上，但只包含 m≥0 |

### 使用示例

```python
from pyope import QuantumNumberCalculator, QuantumNumberGrouper
from fractions import Fraction

# 定义量子数映射
quantum_number_map = {
    w: (Fraction(3, 2), Fraction(0)),  # (m, r)
    g: (Fraction(-1, 2), Fraction(1, 2)),
    # ...
}

# 创建分组器
calculator = QuantumNumberCalculator(quantum_number_map)
grouper = QuantumNumberGrouper(calculator)

# 方法 1: 只按 m 分组（对应 mgroupstates）
states_by_m = grouper.group_by_m_only(operators)
# 返回: {Fraction(0): [op1, op2], Fraction(3, 2): [op3], ...}

# 方法 2: 在 m 扇区内按 r 分组（对应 fullgrouping）
states_by_m_r = grouper.group_by_r_within_m(states_by_m)
# 返回: {(Fraction(0), Fraction(0)): [op1], (Fraction(3, 2), Fraction(0)): [op2], ...}

# 方法 3: 直接按 (m, r) 分组（pyope 特有，更简洁）
states_by_m_r_direct = grouper.group_operators(operators)
# 返回结果与方法 2 相同
```

## 关键改进

### 1. 显式映射 vs 隐式索引

**Mathematica 问题**：
- `mgroupstates[[j]]` 对应的 m 值需要从 `mparam[[j]]` 或重新计算得知
- 平行数组依赖：`mgroupstates` 与 `mparam` 必须保持索引对齐
- 易错且难调试

**Python 解决方案**：
```python
# 字典键直接就是 m 值，无需额外数组
states_by_m = {
    Fraction(0): [op1, op2],
    Fraction(3, 2): [op3]
}
```

### 2. 类型提示和自文档化

```python
def group_by_m_only(
    self,
    operators: List[Any],
    only_non_negative_m: bool = False
) -> Dict[Fraction, List[Any]]:
    """
    只按 m 量子数分组（对应 Mathematica 的 mgroupstates）

    物理意义：m 对应 osp(2|2) 中 gl(1) 子代数的荷，是加性守恒量。
    利用守恒性将系数矩阵自然分块对角化，从而降低计算复杂度。
    """
```

### 3. 灵活的 API 设计

提供三种分组方式：
1. `group_by_m_only()` - 单层 m 分组
2. `group_by_r_within_m()` - 两步法：m → (m, r)
3. `group_operators()` - 直接法：一步得到 (m, r)

用户可根据需要选择最合适的方法。

## 物理意义

### 为什么要分组？

1. **守恒荷**：m 和 r 都是加性守恒量
   - m: osp(2|2) gl(1) 荷
   - r: 外自同构 gl(1) 荷

2. **分块对角化**：由于守恒性，系数矩阵天然是块对角的
   ```
   [A₁  0   0  ]
   [0   A₂  0  ]  ← 每个块对应一个 (m, r) 扇区
   [0   0   A₃ ]
   ```

3. **计算效率**：
   - 不分组：计算一个 N×N 矩阵的秩，复杂度 O(N³)
   - 分组后：计算 k 个小矩阵的秩，复杂度 Σ O(nᵢ³) << O(N³)

### 为什么先按 m，再按 r？

1. **Mathematica 实现考虑**：
   - 利用 m 的对称性优化：只计算 m≥0 扇区（`sfullgrouping`）
   - 在字符求和时用 `z^m + z^{-m}` 补回负 m 贡献

2. **Python 实现**：
   - 可以直接按 (m, r) 分组（`group_operators()`）
   - 也支持两步法以保持与 Mathematica 逻辑的对应

## 测试验证

运行 `tests/test_mgroupstates_demo.py` 验证：

```bash
cd tests
python test_mgroupstates_demo.py
```

验证内容：
- ✓ `group_by_m_only()` 正确按 m 分组
- ✓ `group_by_r_within_m()` 正确在 m 扇区内按 r 分组
- ✓ 两步法与直接法结果一致
- ✓ `only_non_negative_m=True` 正确过滤负 m 扇区

## 相关文件

- 实现：`src/pyope/null_states.py:648-812`
- 测试：`tests/test_mgroupstates_demo.py`
- 参考：`.claude/skills/voa/computations/null_states.wls:454-470`
