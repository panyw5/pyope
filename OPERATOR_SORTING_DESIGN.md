# 算符枚举排序方案设计

## 问题描述

在实现 W-algebra 的 null states 计算时，需要枚举给定 level 的所有独立算符组合。当前实现与 Mathematica 参考实现的排序不一致，导致可能的重复计数或遗漏。

## Mathematica 的排序规则

### 1. 基本算符的定义顺序

在 `null_states.wls` 中，算符按以下顺序定义：

```mathematica
opfields = {ww[0,1,z], jj0[0,1,z], wwb[0,1,z], tt[0,1,z],
            gg[0,1,z], ggt[0,1,z], ggw[0,1,z], ggwb[0,1,z]}
fields = {ww, jj0, wwb, tt, gg, ggt, ggw, ggwb}
```

对应 Python 中的：
- `w` (ww): 权重 3/2, m=3/2, r=0
- `j0` (jj0): 权重 1, m=0, r=0
- `wb` (wwb): 权重 3/2, m=-3/2, r=0
- `t` (tt): 权重 2, m=0, r=0
- `g` (gg): 权重 3/2, m=-1/2, r=1/2
- `gt` (ggt): 权重 3/2, m=1/2, r=-1/2
- `gw` (ggw): 权重 2, m=1, r=1/2
- `gwb` (ggwb): 权重 2, m=-1, r=-1/2

**关键点：这个顺序是固定的，用于确定同一权重下不同算符的优先级。**

### 2. 整数分拆的顺序

对于 level L，生成所有整数分拆，按**降序**排列：

```
level 4 的分拆：
[4]
[7/2, 1/2]
[3, 1]
[5/2, 3/2]
[5/2, 1, 1/2]
...
```

### 3. 分拆到算符的映射

对于分拆 `[w1, w2, ..., wn]`（降序）：

1. **权重大的在左（前），权重小的在右（后）**
2. 对于每个权重 `wi`，找到所有可以达到该权重的算符（包括导数）
3. 按照 `opfields` 的顺序枚举这些算符
4. 构造正规序乘积：`NO(op1, NO(op2, NO(...)))`

### 4. 算符在同一权重下的排序

当多个算符（包括导数）可以达到同一权重时：

**规则：按照基本算符在 `opfields` 中的顺序排序**

例如，权重 3/2 可以由以下算符达到：
- `w` (基础权重 3/2, 导数 0)
- `wb` (基础权重 3/2, 导数 0)
- `g` (基础权重 3/2, 导数 0)
- `gt` (基础权重 3/2, 导数 0)
- `∂j0` (基础权重 1, 导数 1)

排序为：`w, j0', wb, g, gt`（按 opfields 顺序，j0 在 w 和 wb 之间）

**注意：导数不改变基本算符的相对顺序，∂j0 的位置由 j0 在 opfields 中的位置决定。**

### 5. 示例分析

#### 示例 1：分拆 [5/2, 3/2]

权重 5/2 可以由：
- `∂w` (3/2 + 1 = 5/2)
- `∂wb` (3/2 + 1 = 5/2)
- `∂g` (3/2 + 1 = 5/2)
- `∂gt` (3/2 + 1 = 5/2)
- `∂²j0` (1 + 2 = 3)  ❌ 不是 5/2
- `∂t` (2 + 1/2 = 5/2)  ❌ 导数必须是整数

权重 3/2 可以由：
- `w, wb, g, gt`

按 opfields 顺序，生成的算符组合：
```
NO(∂w, w)    # NO[w[1,1,z], w[0,1,z]]
NO(∂w, wb)   # NO[w[1,1,z], wb[0,1,z]]
NO(∂w, g)    # NO[w[1,1,z], g[0,1,z]]  ✓ 你的例子
NO(∂w, gt)   # NO[w[1,1,z], gt[0,1,z]]
NO(∂wb, w)   # NO[wb[1,1,z], w[0,1,z]]
...
```

#### 示例 2：分拆 [3/2, 3/2]

两个位置都是权重 3/2，可以由 `w, wb, g, gt`：

```
NO(w, w)     # NO[w[0,1,z], w[0,1,z]]
NO(w, wb)    # NO[w[0,1,z], wb[0,1,z]]
NO(w, g)     # NO[w[0,1,z], g[0,1,z]]
NO(w, gt)    # NO[w[0,1,z], gt[0,1,z]]
NO(wb, w)    # NO[wb[0,1,z], w[0,1,z]]
NO(wb, wb)   # NO[wb[0,1,z], wb[0,1,z]]
...
NO(g, g)     # NO[g[0,1,z], g[0,1,z]]
NO(g, gt)    # NO[g[0,1,z], gt[0,1,z]]  ✓ 你的例子
...
```

**注意：这里会产生重复！** 例如 `NO(w, w)` 和 `NO(w, w)` 是同一个算符。

### 6. 去重策略

对于玻色算符，`NO(A, B) = NO(B, A)`（交换对称）
对于费米算符，`NO(A, B) = -NO(B, A)`（反交换）

**Mathematica 的策略：**
- 使用 `itertools.product` 生成所有组合
- 使用字符串表示去重：`seen.add(str(op))`
- 依赖 OPE 简化规则自动处理对称性

**推荐策略：**
1. **生成所有组合**（不预先过滤）
2. **使用字符串表示去重**
3. **依赖简化系统处理对称性**

原因：
- 预先判断对称性很复杂（需要考虑嵌套 NO、混合玻色/费米等）
- 简化系统已经实现了完整的对称性处理
- 去重可以有效移除重复项

## Python 实现方案

### 方案 1：显式排序键（推荐）

为每个生成元分配一个排序键（order key），基于 `opfields` 的顺序：

```python
class OperatorEnumerator:
    def __init__(self, generators: Dict[str, Dict[str, Any]]):
        self.generators = generators

        # 定义算符的排序顺序（对应 Mathematica 的 opfields）
        self.operator_order = ['w', 'j0', 'wb', 't', 'g', 'gt', 'gw', 'gwb']

        # 创建排序键映射
        self.sort_key_map = {name: i for i, name in enumerate(self.operator_order)}

    def _get_sort_key(self, gen_name: str) -> int:
        """获取生成元的排序键"""
        return self.sort_key_map.get(gen_name, 999)  # 未知算符放最后

    def _generate_operators_at_weight(
        self,
        weight: Fraction,
        max_derivative_order: int
    ) -> List[Tuple[str, Any, int]]:
        """
        生成指定权重的所有可能算符

        Returns:
            List of (gen_name, operator, sort_key)
        """
        operators = []

        for name, gen_info in self.generators.items():
            base_weight = gen_info['weight']
            deriv_order = weight - base_weight

            if deriv_order >= 0 and deriv_order == int(deriv_order):
                deriv_order_int = int(deriv_order)
                if deriv_order_int <= max_derivative_order:
                    # 构造算符
                    if deriv_order_int == 0:
                        op = gen_info['op']
                    else:
                        op = d(gen_info['op'], deriv_order_int)

                    # 添加排序键
                    sort_key = self._get_sort_key(name)
                    operators.append((name, op, sort_key))

        # 按排序键排序
        operators.sort(key=lambda x: x[2])

        # 返回算符列表（不包含排序键）
        return [op for _, op, _ in operators]
```

### 方案 2：使用 BasisOperator 的内置排序

为 `BasisOperator` 添加排序方法：

```python
class BasisOperator(Operator):
    # ... 现有代码 ...

    # 类变量：全局排序顺序
    _global_order = []

    @classmethod
    def set_global_order(cls, order: List[str]):
        """设置全局算符排序顺序"""
        cls._global_order = order

    def _get_sort_key(self) -> int:
        """获取排序键"""
        try:
            return self._global_order.index(self._base_name)
        except ValueError:
            return 999  # 未知算符放最后

    def __lt__(self, other):
        """小于比较（用于排序）"""
        if not isinstance(other, BasisOperator):
            return NotImplemented
        return self._get_sort_key() < other._get_sort_key()

    def __le__(self, other):
        """小于等于比较"""
        if not isinstance(other, BasisOperator):
            return NotImplemented
        return self._get_sort_key() <= other._get_sort_key()
```

使用方式：

```python
# 在初始化时设置全局顺序
BasisOperator.set_global_order(['w', 'j0', 'wb', 't', 'g', 'gt', 'gw', 'gwb'])

# 然后可以直接排序
operators.sort()  # 自动按照设定的顺序排序
```

### 方案 3：混合方案（最灵活）

结合方案 1 和 2：
- `BasisOperator` 提供内置排序能力
- `OperatorEnumerator` 可以覆盖排序顺序（用于特殊情况）

## 实现步骤

### Step 1: 为 BasisOperator 添加排序支持

修改 `src/pyope/operators.py`：

```python
class BasisOperator(Operator):
    # 类变量：全局排序顺序
    _global_order = []

    @classmethod
    def set_global_order(cls, order: List[str]):
        """设置全局算符排序顺序"""
        cls._global_order = order

    @classmethod
    def get_global_order(cls) -> List[str]:
        """获取全局排序顺序"""
        return cls._global_order.copy()

    def _get_sort_key(self) -> int:
        """获取排序键"""
        if not self._global_order:
            # 如果没有设置全局顺序，使用字母顺序
            return ord(self._base_name[0]) if self._base_name else 999

        try:
            return self._global_order.index(self._base_name)
        except ValueError:
            # 未知算符放最后
            return len(self._global_order) + 999

    def __lt__(self, other):
        """小于比较（用于排序）"""
        if not isinstance(other, (BasisOperator, DerivativeOperator)):
            return NotImplemented

        # 获取基本算符
        self_base = self if isinstance(self, BasisOperator) else self.base
        other_base = other if isinstance(other, BasisOperator) else other.base

        # 比较基本算符
        if isinstance(self_base, BasisOperator) and isinstance(other_base, BasisOperator):
            return self_base._get_sort_key() < other_base._get_sort_key()

        return NotImplemented

    def __le__(self, other):
        """小于等于比较"""
        return self == other or self < other

    def __gt__(self, other):
        """大于比较"""
        if not isinstance(other, (BasisOperator, DerivativeOperator)):
            return NotImplemented
        return not (self <= other)

    def __ge__(self, other):
        """大于等于比较"""
        return self == other or self > other
```

### Step 2: 为 DerivativeOperator 添加排序支持

修改 `src/pyope/operators.py`：

```python
class DerivativeOperator(Operator):
    # ... 现有代码 ...

    def __lt__(self, other):
        """小于比较（用于排序）"""
        if isinstance(other, DerivativeOperator):
            # 两个都是导数算符，比较基本算符
            if self._base == other._base:
                # 同一个基本算符，导数阶数小的在前
                return self._order < other._order
            else:
                # 不同基本算符，使用基本算符的排序
                return self._base < other._base
        elif isinstance(other, BasisOperator):
            # 导数算符 vs 基本算符，使用基本算符的排序
            return self._base < other
        else:
            return NotImplemented

    def __le__(self, other):
        """小于等于比较"""
        return self == other or self < other

    def __gt__(self, other):
        """大于比较"""
        if isinstance(other, (BasisOperator, DerivativeOperator)):
            return not (self <= other)
        return NotImplemented

    def __ge__(self, other):
        """大于等于比较"""
        return self == other or self > other
```

### Step 3: 修改 OperatorEnumerator

修改 `src/pyope/null_states.py` 中的 `OperatorEnumerator` 类：

```python
class OperatorEnumerator:
    def __init__(self, generators: Dict[str, Dict[str, Any]]):
        self.generators = generators
        self.gen_list = list(generators.items())

        # 设置全局算符顺序（对应 Mathematica 的 opfields）
        # 这个顺序应该与 generators 的定义顺序一致
        operator_names = list(generators.keys())
        BasisOperator.set_global_order(operator_names)

    def _generate_operators_at_weight(
        self,
        weight: Fraction,
        max_derivative_order: int
    ) -> List[Any]:
        """
        生成指定权重的所有可能算符（包括导数）

        返回按照 opfields 顺序排序的算符列表
        """
        operators = []

        for name, gen_info in self.generators.items():
            base_op = gen_info['op']
            base_weight = gen_info['weight']

            # 计算需要的导数阶数
            deriv_order = weight - base_weight

            # 检查导数阶数是否有效
            if deriv_order >= 0 and deriv_order <= max_derivative_order:
                # 检查导数阶数是否为整数
                if deriv_order.denominator == 1:
                    deriv_order_int = int(deriv_order)
                    if deriv_order_int == 0:
                        operators.append(base_op)
                    else:
                        operators.append(d(base_op, deriv_order_int))

        # 排序：按照 BasisOperator 的内置排序
        # 这会自动使用 _global_order 中定义的顺序
        operators.sort()

        return operators
```

### Step 4: 测试排序

创建测试文件 `tests/test_operator_sorting.py`：

```python
#!/usr/bin/env python3
"""测试算符排序功能"""

import sys
sys.path.insert(0, '../src')

from pyope import BasisOperator, d
from fractions import Fraction

# 定义算符（按 Mathematica 的 opfields 顺序）
w = BasisOperator('w', bosonic=True, conformal_weight=Fraction(3, 2))
j0 = BasisOperator('j0', bosonic=True, conformal_weight=Fraction(1))
wb = BasisOperator('wb', bosonic=True, conformal_weight=Fraction(3, 2))
t = BasisOperator('t', bosonic=True, conformal_weight=Fraction(2))
g = BasisOperator('g', bosonic=False, conformal_weight=Fraction(3, 2))
gt = BasisOperator('gt', bosonic=False, conformal_weight=Fraction(3, 2))
gw = BasisOperator('gw', bosonic=False, conformal_weight=Fraction(2))
gwb = BasisOperator('gwb', bosonic=False, conformal_weight=Fraction(2))

# 设置全局排序顺序
BasisOperator.set_global_order(['w', 'j0', 'wb', 't', 'g', 'gt', 'gw', 'gwb'])

# 测试 1：基本算符排序
print("测试 1: 基本算符排序")
ops = [g, w, gt, j0, wb]
ops_sorted = sorted(ops)
print(f"原始顺序: {[str(op) for op in ops]}")
print(f"排序后: {[str(op) for op in ops_sorted]}")
print(f"期望: ['w', 'j0', 'wb', 'g', 'gt']")
print()

# 测试 2：导数算符排序
print("测试 2: 导数算符排序")
ops2 = [d(w), w, d(g), g, d(j0)]
ops2_sorted = sorted(ops2)
print(f"原始顺序: {[str(op) for op in ops2]}")
print(f"排序后: {[str(op) for op in ops2_sorted]}")
print(f"期望: ['w', 'd(w)', 'd(j0)', 'g', 'd(g)']")
print()

# 测试 3：权重 3/2 的所有算符
print("测试 3: 权重 3/2 的所有算符")
weight_3_2 = [w, wb, g, gt, d(j0)]
weight_3_2_sorted = sorted(weight_3_2)
print(f"排序后: {[str(op) for op in weight_3_2_sorted]}")
print(f"期望: ['w', 'd(j0)', 'wb', 'g', 'gt']")
print()

# 测试 4：权重 5/2 的所有算符
print("测试 4: 权重 5/2 的所有算符")
weight_5_2 = [d(w), d(wb), d(g), d(gt), d(d(j0)), d(t)]
weight_5_2_sorted = sorted(weight_5_2)
print(f"排序后: {[str(op) for op in weight_5_2_sorted]}")
print(f"期望: ['d(w)', 'd^2(j0)', 'd(wb)', 'd(t)', 'd(g)', 'd(gt)']")
```

## 总结

### 关键点

1. **算符顺序由 `opfields` 定义决定**
   - 在 Python 中，通过 `generators` 字典的键顺序或显式的 `operator_order` 列表定义

2. **分拆按降序排列**
   - 权重大的在左（前），权重小的在右（后）
   - 对应 `NO(op1, op2)` 的嵌套结构

3. **同一权重下按 opfields 顺序枚举**
   - 导数不改变基本算符的相对顺序
   - `∂j0` 的位置由 `j0` 在 opfields 中的位置决定

4. **去重策略**
   - 生成所有组合后使用字符串表示去重
   - 依赖简化系统处理对称性

### 推荐实现

**方案 2（使用 BasisOperator 的内置排序）** 是最优雅和可维护的方案：
- 排序逻辑封装在算符类中
- 全局顺序可以灵活配置
- 代码简洁，易于理解和测试

### 下一步

1. 实现 `BasisOperator` 和 `DerivativeOperator` 的排序方法
2. 修改 `OperatorEnumerator._generate_operators_at_weight` 使用排序
3. 编写测试验证排序正确性
4. 与 Mathematica 结果对比验证
