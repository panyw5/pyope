# 自动类型转换功能

## 概述

从版本 0.1.0 开始，`pyope` 自动将用户输入的 Python 数字（整数、浮点数）转换为 SymPy 的精确类型（Integer、Rational），确保符号计算的精确性。

## 问题背景

在 Python 中，`3/2` 会计算为浮点数 `1.5`，而不是有理数。这在符号计算中会导致精度问题：

```python
# Python 的行为
>>> 3/2
1.5  # 浮点数

>>> 26/2
13.0  # 浮点数
```

在 SymPy 中，我们希望保持精确的有理数表示：

```python
# SymPy 的期望
>>> sp.Rational(3, 2)
3/2  # 精确的有理数

>>> sp.Integer(13)
13  # 精确的整数
```

## 解决方案

`pyope` 现在会自动转换所有系数，用户可以直接使用 Python 的数字：

```python
from pyope import BasisOperator, MakeOPE, Bosonic, OPE, One

T = BasisOperator("T", conformal_weight=2)
Bosonic(T)

# 方式 1: 直接使用 Python 除法（推荐！）
OPE[T, T] = MakeOPE([26/2 * One, 0, 2*T])
# 自动转换为: 13*One（精确整数）

# 方式 2: 使用浮点数
OPE[T, T] = MakeOPE([0.5 * One, 0, 2*T])
# 自动转换为: One/2（精确有理数）

# 方式 3: 使用 SymPy 类型（也可以）
from sympy import Rational
OPE[T, T] = MakeOPE([Rational(1, 2) * One, 0, 2*T])
# 保持为: One/2
```

## 转换规则

`sympify_coefficient` 函数会自动应用以下转换规则：

1. **Python int** → `sympy.Integer`
   ```python
   13 → Integer(13)
   ```

2. **Python float** → `sympy.Rational` 或 `sympy.Integer`
   ```python
   13.0 → Integer(13)
   0.5 → Rational(1, 2)
   1.5 → Rational(3, 2)
   ```

3. **SymPy 表达式中的 Float** → 递归转换
   ```python
   13.0 * One → 13 * One
   0.5 * T → T/2
   ```

4. **算符对象** → 保持不变
   ```python
   T → T（不转换）
   ```

## 使用示例

### 示例 1: Virasoro 代数

```python
from pyope import BasisOperator, MakeOPE, Bosonic, OPE, One, d
from sympy import Symbol

T = BasisOperator("T", conformal_weight=2)
Bosonic(T)

c = Symbol("c")

# 定义 Virasoro OPE（使用 Python 除法）
OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])

# 计算 OPE
result = OPE(T, T)
print(result.pole(4))  # c/2*One（精确有理数）
```

### 示例 2: 特定中心荷

```python
# c = 26（临界弦理论）
c_val = 26
OPE[T, T] = MakeOPE([c_val/2 * One, 0, 2*T, d(T)])

result = OPE(T, T)
assert result.pole(4) == 13 * One  # ✓ 通过！

# c = 1（自由玻色子）
c_val = 1
OPE[T, T] = MakeOPE([c_val/2 * One, 0, 2*T, d(T)])

result = OPE(T, T)
assert result.pole(4) == Rational(1, 2) * One  # ✓ 通过！
```

### 示例 3: β-γ 系统

```python
from pyope import BasisOperator, MakeOPE, Fermionic, OPE, d

β = BasisOperator("β", conformal_weight=3/2)  # 直接使用 3/2
γ = BasisOperator("γ", conformal_weight=-1/2)  # 直接使用 -1/2
Fermionic(β, γ)

# 定义 OPE
OPE[β, γ] = MakeOPE([One])  # β(z)γ(w) ~ 1/(z-w)

# 应力张量
T = -3/2 * NO(β, d(γ)) + 1/2 * NO(d(β), γ)
# 系数 -3/2 和 1/2 会自动转换为精确有理数
```

## 技术细节

### 实现位置

- **核心函数**: `src/pyope/utils.py::sympify_coefficient()`
- **应用位置**: 
  - `OPEData.__init__()` - 创建 OPE 数据时
  - `OPEData.from_list()` - 从列表创建时

### 转换算法

使用 SymPy 的 `nsimplify()` 函数，它能够：
- 识别浮点数的精确有理数表示
- 递归处理复杂表达式（Mul、Add）
- 保持算符对象不变

```python
def sympify_coefficient(coeff):
    if isinstance(coeff, float):
        # 13.0 → 13, 0.5 → 1/2
        return sp.nsimplify(coeff)
    elif isinstance(coeff, sp.Mul):
        # 13.0 * One → 13 * One
        return sp.nsimplify(coeff)
    # ... 其他情况
```

## 最佳实践

### ✅ 推荐做法

```python
# 1. 直接使用 Python 数字（最简洁）
OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])

# 2. 使用变量
c_val = 26
OPE[T, T] = MakeOPE([c_val/2 * One, 0, 2*T, d(T)])

# 3. 使用 SymPy 符号（参数化）
from sympy import Symbol
c = Symbol("c")
OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
```

### ⚠️ 注意事项

```python
# 避免使用字符串
# ❌ 不推荐
OPE[T, T] = MakeOPE(["c/2", 0, "2*T", "d(T)"])

# ✅ 推荐
from sympy import Symbol
c = Symbol("c")
OPE[T, T] = MakeOPE([c/2 * One, 0, 2*T, d(T)])
```

## 测试验证

所有相关测试都已通过：

```bash
$ pytest tests/test_virasoro_voa.py -v
# 17 passed in 0.10s

$ pytest tests/test_constants.py -v
# All tests passed
```

## 相关问题

- Issue: "浮点数 vs 整数比较问题"
- 修复: 自动类型转换功能
- 测试: `tests/test_virasoro_voa.py::TestVirasoroNumerical`

## 参考资料

- SymPy 文档: [nsimplify](https://docs.sympy.org/latest/modules/simplify/simplify.html#sympy.simplify.simplify.nsimplify)
- Python 除法: [PEP 238](https://www.python.org/dev/peps/pep-0238/)
