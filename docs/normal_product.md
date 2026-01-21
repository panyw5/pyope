# normal_product() 函数

## 功能说明

`normal_product()` 是一个辅助函数，用于将多个算符的乘积转化为嵌套的正规序（normal ordered product）。

## 语法

```python
from pyope import normal_product

result = normal_product(A, B, C, D, ...)
```

## 转换规则

`normal_product(A, B, C, D)` 等价于 `NO(A, NO(B, NO(C, D)))`

即从右向左构建嵌套的正规序：

```
A * B * C * D  →  (A(B(CD)))
```

## 使用示例

### 基本用法

```python
from pyope import BasisOperator, NO, normal_product

T = BasisOperator("T", 2)
J = BasisOperator("J", 1)
W = BasisOperator("W", 3)

# 两个算符
result = normal_product(T, J)
# 等价于: NO(T, J)

# 三个算符
result = normal_product(T, J, W)
# 等价于: NO(T, NO(J, W))

# 四个算符
L = BasisOperator("L", 1)
result = normal_product(T, J, W, L)
# 等价于: NO(T, NO(J, NO(W, L)))
```

### 特殊情况

```python
# 空参数列表返回单位算符
normal_product()  # → One

# 单个算符返回自身
normal_product(T)  # → T

# 包含零算符返回零
normal_product(T, Zero, J)  # → Zero

# 包含单位算符会被简化
normal_product(T, One, J)  # → NO(T, J)
```

### 包含导数算符

```python
from pyope import d

result = normal_product(d(T), J, T)
# 等价于: NO(∂T, NO(J, T))
```

### 包含标量系数

```python
import sympy as sp

c = sp.Symbol("c")
k = sp.Symbol("k")

result = normal_product(c * T, k * J)
# 标量系数会自动提取: c*k*NO(T, J)
```

### 多个算符

```python
# 使用解包语法处理多个算符
ops = [BasisOperator(f"O{i}", 1) for i in range(6)]
result = normal_product(*ops)
# 等价于: NO(O0, NO(O1, NO(O2, NO(O3, NO(O4, O5)))))
```

## 应用场景

1. **构建复合算符**：简化多个算符的正规序乘积表达式
2. **W-代数**：构建高阶生成元
3. **Sugawara 构造**：构建应力张量
4. **任何需要多个算符正规序乘积的场合**

## 优势

相比手动嵌套 `NO()` 调用：

- ✅ **更简洁**：`normal_product(A, B, C, D)` vs `NO(A, NO(B, NO(C, D)))`
- ✅ **更易读**：直接表达算符乘积的意图
- ✅ **更灵活**：支持任意数量的算符
- ✅ **自动处理**：自动处理特殊情况（零、单位算符、标量系数）

## 实现细节

函数从右向左构建嵌套的 `NO`：

```python
def normal_product(*operators):
    if len(operators) == 0:
        return One
    if len(operators) == 1:
        return operators[0]
    
    result = operators[-1]
    for op in reversed(operators[:-1]):
        result = NO(op, result)
    return result
```

## 相关函数

- `NO(A, B)`: 计算两个算符的正规序乘积
- `OPE(A, B)`: 计算算符积展开
- `bracket(A, B, n)`: 计算 bracket {AB}_n

## 参考

- 测试文件: `tests/test_normal_product.py`
- 演示 notebook: `demo/normal_product_demo.ipynb`
