# PyOPE LaTeX 显示测试

## 在 Jupyter Notebook 中显示 LaTeX 格式的 OPE

有三种方式显示 OPE：

### 方法 1：直接显示对象（推荐）

```python
from pyope.api import OPE
from pyope.operators import BasisOperator, d
from pyope.constants import One
import sympy as sp

# 定义 Virasoro OPE
T = BasisOperator('T', bosonic=True)
c = sp.Symbol('c')
OPE[T, T] = OPE.make([c/2 * One, 0, 2*T, d(T)])

# 方法 1：直接显示（会自动渲染 LaTeX）
result = OPE(T, T)
result  # 在 Jupyter 中，这会自动调用 _repr_latex_() 并渲染
```

### 方法 2：使用 display() 方法

```python
# 方法 2：使用 display() 方法（显式渲染）
result = OPE(T, T)
result.display()  # 显式调用 IPython.display.Math
```

### 方法 3：获取 LaTeX 字符串

```python
# 方法 3：获取 LaTeX 字符串（用于自定义显示）
result = OPE(T, T)
latex_str = result.to_latex()
print(latex_str)  # 打印 LaTeX 代码

# 或者手动渲染
from IPython.display import display, Math
display(Math(latex_str))
```

## 注意事项

⚠️ **不要使用 `print(result)`**
- `print()` 会调用 `__str__()` 方法，不会渲染 LaTeX
- 应该直接输入 `result` 或使用 `result.display()`

✅ **正确用法：**
```python
result = OPE(T, T)
result  # 自动 LaTeX 渲染
```

或

```python
result = OPE(T, T)
result.display()  # 显式 LaTeX 渲染
```

❌ **错误用法：**
```python
result = OPE(T, T)
print(result)  # 只会打印文本，不会渲染 LaTeX
```

## 示例输出

当正确使用时，你会看到漂亮的数学公式：

$$\frac{\frac{c}{2}}{(z-w)^{4}} + \frac{2 T}{(z-w)^{2}} + \frac{\partial T}{(z-w)}$$

而不是原始的 LaTeX 代码。
