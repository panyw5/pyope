# BasisOperator API 改进

## 改进内容

为了让 API 更直观，`BasisOperator` 现在支持 `fermionic` 参数。

## 新的推荐用法

```python
from pyope import BasisOperator

# 玻色子（默认）
T = BasisOperator("T", conformal_weight=2)
J = BasisOperator("J", conformal_weight=1)

# 费米子（推荐用法）
b = BasisOperator("b", fermionic=True, conformal_weight=2)
c = BasisOperator("c", fermionic=True, conformal_weight=-1)
psi = BasisOperator("ψ", fermionic=True, conformal_weight=3/2)
```

## 旧的用法（仍然支持）

```python
# 旧的方式（仍然有效）
b = BasisOperator("b", bosonic=False, conformal_weight=2)
c = BasisOperator("c", bosonic=False, conformal_weight=-1)
```

## 参数优先级

- 如果只指定 `fermionic=True`，算符为费米子
- 如果只指定 `bosonic=False`，算符为费米子
- 如果同时指定，`fermionic` 参数优先
- 如果都不指定，默认为玻色子

## 示例对比

### 之前（不够直观）

```python
# 需要记住 bosonic=False 表示费米子
b = BasisOperator("b", bosonic=False)  # 费米子
c = BasisOperator("c", bosonic=False)  # 费米子
```

### 现在（更直观）

```python
# 直接说明是费米子
b = BasisOperator("b", fermionic=True)  # 费米子
c = BasisOperator("c", fermionic=True)  # 费米子

# 玻色子不需要指定（默认）
T = BasisOperator("T")  # 玻色子
J = BasisOperator("J")  # 玻色子
```

## 完整示例

```python
from pyope import BasisOperator, OPE, NO, Fermionic, One
from fractions import Fraction

# bc ghost 系统（费米子）
b = BasisOperator("b", fermionic=True, conformal_weight=2)
c = BasisOperator("c", fermionic=True, conformal_weight=-1)

# 注册（可选，但推荐）
Fermionic(b, c)

# 定义 OPE
OPE[b, c] = OPE.make([One])

# 计算
print(NO(b, b))  # 输出: 0（费米子自身正规序）
```

## 向后兼容性

所有旧代码仍然可以正常工作，不需要修改。新的 `fermionic` 参数只是提供了更直观的替代方式。
