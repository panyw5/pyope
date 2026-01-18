# pyope 快速参考：费米子处理

## 定义算符

### 推荐方式（新）

```python
from pyope import BasisOperator

# 玻色子（默认）
T = BasisOperator("T", conformal_weight=2)
J = BasisOperator("J", conformal_weight=1)

# 费米子（推荐）
b = BasisOperator("b", fermionic=True, conformal_weight=2)
c = BasisOperator("c", fermionic=True, conformal_weight=-1)
psi = BasisOperator("ψ", fermionic=True, conformal_weight=3/2)
```

### 旧方式（仍然支持）

```python
# 费米子
b = BasisOperator("b", bosonic=False, conformal_weight=2)
c = BasisOperator("c", bosonic=False, conformal_weight=-1)
```

## 费米子自身正规序

### 自动化简

```python
from pyope import NO, OPE, Fermionic, One

b = BasisOperator("b", fermionic=True)
c = BasisOperator("c", fermionic=True)

Fermionic(b, c)
OPE[b, c] = OPE.make([One])

# 费米子自身正规序自动化简
NO(b, b)  # 返回: 0 ✓
```

### 从 OPE 计算

```python
import sympy as sp

G = BasisOperator("G", fermionic=True, conformal_weight=3/2)
T = BasisOperator("T", conformal_weight=2)

c = sp.Symbol("c")
OPE[G, G] = OPE.make([2*c/3*One, 0, 2*T])

# 根据公式 (2.3.17) 从奇数阶极点计算
NO(G, G)  # 返回: ∂T ✓
```

## 关键点

✅ **费米子自身正规序**：通过 OPE 自动计算，不需要手动处理  
✅ **玻色子自身正规序**：返回 `NormalOrderedOperator(A, A)`  
✅ **向后兼容**：旧代码无需修改  
✅ **与 OPEdefs.m 一致**：行为完全匹配

## 理论公式

对于费米子 $\psi$：

$$[\psi\psi]_0 = \sum_{k \geq 0} a_k \partial^{2k+1} [\psi\psi]_{2k+1}$$

其中：
- $a_0 = 1/2$
- $a_1 = -1/24$
- $a_2 = 1/240$
- $a_3 = -17/40320$

## 示例系统

### bc ghost 系统

```python
b = BasisOperator("b", fermionic=True, conformal_weight=2)
c = BasisOperator("c", fermionic=True, conformal_weight=-1)
OPE[b, c] = OPE.make([One])
```

### βγ 系统

```python
beta = BasisOperator("β", fermionic=True, conformal_weight=3/2)
gamma = BasisOperator("γ", fermionic=True, conformal_weight=-1/2)
OPE[beta, gamma] = OPE.make([-One])
```

### N=1 超共形代数

```python
G = BasisOperator("G", fermionic=True, conformal_weight=3/2)
T = BasisOperator("T", conformal_weight=2)
c = sp.Symbol("c")
OPE[G, G] = OPE.make([2*c/3*One, 0, 2*T])
```
