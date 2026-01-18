# 正规序的正则化（Canonicalization）

根据 Thielemans 论文中的 eq (2.3.16)，正规序应该按照算符的正则顺序自动排列。

## 问题描述

在之前的实现中，`NO(β, γ)` 和 `NO(γ, β)` 被视为不同的对象：

```python
from fractions import Fraction
from pyope import BasisOperator, NO

beta = BasisOperator('β', conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', conformal_weight=Fraction(-1, 2))

print(NO(beta, gamma))  # 输出: NO(β,γ)
print(NO(gamma, beta))  # 输出: NO(γ,β)
```

但是，当这两个算符之间没有定义 OPE 时（即 OPE = 0），根据 eq (2.3.16)：

$$[BA]_q = (-1)^{|A||B|} \sum_{l \geq q} \frac{(-1)^l}{(l-q)!} \partial^{(l-q)} [AB]_l$$

对于 $q=0$（正规序）且所有 $[AB]_l = 0$：

$$[BA]_0 = (-1)^{|A||B|} [AB]_0$$

对于玻色算符（$|A||B| = 0$），这意味着 `NO(B,A) = NO(A,B)`。

## 解决方案

我们实现了 `simplify()` 函数，它会自动将正规序算符按照正则顺序排列：

```python
from pyope import simplify

# 化简后，两者变成相同的形式
print(simplify(NO(beta, gamma)))  # 输出: NO(β,γ)
print(simplify(NO(gamma, beta)))  # 输出: NO(β,γ)
```

## 算符排序规则

算符按照以下规则排序：

1. **注册顺序优先**：如果算符已经通过 `Bosonic()` 或 `Fermionic()` 注册，按注册顺序排列
2. **字典序作为后备**：如果算符未注册，按字符串的字典序排列
3. **导数阶数**：对于同一个算符的不同导数，阶数小的在前（例如 `J` 在 `∂J` 前面）

## 示例

### 1. 玻色算符（未定义 OPE）

```python
from pyope import BasisOperator, NO, simplify

A = BasisOperator('A', conformal_weight=1)
B = BasisOperator('B', conformal_weight=1)

# 未化简
print(NO(B, A))  # NO(B,A)

# 化简后按字典序排列
print(simplify(NO(B, A)))  # NO(A,B)
```

### 2. 费米算符（未定义 OPE）

```python
psi = BasisOperator('ψ', conformal_weight=1/2, bosonic=False)
chi = BasisOperator('χ', conformal_weight=1/2, bosonic=False)

# 费米算符有额外的符号因子 (-1)^{|A||B|} = -1
print(simplify(NO(chi, psi)))  # -NO(ψ,χ)
```

### 3. 有定义 OPE 的算符

当算符之间有定义 OPE 时，`simplify()` 会应用完整的 eq (2.3.16) 公式：

```python
from pyope import BasisOperator, NO, OPE, MakeOPE, d, simplify

T = BasisOperator('T', conformal_weight=2)
c = BasisOperator('c', conformal_weight=0)

# 定义 Virasoro OPE
OPE[T, T] = MakeOPE([c/2, 0, 2*T, d(T)])

# 计算 NO(T, T)
print(simplify(NO(T, T)))  # NO(T,T) - 已经是正确顺序
```

### 4. 导数算符

```python
J = BasisOperator('J', conformal_weight=1)

# 导数阶数小的在前
print(simplify(NO(d(J), J)))  # NO(J,∂J)
```

### 5. 列表化简

```python
beta = BasisOperator('β', conformal_weight=3/2)
gamma = BasisOperator('γ', conformal_weight=-1/2)

# 化简列表中的所有项
no_list = [NO(beta, gamma), NO(gamma, beta)]
simplified = [simplify(no) for no in no_list]

print(simplified)  # [NO(β,γ), NO(β,γ)]
```

## 数学背景

### Thielemans 论文 eq (2.3.16)

交换公式：

$$[BA]_q = (-1)^{|A||B|} \sum_{l \geq q} \frac{(-1)^l}{(l-q)!} \partial^{(l-q)} [AB]_l$$

其中：
- $[AB]_n$ 是 OPE $A(z)B(w)$ 中 $(z-w)^{-n}$ 的系数
- $|A|, |B|$ 是算符的 parity（0 表示玻色子，1 表示费米子）
- $\partial$ 表示对 $w$ 求导

### 正规序的定义

正规序定义为 OPE 的 0 阶极点：

$$NO(A,B) = [AB]_0$$

### 特殊情况

当 OPE 未定义（所有极点为 0）时：

$$[BA]_0 = (-1)^{|A||B|} [AB]_0$$

这意味着：
- **玻色算符**：$NO(B,A) = NO(A,B)$
- **费米算符**：$NO(B,A) = -NO(A,B)$

## 使用建议

1. **总是使用 `simplify()`**：在比较或输出正规序表达式之前，使用 `simplify()` 将其规范化
2. **注册算符**：使用 `Bosonic()` 或 `Fermionic()` 注册算符，以控制排序顺序
3. **定义 OPE**：为基本算符定义 OPE，以获得正确的交换关系

## 实现细节

`simplify()` 函数的工作流程：

1. 检查算符顺序（使用 `compare_operators()`）
2. 如果顺序不正确（`order < 0`），应用 eq (2.3.16)
3. 计算 OPE(A, B)
4. 如果所有极点为 0，返回 `swap_sign * NO(A, B)`
5. 否则，应用完整公式计算修正项

## 相关函数

- `NO(A, B)`: 创建正规序算符
- `simplify(expr)`: 化简表达式，包括正规序的正则化
- `OPE(A, B)`: 计算算符积展开
- `compare_operators(A, B)`: 比较算符顺序

## 参考文献

- Thielemans, K. (1994). "An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings", PhD thesis, eq (2.3.16)
