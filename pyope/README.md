这是一个 `python` 版本的 OPEdef，用于计算 vertex operator algebras (VOA) 中的算符积展开 (OPE)。

# 核心功能初步设想 (待完善)
- 公共函数/算符: `OPE`，`NO`，`bracket`，`d`，`dn`
- VOA 元素类: `Operator` 或者 `Op`
- `OPE(A, B)`: 计算两个 VOA 元素 A 和 B 的算符积展开 $A(z)B(w)$。
- `NO(A, B)`: 计算 A 和 B 的正则序列化 (normal ordering) $(AB)(z)$。
- `bracket(A, B)`: 计算 A 和 B 的 $\{AB\}(z)$
- `d(A)` 或者 `dn(n, A)`: 计算 VOA 元素 A 的导数 $∂A(z)$ 或者 $∂^n A(z)$。
- `A.mode(n)`: 获取 VOA 元素 A 的第 n 个模态 (mode) $A_n$。
