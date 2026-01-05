本项目的目标是开发 `python` 版本的 `OPEdefs.m`，用于计算 vertex operator algebras (VOA) 中的算符积展开 (OPE)。

# 核心功能初步设想 (待完善)
- 公共函数/算符: `OPE`，`NO`，`bracket`，`d`，`dn`
- 父类: `Operator`
- VOA 元素类: `LocalOperator`
  - 能够做加法 `A+B`，减法 `A-B`，结果还是 `LocalOperator`
    并且能够自动做基本的合并同类项，比如 `A + (-1)*A = Zero`
  - 能够做数乘 $3 * A$ ，结果还是 `LocalOperator`
  - 能够求导，结果还是 `LocalOperator`
- 针对不同的应用场景设置一些 `LocalOperator` 子类，添加更多的方便的性质与方法
  - 特殊的 `LocalOperator`：`One`，`Zero` (`0 * A = Zero`)
  - 比如 `BasicLocalOperator` 为用户指定的基本的算符
  - 算符的各种复合，比如 normal order product $(AB)$，bracket $\{AB\}$ 等
- `OPE(A, B)`: 计算两个 VOA 元素 A 和 B 的算符积展开 $A(z)B(w)$，结果是 `OPEData`，包含各个极点的信息
  - 从 `OPE(A, B)` 可以提取 `.pole(n)` 得到 $(z-w)^{-n}$ 的系数 $\{AB\}_n$ (`LocalOperator` 类)
- `NO(A, B)` (`LocalOperator` 类): 计算 A 和 B 的正则序列化 (normal ordering) $(AB)(z)$。
- `bracket(A, B)(n)` (`LocalOperator` 类): 计算 A 和 B 的 $\{AB\}(z)$
- 求导：`d(A)` 或者 `dn(n, A)`  (`LocalOperator` 类): 计算 VOA 元素 A 的导数 $∂A(z)$ 或者 $∂^n A(z)$。
  - 可加性：$∂(A+B) = ∂A + ∂B$
  - 线性性：$∂(cA) = c∂A$
- `A[n]` (`Operator`): 获取 `LocalOperator` `A` 的第 n 个模态 (mode) $A_n$。
  - 模态的加减和数乘依然是模态

