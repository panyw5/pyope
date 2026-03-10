# PyOPE

PyOPE 是一个面向顶点算符代数 (VOA) / 共形场论 (CFT) 的符号计算库，用 Python + SymPy 实现算符积展开 (OPE) 的自动推导与化简。

项目目标是复现并对齐 Mathematica 参考实现 `OPEdefs.m`（仓库内 `OPEdefs/` 与 `tests/*.wls` 提供了大量对照）。

## 特性

- OPE 计算：`OPE(A, B)` 返回奇异部分（各阶极点系数）
- 正规序：`NO(A, B)` / `normal_product(A, B, ...)`
- bracket：`bracket(A, B, n)`（`n=0` 约定为 `NO(A, B)`）
- 导数：`d(A)` / `dn(n, A)`，并在 OPE 里自动应用导数规则
- 代数化简：`simplify(expr)` 将表达式规范化、重排并合并同类项
- Jacobi 恒等式检查：`check_jacobi_identity` / `verify_jacobi_identity`
- Null states（实验性）：`calculate_null_states` 等工具

## 安装

本仓库使用标准 `pyproject.toml` (setuptools)。

```bash
pip install -e .
```

安装开发依赖（测试/格式化等）：

```bash
pip install -e ".[dev]"
```

安装 Notebook 相关依赖：

```bash
pip install -e ".[jupyter]"
```

## 快速上手

最小例子：定义 Virasoro 的 `T(z)T(w)` OPE。

```python
import sympy as sp

from pyope import BasisOperator, Bosonic, OPE, MakeOPE
from pyope import One, Zero, d, NO, bracket

T = BasisOperator("T", conformal_weight=2)
Bosonic(T)

c = sp.Symbol("c")

# MakeOPE([...]) 使用 Mathematica 风格：从最高极点到 (z-w)^-1
OPE[T, T] = MakeOPE(
    [
        sp.Rational(1, 2) * c * One,  # (z-w)^-4
        Zero,                          # (z-w)^-3
        2 * T,                         # (z-w)^-2
        d(T),                          # (z-w)^-1
    ]
)

print(OPE(T, T))
print("{TT}_2 =", bracket(T, T, 2))
print("(TT) =", NO(T, T))
```

注意：`Operator` 之间的普通乘法 `A * B` 会抛错（避免把 VOA 语义误写成普通乘法）；复合算符请用 `NO(A, B)` 或 `normal_product(...)`。

更多快速上手例子：

1) 提取极点系数 / 导数

```python
import sympy as sp

from pyope import BasisOperator, Bosonic, OPE, MakeOPE
from pyope import One, Zero, d, dn, bracket

T = BasisOperator("T", conformal_weight=2)
Bosonic(T)
c = sp.Symbol("c")

OPE[T, T] = MakeOPE([
    sp.Rational(1, 2) * c * One,
    Zero,
    2 * T,
    d(T),
])

tt = OPE(T, T)
print("max pole =", tt.max_pole)
print("{TT}_4 =", bracket(T, T, 4))
print("{TT}_1 =", bracket(T, T, 1))
print("{(dT)T}_3 =", bracket(d(T), T, 3))
print("d^2 T =", dn(2, T))
```

2) 复合算符与化简（使用 `normal_product` 构造嵌套正规序）

```python
from pyope import BasisOperator, Bosonic
from pyope import NO, normal_product, simplify

A = BasisOperator("A")
B = BasisOperator("B")
Bosonic(A, B)

expr = normal_product(B, A, B)  # NO(B, NO(A, B))

print("raw =", expr)
print("simplified =", simplify(expr))
```

3) W3 代数（最小骨架）

```python
import sympy as sp

from pyope import BasisOperator, Bosonic, OPE, MakeOPE
from pyope import One, Zero, d, dn, NO

T = BasisOperator("T", conformal_weight=2)
W = BasisOperator("W", conformal_weight=3)
Bosonic(T, W)

c = sp.Symbol("c")
beta = sp.Symbol("beta")

OPE[T, T] = MakeOPE([
    sp.Rational(1, 2) * c * One,
    Zero,
    2 * T,
    d(T),
])

OPE[T, W] = MakeOPE([
    3 * W,
    d(W),
])

Lambda = NO(T, T) - sp.Rational(3, 10) * dn(2, T)
OPE[W, W] = MakeOPE([
    c * One,
    Zero,
    2 * T,
    d(T),
    2 * beta * Lambda + sp.Rational(3, 10) * dn(2, T),
    beta * d(Lambda) + sp.Rational(1, 15) * dn(3, T),
])

ww = OPE(W, W)
print("max pole =", ww.max_pole)
print("(z-w)^-6 coeff =", ww.pole(6))
print("(z-w)^-2 coeff =", ww.pole(2))
```

4) Jacobi 恒等式快速验证

```python
import sympy as sp

from pyope import BasisOperator, Bosonic, OPE, MakeOPE
from pyope import One, Zero, d, verify_jacobi_identity

T = BasisOperator("T", conformal_weight=2)
Bosonic(T)
c = sp.Symbol("c")

OPE[T, T] = MakeOPE([
    sp.Rational(1, 2) * c * One,
    Zero,
    2 * T,
    d(T),
])

print("Jacobi(T,T,T) =", verify_jacobi_identity(T, T, T))
```

## API 概览

- `pyope.api`：`OPE`, `NO`, `MakeOPE`, `bracket`, `normal_product`
- `pyope.operators`：`BasisOperator`, `d`, `dn`
- `pyope.simplify`：`simplify`
- `pyope.jacobi`：Jacobi 恒等式验证

`OPE(A, B)` 的返回值是 `OPEData`，可用 `.pole(n)` 取出 $(z-w)^{-n}$ 系数。

## 示例与文档

- Notebooks：`demo/`（Virasoro, Kac-Moody, Jacobi identity 等）
- 文档碎片：`docs/normal_product.md`
- 测试快速参考：`tests/QUICK_REFERENCE.md`

## 运行测试

```bash
python -m pytest
```

如果你想只跑参考对照类的测试：

```bash
python -m pytest -m mathematica_ref
```

更多测试结构说明见：`tests/TEST_FRAMEWORK.md`。

## 参考

- K. Thielemans, "An Algorithmic Approach to Operator Product Expansions, W-algebras and W-strings", arXiv:hep-th/9506159

## License

MIT
