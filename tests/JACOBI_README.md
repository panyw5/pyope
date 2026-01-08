# Jacobi 恒等式测试 - 快速开始

## 一句话总结

✅ pyope 的 Jacobi 恒等式实现已通过完整测试，与 Mathematica OPEdefs.m 完全一致。

---

## 快速运行

```bash
# Python 测试
pytest tests/test_jacobi_virasoro.py -v

# Mathematica 参考
wolframscript tests/ref_jacobi_virasoro.wls
```

---

## 测试结果

```
6 passed in 0.15s ✅
```

---

## 核心 API

```python
from pyope import check_jacobi_identity, verify_jacobi_identity

# 方式 1: 获取完整矩阵
result = check_jacobi_identity(T, T, T)
# [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], ...]

# 方式 2: 布尔判断
is_valid = verify_jacobi_identity(T, T, T)
# True
```

---

## 示例: Virasoro 代数

```python
import sympy as sp
from pyope import BasisOperator, OPE, d, One

# 定义能动张量 T
c = sp.Symbol('c')
T = BasisOperator("T", bosonic=True, conformal_weight=2)

# 定义 OPE: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])

# 验证 Jacobi 恒等式
assert verify_jacobi_identity(T, T, T) == True  ✅
```

---

## 文档

- **快速总结**: `JACOBI_TEST_SUMMARY.md`
- **详细报告**: `JACOBI_TEST_REPORT.md`
- **测试说明**: `JACOBI_IDENTITY_TEST.md`
- **测试索引**: `TEST_INDEX.md`

---

## 数学背景

Jacobi 恒等式是顶点算符代数的基本恒等式：

```
[A, [B, C]_q]_m - (-1)^(|A||B|) [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p} = 0
```

保证了 OPE 计算的结合律和对称性。

---

## 参考

- **实现**: `src/pyope/jacobi.py`
- **参考**: OPEdefs.m 第 1601-1637 行
- **理论**: V. Kac, "Vertex Algebras for Beginners", Chapter 3

---

**状态**: ✅ 生产就绪
**版本**: pyope 0.1.0
**日期**: 2026-01-07
