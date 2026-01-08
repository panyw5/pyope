# Jacobi 恒等式测试总结

## 快速概览

✅ **状态**: 所有测试通过
📊 **测试数量**: 6 个测试用例
⏱️ **执行时间**: ~0.15 秒
🎯 **覆盖率**: Virasoro 代数 + 简单流代数

---

## 测试文件

| 文件 | 类型 | 功能 |
|-----|------|------|
| `src/pyope/jacobi.py` | 实现 | Jacobi 恒等式核心函数 |
| `tests/ref_jacobi_virasoro.wls` | 参考 | Mathematica 参考实现 |
| `tests/test_jacobi_virasoro.py` | 测试 | Python 测试套件 |

---

## 快速运行

### Mathematica 参考测试
```bash
wolframscript tests/ref_jacobi_virasoro.wls
```

**预期输出**: 5×5 全零矩阵，Jacobi 恒等式成立

### Python 测试
```bash
# pytest 方式
pytest tests/test_jacobi_virasoro.py -v

# 独立运行
python tests/test_jacobi_virasoro.py
```

**预期结果**: 6 passed in ~0.15s

---

## 测试内容

### 1. Virasoro OPE 结构验证 ✓
验证 `T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)`

### 2. Jacobi 恒等式计算 ✓
计算 `check_jacobi_identity(T, T, T)` 返回 5×5 全零矩阵

### 3. 便捷函数验证 ✓
`verify_jacobi_identity(T, T, T)` 返回 `True`

### 4. 维度验证 ✓
结果矩阵维度为 5×5，与 Mathematica 一致

### 5. Mathematica 对比 ✓
pyope 结果与 Mathematica 逐元素对比，完全一致

### 6. 流代数测试 ✓
简单流代数 `J(z)J(w) = k/(z-w)^2` 满足 Jacobi 恒等式

---

## 核心 API

### check_jacobi_identity
```python
from pyope import check_jacobi_identity

result = check_jacobi_identity(A, B, C, simplify_func=sp.expand)
# 返回: List[List[Any]] - Jacobi 恒等式矩阵
```

### verify_jacobi_identity
```python
from pyope import verify_jacobi_identity

is_valid = verify_jacobi_identity(A, B, C, simplify_func=sp.expand)
# 返回: bool - True 表示恒等式成立
```

---

## 测试示例

### Virasoro 代数
```python
import sympy as sp
from pyope import BasisOperator, OPE, d, One, check_jacobi_identity

# 定义算符
c = sp.Symbol('c')
T = BasisOperator("T", bosonic=True, conformal_weight=2)

# 定义 OPE
OPE[T, T] = OPE.make([c/2*One, 0, 2*T, d(T)])

# 检查 Jacobi 恒等式
result = check_jacobi_identity(T, T, T)
# result = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], ...]
```

### 简单流代数
```python
# 定义流算符
k = sp.Symbol('k')
J = BasisOperator("J", bosonic=True, conformal_weight=1)

# 定义 OPE
OPE[J, J] = OPE.make([k*One, 0])

# 验证 Jacobi 恒等式
assert verify_jacobi_identity(J, J, J) == True
```

---

## 数学背景

### Jacobi 恒等式
```
[A, [B, C]_q]_m - (-1)^(|A||B|) [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p} = 0
```

**物理意义**: 保证顶点算符代数的结合律和对称性

### Virasoro 代数
```
T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w) + ...
```

**物理意义**: 共形场论的核心对称性

---

## 实现对比

| 特性 | OPEdefs.m | pyope | 状态 |
|-----|-----------|-------|------|
| Parity 符号 | `(-1)^(OPEParity[A] OPEParity[B])` | `(-1)**(parity_A * parity_B)` | ✓ |
| 嵌套 OPE | `OPE[A, OPEPole[n][BC]]` | `OPE(A, bracket(B, C, n))` | ✓ |
| A==B 优化 | `If[SameQ[A,B], ...]` | `if A == B: ...` | ✓ |
| 二项式系数 | `Binomial[n-1,p-1]` | `binomial(n-1, p-1)` | ✓ |
| 结果矩阵 | `Table[..., {m, maxm}, {n, maxn}]` | `for m in range(1, max_m+1): ...` | ✓ |

---

## 参考资料

### 代码参考
- **OPEdefs.m**: 第 1601-1637 行 - `OPEJacobi` 实现
- **voa-manual.md**: Section 3.3 - Implementation details

### 理论参考
- Di Francesco et al., "Conformal Field Theory" (1997)
- V. Kac, "Vertex Algebras for Beginners" (1998)

---

## 下一步

### 建议的扩展测试
- [ ] W 代数（W_3, W_4）
- [ ] 超对称代数（N=1, N=2）
- [ ] 费米算符
- [ ] 混合玻色-费米算符

### 性能优化
- [ ] 缓存中间结果
- [ ] 并行化计算
- [ ] 优化符号简化

---

**最后更新**: 2026-01-07
**测试版本**: pyope 0.1.0
**参考版本**: OPEdefs 3.1 (beta 4)
