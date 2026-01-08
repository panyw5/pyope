# Jacobi 恒等式测试文档

## 测试概述

本测试验证了 pyope 库中 Jacobi 恒等式的实现，使用 Virasoro 代数作为测试对象。

## 测试文件

### 1. Mathematica 参考测试
- **文件**: `/Users/lelouch/pyope/tests/ref_jacobi_virasoro.wls`
- **功能**: 使用 OPEdefs.m 计算 Virasoro 算符 T 的 Jacobi 恒等式
- **结果**: 5x5 全零矩阵（Jacobi 恒等式成立）

### 2. Python 实现
- **文件**: `/Users/lelouch/pyope/src/pyope/jacobi.py`
- **功能**: 实现 `check_jacobi_identity` 和 `verify_jacobi_identity` 函数
- **参考**: OPEdefs.m 第 1601-1637 行的 OPEJacobi 实现

### 3. Python 测试
- **文件**: `/Users/lelouch/pyope/tests/test_jacobi_virasoro.py`
- **功能**: 测试 Jacobi 恒等式在 Virasoro 代数上的应用
- **结果**: 所有测试通过

## Jacobi 恒等式

Jacobi 恒等式是顶点算符代数的基本恒等式之一：

```
[A, [B, C]_q]_m - (-1)^(|A||B|) [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p} = 0
```

其中：
- `[A, B]_n` 表示 OPE(A, B) 的第 n 阶极点
- `|A|` 表示算符 A 的 parity（0 表示玻色子，1 表示费米子）
- `C(n, k)` 表示二项式系数

## Virasoro 代数

Virasoro 代数是共形场论中的核心代数结构，其能动张量 T 满足：

```
T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w) + ...
```

其中 c 是中心荷（central charge）。

## 测试结果

### Mathematica 参考结果

```
Loading OPEdefs.m...
OPEdefs Version 3.1 (beta 4) by Kris Thielemans
Defining Virasoro operator T with central charge c...
OPE[T,T] = OPEData[{(c*One)/2, 0, 2*T, Derivative[1][T]}]

Verifying T(z)T(w) OPE structure:
  Pole 4: (c*One)/2 (should be c/2*One)
  Pole 3: 0 (should be 0)
  Pole 2: 2*T (should be 2*T)
  Pole 1: Derivative[1][T] (should be T')

Computing Jacobi identity OPEJacobi[T, T, T]...
Jacobi identity result:
  Type: List
  Dimensions: {5, 5}
  Full result: {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}

SUCCESS: All entries are zero. Jacobi identity holds!
Result: PASS
```

### Python 测试结果

```
======================================================================
Jacobi Identity Test for Virasoro Algebra
======================================================================

Test 1: Virasoro OPE Structure
  Pole 4: One*c/2 ✓
  Pole 3: 0 ✓
  Pole 2: 2*T ✓
  Pole 1: ∂T ✓

Test 2: Jacobi Identity check_jacobi_identity(T, T, T)
  Result dimensions: 5 x 5
  All entries are zero! ✓

Test 3: verify_jacobi_identity(T, T, T)
  Result: True ✓

Test 4: Jacobi Identity Dimensions
  Expected: 5 x 5
  Actual: 5 x 5 ✓

Test 5: Comparison with Mathematica
  All entries match! ✓

Test 6: Simple Current Algebra
  Jacobi identity holds for J(z)J(w) = k/(z-w)^2 ✓

All tests passed! ✓
```

### pytest 输出

```
tests/test_jacobi_virasoro.py::TestJacobiIdentityVirasoro::test_virasoro_ope_structure PASSED
tests/test_jacobi_virasoro.py::TestJacobiIdentityVirasoro::test_jacobi_identity_TTT PASSED
tests/test_jacobi_virasoro.py::TestJacobiIdentityVirasoro::test_verify_jacobi_identity_TTT PASSED
tests/test_jacobi_virasoro.py::TestJacobiIdentityVirasoro::test_jacobi_identity_dimensions PASSED
tests/test_jacobi_virasoro.py::TestJacobiIdentityVirasoro::test_comparison_with_mathematica PASSED
tests/test_jacobi_virasoro.py::TestJacobiIdentityProperties::test_jacobi_identity_for_simple_operators PASSED

6 passed, 4 warnings in 0.25s
```

## 测试覆盖

### 1. OPE 结构验证
- ✓ 验证 Virasoro OPE 的各阶极点
- ✓ 确认 T(z)T(w) 的正确形式

### 2. Jacobi 恒等式计算
- ✓ 实现 `check_jacobi_identity(A, B, C)` 函数
- ✓ 计算 Jacobi 恒等式的三项
- ✓ 验证结果为全零矩阵

### 3. 便捷函数
- ✓ 实现 `verify_jacobi_identity(A, B, C)` 函数
- ✓ 返回布尔值表示恒等式是否成立

### 4. 维度验证
- ✓ 验证结果矩阵的维度（5x5）
- ✓ 与 Mathematica 参考结果一致

### 5. 对比验证
- ✓ pyope 结果与 Mathematica 完全一致
- ✓ 所有矩阵元素都为 0

### 6. 扩展测试
- ✓ 测试简单流代数（current algebra）
- ✓ 验证 Jacobi 恒等式的普遍性

## 实现细节

### check_jacobi_identity 函数

```python
def check_jacobi_identity(A: Any, B: Any, C: Any, simplify_func=None) -> List[List[Any]]:
    """
    检查三个算符的 Jacobi 恒等式

    计算：
    [A, [B, C]_q]_m - sign * [B, [A, C]_m]_q - Σ_p C(n-1, p-1) [[A,B]_p, C]_{m+n-p}

    其中 sign = (-1)^(|A||B|)
    """
```

### 关键步骤

1. **计算基本 OPE**
   - OPE(A, B)
   - OPE(B, C)
   - OPE(A, C)

2. **计算嵌套 OPE**
   - AnBC[n] = OPE[A, {BC}_n]
   - BnAC[n] = OPE[B, {AC}_n]
   - ABnC[n] = OPE[{AB}_n, C]

3. **构建 Jacobi 矩阵**
   - 对每个 (m, n) 计算三项之和
   - 简化结果
   - 验证是否为 0

## 数学意义

Jacobi 恒等式的成立验证了：

1. **代数结构的一致性**: OPE 计算满足基本的代数恒等式
2. **正规序的正确性**: 正规序算符的 OPE 计算正确
3. **导数规则的正确性**: 导数算符的 OPE 规则正确
4. **对称性的正确性**: 算符交换的对称性公式正确

## 参考文献

1. **OPEdefs.m**: Kris Thielemans, "A Mathematica package for computing operator product expansions"
   - 第 1601-1637 行：OPEJacobi 实现

2. **共形场论**: P. Di Francesco, P. Mathieu, D. Sénéchal, "Conformal Field Theory"
   - Chapter 5: Operator Product Expansion
   - Chapter 6: Virasoro Algebra

3. **顶点算符代数**: V. Kac, "Vertex Algebras for Beginners"
   - Chapter 3: Jacobi Identity

## 结论

本测试成功验证了 pyope 库中 Jacobi 恒等式的实现：

1. ✓ Python 实现与 Mathematica 参考完全一致
2. ✓ Virasoro 代数满足 Jacobi 恒等式
3. ✓ 简单流代数满足 Jacobi 恒等式
4. ✓ 所有测试用例通过

这表明 pyope 库的 OPE 计算引擎在数学上是正确的，可以用于更复杂的顶点算符代数计算。

## 未来工作

1. 测试更多代数结构（如 W 代数、超对称代数）
2. 测试费米算符的 Jacobi 恒等式
3. 测试混合玻色-费米算符的情况
4. 性能优化和大规模计算测试
