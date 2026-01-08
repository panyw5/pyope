# N=1 超共形代数 OPE 公式修正

## 问题发现

在测试 N=1 超共形代数的 Jacobi 恒等式时，发现 `OPEJacobi[T, G, G]` 测试失败，出现非零项 `c/2 * One`。

## 根本原因

通过查阅 Thielemans 论文 "An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings"（第 2.4.13 节），发现我们使用的 G(z)G(w) OPE 公式中的中心荷系数不正确。

### 错误的公式

```
G(z)G(w) = (c/3)/(z-w)^3 + 2T(w)/(z-w)
```

### 正确的公式（来自论文 eq. 2.4.13）

```
G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)
```

**关键差异**：中心荷项的系数应该是 `2c/3` 而不是 `c/3`，差了一个因子 2。

## 论文引用

来自 Thielemans 论文第 2.4.13 节：

> As an example, the N = 1 superconformal algebra consists of a Virasoro operator T and a fermionic dimension 3/2 primary operator G with OPE:
>
> $$G \times G = \ll \frac{2c}{3} \mathbb{1} | 0 | 2T \gg$$

这里使用的是 OPEdefs.m 的简写记法，展开后就是：

```
G(z)G(w) = (2c/3)/(z-w)^3 + 0/(z-w)^2 + 2T(w)/(z-w)
```

## 修正内容

### 1. Python 测试文件

文件：`tests/test_jacobi_superconformal.py`

修改：
```python
# 修改前
OPE[self.G, self.G] = OPE.make([self.c/3 * One, 0, 2*self.T])

# 修改后
OPE[self.G, self.G] = OPE.make([2*self.c/3 * One, 0, 2*self.T])
```

### 2. Mathematica 参考测试

文件：`tests/ref_jacobi_superconformal.wls`

修改：
```mathematica
(* 修改前 *)
OPE[G, G] = OPEData[{0, c/3 One, 2 T}];

(* 修改后 *)
OPE[G, G] = OPEData[{2*c/3 One, 0, 2 T}];
```

注意：Mathematica 的 OPEData 列表顺序是从高阶极点到低阶极点，所以 `{2*c/3 One, 0, 2 T}` 对应：
- 3 阶极点：`2*c/3 One`
- 2 阶极点：`0`
- 1 阶极点：`2 T`

## 测试结果

修正后，所有 Jacobi 恒等式测试都通过：

### Python (pyope) 测试结果

```
✓ G-G-G: PASS (fermionic sign factor works correctly)
✓ T-G-G: PASS (corrected OPE formula with 2c/3)
✓ T-T-G: PASS
```

### Mathematica (OPEdefs.m) 测试结果

```
OPEJacobi[G, G, G]: PASS
OPEJacobi[T, G, G]: PASS
OPEJacobi[T, T, G]: PASS

ALL TESTS PASSED!
```

### pytest 结果

```
8 passed, 14 warnings in 0.14s
```

所有 8 个测试用例全部通过。

## 完整的 N=1 超共形代数 OPE

修正后的完整 OPE 定义：

```
1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)

2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)

3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)
```

其中：
- T: 能动张量，共形权重 h=2，玻色子（parity=0）
- G: 超流，共形权重 h=3/2，费米子（parity=1）
- c: 中心荷

## 验证方法

1. 运行 Mathematica 参考测试：
   ```bash
   cd tests
   wolframscript ref_jacobi_superconformal.wls
   ```

2. 运行 Python 测试：
   ```bash
   python tests/test_jacobi_superconformal.py
   ```

3. 运行 pytest：
   ```bash
   pytest tests/test_jacobi_superconformal.py -v
   ```

## 结论

通过查阅原始论文并修正 OPE 公式，我们成功解决了 Jacobi 恒等式测试失败的问题。这证明：

1. **pyope 的实现是正确的**：费米子符号处理、OPE 计算、Jacobi 恒等式检查都正确实现
2. **问题在于输入数据**：使用了错误的 OPE 公式
3. **文献查证的重要性**：在遇到测试失败时，应该首先查证理论公式是否正确

这次修正也为 pyope 库提供了一个重要的测试用例，验证了库在处理费米子算符时的正确性。

## 参考文献

Kris Thielemans, "An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings", PhD thesis KULeuven (1994), hep-th/9506159, Section 2.4, Equation 2.4.13.
