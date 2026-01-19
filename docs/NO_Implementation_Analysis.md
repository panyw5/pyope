# OPEdefs.m 正规排序实现分析与 Python 移植方案

## 执行摘要

本文档分析了 Mathematica 包 OPEdefs.m 中正规排序（Normal Ordering）的实现机制，并提出了 Python 版本的实现方案。

**核心发现**：
1. OPEdefs.m 使用**模式匹配规则**递归处理嵌套 NO
2. 多参数 NO 通过右结合转换：`NO[A,B,C]` → `NO[A,NO[B,C]]`
3. 嵌套 NO 通过 **Jacobi 恒等式**展开，而非简单的 flatten
4. 算符排序基于声明顺序和导数阶数

---

## 1. OPEdefs.m 中的 NO 实现

### 1.1 多参数 NO 的处理（Line 1437）

```mathematica
Literal[NO[A_,B_,C__]] := NO[A,NO[B,C]]
```

**关键点**：
- 使用**右结合**：`NO[A,B,C]` 自动转换为 `NO[A,NO[B,C]]`
- 这是递归定义的基础，所有多参数情况最终归约为二元 NO

**Python 对应**：
```python
def NO(*args):
    if len(args) > 2:
        # 右结合：NO(A,B,C) -> NO(A, NO(B,C))
        return NO(args[0], NO(*args[1:]))
    elif len(args) == 2:
        return NormalOrderedOperator(args[0], args[1])
```

### 1.2 嵌套 NO 的展开规则

#### 规则 1：`NO[NO[A,B], C]` 的处理（Line 1467-1471）

```mathematica
Literal[NO[NO[A_,B_],C_]] :=
    If[NOOrder[A,C]>0,
       CallAndSave[NOCompositeHelpL,A,B,C],       (* -> {A,{B,C}}*)
       SwapSign[C,NO[A,B]] (NO[C,NO[A,B]] - NOCommuteHelp[C,NO[A,B]])
    ]
```

**逻辑**：
- 如果 `A` 应该排在 `C` 之后（`NOOrder[A,C]>0`），调用 `NOCompositeHelpL`
- 否则，使用对易关系将 `C` 移到左边

#### 规则 2：`NO[B, NO[A,C]]` 的处理（Line 1473-1475）

```mathematica
Literal[NO[B_,NO[A_,C_]]] :=
    CallAndSave[NOCompositeHelpR,B,A,C] /;
    NOOrder[A,B]>0
```

**逻辑**：
- 如果 `A` 应该排在 `B` 之后，调用 `NOCompositeHelpR`

#### 规则 3：`NO[NO[A,B], NO[C,D]]` 的处理（Line 1461-1465）

```mathematica
Literal[NO[NO[A_,B_],NO[C_,D_]]] :=
    If[NOOrder[A,C]<=0,
       CallAndSave[NOCompositeHelpR,NO[A,B],C,D] (* -> {C,{{A,B},D}}*),
       CallAndSave[NOCompositeHelpL,A,B,NO[C,D]] (* -> {A,{B,{C,D}}}*)
    ]
```

### 1.3 核心展开函数：`NOCompositeHelpLQ`（Line 1553-1569）

这是量子情况下 `NO[NO[A,B],C]` 的展开公式：

```mathematica
NOCompositeHelpLQ[A_,B_,C_] :=
    Block[{AC, BC, l,sign = SwapSign[A,B],
             maxAC, maxBC, maxl},
        AC = OPE[A,C];
        BC = If[ SameQ[A,B], AC, OPE[B,C]];
        maxAC = MaxPole[AC]; maxBC = MaxPole[BC];
        PoleSimplify[
            NO[A,NO[B,C]] +                                    (* 第1项 *)
            Sum[ NO[Derivative[l][A], opepole[l][BC]] /l!,    (* 第2项 *)
                {l,1, maxBC}
            ] +
            Sum[ sign * NO[Derivative[l][B], opepole[l][AC]] /l!,  (* 第3项 *)
                {l,1, maxAC}
            ]
        ]
    ]
```

**公式解释**（对应 Thielemans 论文 eq 3.3.4）：

$$
(AB)C = A(BC) + \sum_{l=1}^{\infty} \frac{1}{l!} (\partial^l A \{BC\}_l) + (-1)^{|A||B|} \sum_{l=1}^{\infty} \frac{1}{l!} (\partial^l B \{AC\}_l)
$$

**三项含义**：
1. **第1项**：`NO[A,NO[B,C]]` - 重新结合
2. **第2项**：A 与 C 的收缩通过 B 传递
3. **第3项**：B 与 C 的收缩通过 A 传递（带符号因子）

### 1.4 算符排序规则（Line 1424-1431）

```mathematica
NOOrder[A_,B_] :=
    Block[{res = NOOrderHelp[A,B]},
        If[res == 0,
           NOOrderingValue * (NOOrderHelp2[B]-NOOrderHelp2[A]),
           res]
    ]
NOOrderHelp2[Derivative[i_][_]] := i
NOOrderHelp2[_] = 0;
```

**排序逻辑**：
1. 首先按 `OPEOrder`（声明顺序）排序
2. 如果基础算符相同，按导数阶数排序
3. `NOOrderingValue` 控制导数排序方向：
   - `< 0`：高阶导数在左（默认）
   - `> 0`：低阶导数在左

---

## 2. Python 当前实现的问题

### 2.1 当前代码（`src/pyope/simplify.py:224-234`）

```python
# 处理嵌套的 NO
# 如果左侧是 NO，展开为 NO(NO(A,B), C)
if isinstance(left, NormalOrderedOperator):
    # 这需要使用 OPE 来展开
    # 暂时保持原样，因为完整展开需要 OPE 信息
    pass

# 如果右侧是 NO，展开为 NO(A, NO(B,C))
if isinstance(right, NormalOrderedOperator):
    # 同样暂时保持原样
    pass
```

**问题**：
- 明确标记为 TODO，未实现任何逻辑
- 导致 `simplify(NO(NO(J_zero, J_minus), J_plus))` 原样返回

### 2.2 触发条件限制（Line 239-241）

```python
if isinstance(left, (BasisOperator, DerivativeOperator)) and isinstance(
    right, (BasisOperator, DerivativeOperator)
):
```

**问题**：
- 只处理两个基本算符的交换
- 遇到 `NormalOrderedOperator` 就跳过
- 无法触发嵌套 NO 的展开

---

## 3. 方案 A 的实现建议

### 3.1 总体策略

**不采用 flatten 方式**，而是**递归应用 Jacobi 恒等式**：

```
NO(NO(A,B),C) → 使用 NOCompositeHelpL 展开
              → A(BC) + 收缩项
              → 递归 simplify 每一项
```

### 3.2 实现步骤

#### 步骤 1：在 `_simplify_normal_ordered` 中添加嵌套处理

```python
def _simplify_normal_ordered(left, right, expand_derivatives=True):
    # ... 现有的线性性、标量提取等逻辑 ...
    
    # 处理嵌套 NO：NO(NO(A,B), C)
    if isinstance(left, NormalOrderedOperator):
        return _expand_nested_no_left(left.left, left.right, right)
    
    # 处理嵌套 NO：NO(A, NO(B,C))
    if isinstance(right, NormalOrderedOperator):
        return _expand_nested_no_right(left, right.left, right.right)
    
    # ... 现有的交换逻辑 ...
```

#### 步骤 2：实现 `_expand_nested_no_left`（对应 `NOCompositeHelpLQ`）

```python
def _expand_nested_no_left(A, B, C):
    """
    展开 NO(NO(A,B), C) 使用 Jacobi 恒等式
    
    公式：(AB)C = A(BC) + Σ (1/l!) (∂^l A {BC}_l) + sign * Σ (1/l!) (∂^l B {AC}_l)
    """
    from .api import OPE
    from .operators import d as derivative_operator
    from .local_operator import get_operator_parity
    
    # 1. 计算符号因子
    parity_A = get_operator_parity(A)
    parity_B = get_operator_parity(B)
    sign = (-1) ** (parity_A * parity_B)
    
    # 2. 计算 OPE
    ope_AC = OPE(A, C)
    ope_BC = OPE(B, C) if A != B else ope_AC
    
    # 3. 第1项：A(BC)
    result = NO(A, NO(B, C))
    
    # 4. 第2项：Σ (1/l!) (∂^l A {BC}_l)
    for l in range(1, ope_BC.max_pole + 1):
        pole_l = ope_BC.pole(l)
        if pole_l != 0:
            deriv_A = derivative_operator(A, l)
            coeff = 1 / sp.factorial(l)
            result = result + coeff * NO(deriv_A, pole_l)
    
    # 5. 第3项：sign * Σ (1/l!) (∂^l B {AC}_l)
    for l in range(1, ope_AC.max_pole + 1):
        pole_l = ope_AC.pole(l)
        if pole_l != 0:
            deriv_B = derivative_operator(B, l)
            coeff = sign / sp.factorial(l)
            result = result + coeff * NO(deriv_B, pole_l)
    
    # 6. 递归化简
    return simplify(result)
```

#### 步骤 3：实现 `_expand_nested_no_right`（对应 `NOCompositeHelpRQ`）

```python
def _expand_nested_no_right(B, A, C):
    """
    展开 NO(B, NO(A,C)) 使用 Jacobi 恒等式
    
    公式：B(AC) = (-1)^{|A||B|} A(BC) + (BA-(-1)^{|A||B|}AB)C
    """
    from .api import OPE
    from .local_operator import get_operator_parity
    
    # 1. 计算符号因子
    parity_A = get_operator_parity(A)
    parity_B = get_operator_parity(B)
    sign = (-1) ** (parity_A * parity_B)
    
    # 2. 第1项：sign * A(BC)
    result = sign * NO(A, NO(B, C))
    
    # 3. 第2项：NOCommuteHelp[B,A] C
    # NOCommuteHelp[B,A] = NO[B,A] - sign * NO[A,B]
    commute_term = _compute_no_commute_help(B, A)
    result = result + NO(commute_term, C)
    
    # 4. 递归化简
    return simplify(result)
```

#### 步骤 4：实现 `_compute_no_commute_help`（对应 `NOCommuteHelpQ`）

```python
def _compute_no_commute_help(A, B):
    """
    计算 NOCommuteHelp[A,B] = NO[A,B] - sign * NO[B,A]
    
    使用公式：Σ (-1)^m / m! ∂^m {AB}_m
    """
    from .api import OPE
    from .operators import d as derivative_operator
    from .constants import Zero
    
    ope_AB = OPE(A, B)
    result = Zero
    
    for m in range(1, ope_AB.max_pole + 1):
        pole_m = ope_AB.pole(m)
        if pole_m != 0:
            deriv_pole = derivative_operator(pole_m, m)
            coeff = (-1) ** m / sp.factorial(m)
            result = result + coeff * deriv_pole
    
    return simplify(result)
```

### 3.3 处理算符排序

在展开之前，先检查是否需要重新排序：

```python
def _simplify_normal_ordered(left, right, expand_derivatives=True):
    # ... 前面的逻辑 ...
    
    # 处理嵌套 NO：NO(NO(A,B), C)
    if isinstance(left, NormalOrderedOperator):
        A, B = left.left, left.right
        # 检查是否需要重新结合
        order_AC = ope_registry.compare_operators(A, C)
        if order_AC > 0:
            # A 应该在 C 之后，需要展开
            return _expand_nested_no_left(A, B, C)
        else:
            # A 应该在 C 之前或相同，使用对易关系
            # NO(NO(A,B),C) = sign * (NO(C,NO(A,B)) - NOCommuteHelp(C,NO(A,B)))
            sign = _compute_swap_sign(C, left)
            commute_term = _compute_no_commute_help(C, left)
            return sign * (NO(C, left) - commute_term)
    
    # ... 类似处理 NO(A, NO(B,C)) ...
```

---

## 4. 关键设计决策

### 4.1 为什么不用 flatten？

**OPEdefs.m 的设计哲学**：
- NO 是**二元运算**，不是 n-ary 列表
- 嵌套 NO 的展开**不是简单的重新排列**，而是涉及 OPE 收缩
- 使用 Jacobi 恒等式保证了**代数一致性**

**如果强行 flatten**：
```python
# 错误的做法
NO(NO(A,B),C) → [A,B,C] → 排序 → NO(NO(A,C),B)
```
这会**丢失 OPE 收缩项**，导致结果错误。

### 4.2 递归终止条件

递归会在以下情况终止：
1. 所有 NO 都是二元且已排序
2. 所有嵌套 NO 都被展开
3. 所有 OPE 收缩项都被计算

### 4.3 性能优化

OPEdefs.m 使用 `CallAndSave` 缓存中间结果：

```mathematica
CallAndSave[NOCompositeHelpL,A,B,C]
```

Python 可以使用 `functools.lru_cache`：

```python
@lru_cache(maxsize=1024)
def _expand_nested_no_left_cached(A, B, C):
    return _expand_nested_no_left(A, B, C)
```

---

## 5. 测试策略

### 5.1 单元测试

```python
def test_nested_no_expansion():
    """测试嵌套 NO 的展开"""
    # 定义算符
    A = BasisOperator("A", conformal_weight=1)
    B = BasisOperator("B", conformal_weight=1)
    C = BasisOperator("C", conformal_weight=1)
    
    # 定义 OPE
    OPE[A,B] = MakeOPE([One, A])  # {AB}_2 = 1, {AB}_1 = A
    OPE[A,C] = MakeOPE([B])       # {AC}_1 = B
    OPE[B,C] = MakeOPE([One])     # {BC}_1 = 1
    
    # 测试展开
    result = simplify(NO(NO(A,B),C))
    
    # 验证结果不再是嵌套形式
    assert not isinstance(result, NormalOrderedOperator) or \
           not isinstance(result.left, NormalOrderedOperator)
```

### 5.2 与 Mathematica 对比

```python
def test_compare_with_mathematica():
    """与 OPEdefs.m 结果对比"""
    # 使用 Kac-Moody 代数
    J_plus = BasisOperator("J⁺", conformal_weight=1)
    J_zero = BasisOperator("J⁰", conformal_weight=1)
    J_minus = BasisOperator("J⁻", conformal_weight=1)
    
    # 定义 OPE（与 Mathematica 相同）
    k = 1
    OPE[J_plus, J_zero] = MakeOPE([-2 * J_plus])
    OPE[J_plus, J_minus] = MakeOPE([k * One, J_zero])
    OPE[J_zero, J_minus] = MakeOPE([-2 * J_minus])
    OPE[J_zero, J_zero] = MakeOPE([2 * k * One, 0])
    
    # 计算
    result = simplify(NO(NO(J_zero, J_minus), J_plus))
    
    # 在 Mathematica 中运行相同计算，对比结果
    # （需要手动验证或使用 wolframclient）
```

---

## 6. 实现优先级

### 阶段 1：基础实现（1-2天）
- [ ] 实现 `_expand_nested_no_left`
- [ ] 实现 `_compute_no_commute_help`
- [ ] 修改 `_simplify_normal_ordered` 调用新函数
- [ ] 基本单元测试

### 阶段 2：完整实现（2-3天）
- [ ] 实现 `_expand_nested_no_right`
- [ ] 处理 `NO(NO[A,B], NO[C,D])` 情况
- [ ] 添加算符排序检查
- [ ] 完整测试套件

### 阶段 3：优化与验证（1-2天）
- [ ] 添加缓存机制
- [ ] 性能测试
- [ ] 与 Mathematica 对比验证
- [ ] 文档完善

---

## 7. 潜在问题与解决方案

### 问题 1：无限递归

**场景**：`simplify(NO(A,B))` 调用 `simplify(NO(A,B))`

**解决**：
- 在递归前检查是否已经是最简形式
- 使用 `visited` 集合追踪已处理的表达式

### 问题 2：表达式膨胀

**场景**：展开后产生大量项

**解决**：
- 每步都调用 `collect_normal_ordered_terms` 合并同类项
- 使用 `expand_derivatives=False` 选项延迟导数展开

### 问题 3：符号因子错误

**场景**：费米子算符的符号计算错误

**解决**：
- 严格按照 Thielemans 公式实现
- 添加专门的符号因子测试

---

## 8. 参考资料

1. **OPEdefs.m 源码**：
   - Line 1437: 多参数 NO 定义
   - Line 1461-1491: 嵌套 NO 规则
   - Line 1553-1569: `NOCompositeHelpLQ` 实现

2. **Thielemans 论文**：
   - Section 2.3: 正规排序定义
   - Equation 2.3.16: 交换公式
   - Equation 3.3.4: Jacobi 恒等式

3. **VOA Manual**：
   - Section 3.3: Implementation details

---

## 9. 总结

**核心洞察**：
1. 嵌套 NO 的展开**不是语法重写**，而是**代数计算**
2. 必须使用 **Jacobi 恒等式**保证一致性
3. 递归 + 缓存是处理复杂嵌套的关键

**实现路径**：
- 不要试图 flatten，而是递归应用 Jacobi 恒等式
- 每一步都化简，避免表达式膨胀
- 充分测试，与 Mathematica 对比验证

**预期效果**：
```python
simplify(NO(NO(J_zero, J_minus), J_plus))
# 输出：展开后的线性组合，不再有嵌套 NO
```
