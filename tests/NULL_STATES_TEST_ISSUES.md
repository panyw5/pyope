# Null States 测试文件问题报告

## 发现的问题

### 1. test_null_states_stage1.py

#### 问题 1: 统计性声明错误（第 33 行）
```python
# 错误代码
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
Bosonic(b, c)  # ❌ 错误！b 和 c 是费米子（bosonic=False）
```

**修复**：
```python
Fermionic(b, c)  # ✓ 正确
```

#### 问题 2: 统计性声明错误（第 98 行）
```python
# 错误代码
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))
Fermionic(beta, gamma)  # ❌ 错误！beta 和 gamma 是玻色子（bosonic=True）
```

**修复**：
```python
Bosonic(beta, gamma)  # ✓ 正确
```

### 2. 缺少的测试

#### 缺少测试 A: 费米子自身 NO 乘积
**问题**：没有测试 `NO(ψ, ψ) = 0` 的情况（当权重相同时）

**建议添加**：
```python
def test_fermionic_self_product():
    """测试费米子与自身的 NO 乘积为零"""
    b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
    Fermionic(b)

    # 定义 OPE（无极点）
    OPE[b, b] = OPE.make([])

    # 测试 NO(b, b) = 0
    result = simplify(NO(b, b))
    assert result == 0, f"Expected NO(b, b) = 0, got {result}"
```

#### 缺少测试 B: 不同权重的相同费米子
**问题**：没有测试 `NO(ψ, ∂ψ) ≠ 0` 的情况

**建议添加**：
```python
def test_fermionic_different_weights():
    """测试不同权重的相同费米子的 NO 乘积非零"""
    g = BasisOperator('g', bosonic=False, conformal_weight=Fraction(3, 2))
    Fermionic(g)

    # g 和 ∂g 权重不同
    dg = d(g)  # 权重 5/2

    # NO(g, ∂g) 应该非零
    result = simplify(NO(g, dg))
    assert result != 0, f"Expected NO(g, ∂g) ≠ 0"
```

#### 缺少测试 C: 算符排序
**问题**：没有测试算符在 NO 中的排序规则

**建议添加**：
```python
def test_operator_ordering():
    """测试算符排序规则"""
    b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
    c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
    Fermionic(b, c)

    OPE[b, c] = OPE.make([One])

    # 测试 NO(b, c) vs NO(c, b)
    no_bc = simplify(NO(b, c))
    no_cb = simplify(NO(c, b))

    # 对于费米子：NO(c, b) = -NO(b, c) + {cb}_{≥1}
    # 如果 {cb}_{≥1} = 0，则应该满足反对称性
    print(f"NO(b, c) = {no_bc}")
    print(f"NO(c, b) = {no_cb}")
```

#### 缺少测试 D: 导数的莱布尼兹律
**问题**：没有测试导数运算的正确性

**建议添加**：
```python
def test_derivative_leibniz():
    """测试导数的莱布尼兹律"""
    b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
    c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
    Fermionic(b, c)

    # ∂(NO(b, c)) = NO(∂b, c) + NO(b, ∂c)
    lhs = d(NO(b, c))
    rhs = NO(d(b), c) + NO(b, d(c))

    diff = simplify(lhs - rhs)
    assert diff == 0, f"Leibniz rule failed: {diff}"
```

### 3. 其他需要检查的文件

需要检查以下文件是否有类似问题：
- [ ] test_null_states_stage2.py
- [ ] test_null_states_verification.py
- [ ] test_null_states_grouped.py
- [ ] test_operator_sorting.py
- [ ] test_partition_enumerator.py
- [ ] compare_enumeration_methods.py

## 优先级

1. **高优先级**：修复统计性声明错误（会导致计算错误）
2. **中优先级**：添加费米子自身 NO 乘积测试
3. **低优先级**：添加其他微妙情况的测试

## 建议的修复顺序

1. 立即修复 test_null_states_stage1.py 中的两个统计性错误
2. 添加费米子自身 NO 乘积的测试
3. 检查其他测试文件是否有类似问题
4. 逐步添加更全面的测试覆盖
