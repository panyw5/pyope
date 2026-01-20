# Bug 修复：NO(∂J⁰, J⁰) 错误简化为 Zero

## 问题描述

在 Kac-Moody 代数的测试中发现，`NO(∂J⁰, J⁰)` 被错误地简化为 `Zero`。

```python
J_zero = BasisOperator("J⁰", conformal_weight=1)
OPE[J_zero, J_zero] = OPEData({2: 2 * k_val * One})

expr = NO(d(J_zero), J_zero)
result = simplify(expr)
# 错误结果: Zero
```

## 为什么这是错的

在 VOA 理论中，`NO(∂J⁰, J⁰)` 对应态 `J⁰_{-2} J⁰_{-2}|0>`，这是一个**普通态**（非零态）。

正规序 `NO(∂J⁰, J⁰)` 不应该简化为零，除非有特殊的 OPE 关系导致约简。

## Bug 根源

问题出在 `src/pyope/simplify.py` 的第 471-531 行，当需要交换 `NO(B, A)` 为 `NO(A, B)` 时，代码错误地使用了 `ope_AB.pole(0)` 来获取 Thielemans 公式 (2.3.16) 的 l=0 项：

```python
# 错误的实现
result = ZeroConst

# l=0 项: [AB]_0
pole_0 = ope_AB.pole(0)  # ❌ 这总是返回 Zero！
if pole_0 != 0 and pole_0 != ZeroConst:
    result = swap_sign * pole_0
```

**问题**：
1. `OPE(A, B)` 不包含 0 阶极点（0 阶由 `NO(A, B)` 定义）
2. `ope_AB.pole(0)` 总是返回 `Zero`
3. 导致 Thielemans 公式 (2.3.16) 的主项 `NO(A, B)` 被丢失

## 修复方案

修改 `src/pyope/simplify.py`，直接使用 `NO(A, B)` 作为 l=0 项：

```python
# 正确的实现
# FIX: l=0 项是 NO(A,B)，不是 ope_AB.pole(0)
# 因为 OPE(A,B) 只包含奇异部分（l>=1），0 阶由 NO 定义
result = swap_sign * NO(right, left)

# l >= 1 项: \frac{(-1)^l}{l!} \partial^l [AB]_l
for l in range(1, max_pole + 1):
    pole_l = ope_AB.pole(l)
    if pole_l != 0 and pole_l != ZeroConst:
        deriv_pole = derivative_operator(pole_l, l)
        coeff = swap_sign * ((-1) ** l) / sp.factorial(l)
        result = result + coeff * deriv_pole
```

## 额外修复：导数排序规则

Oracle 指出当前的排序倾向于把带导数的放右边，但通常应该是**导数越多越靠左**。

修改 `src/pyope/registry.py` 的 `compare_operators` 方法：

```python
# 修改前（错误）
# 基础算符相同，比较导数阶数
# 阶数小的在前
if left_order != right_order:
    return right_order - left_order  # ❌

# 修改后（正确）
# 基础算符相同，比较导数阶数
# 导数越多越靠左（阶数大的在前）
if left_order != right_order:
    return left_order - right_order  # ✅
```

## 测试结果

### 修复前

```python
NO(∂J⁰, J⁰) -> Zero  # ❌ 错误
NO(J⁰, ∂J⁰) -> NO(∂J⁰, J⁰)  # ❌ 排序错误
```

### 修复后

```python
NO(∂J⁰, J⁰) -> NO(∂J⁰, J⁰)  # ✅ 正确（保持不变）
NO(J⁰, ∂J⁰) -> NO(∂J⁰, J⁰)  # ✅ 正确（重排）
```

### Kac-Moody 一致性测试

```python
expr = NO(J_minus, NO(J_plus, J_zero))

# 修复前
simplify(expr) = 2*NO(J⁺,∂J⁻)  # ❌ 丢失了项
expand_nested_no(expr) = NO(J⁺,NO(J⁰,J⁻)) + 2*NO(J⁺,∂J⁻) - NO(∂J⁰,J⁰)

# 修复后
simplify(expr) = NO(J⁺,NO(J⁰,J⁻)) + 2*NO(J⁺,∂J⁻) - NO(∂J⁰,J⁰)  # ✅
expand_nested_no(expr) = NO(J⁺,NO(J⁰,J⁻)) + 2*NO(J⁺,∂J⁻) - NO(∂J⁰,J⁰)  # ✅
# 两者相等！
```

## 新增测试

创建了 `tests/test_no_derivative_bug_fix.py`，包含3个测试：

1. ✅ `test_no_derivative_not_zero` - 验证 `NO(∂J⁰, J⁰)` 不简化为 Zero
2. ✅ `test_kac_moody_consistency` - 验证 `simplify` 和 `expand_nested_no` 的一致性
3. ✅ `test_derivative_ordering` - 验证导数排序规则

## 影响范围

- **修复的文件**：
  - `src/pyope/simplify.py` - 修复 l=0 项的获取
  - `src/pyope/registry.py` - 修复导数排序规则
  - `tests/test_simplify.py` - 更新测试以反映新的排序规则

- **测试结果**：
  - 所有现有测试通过 ✅
  - 新增3个测试全部通过 ✅

## 总结

这是一个严重的 bug，会导致某些正规序项被错误地消除。修复后：

1. ✅ `NO(∂J⁰, J⁰)` 不再错误地简化为 `Zero`
2. ✅ 导数排序规则符合 VOA 惯例（导数越多越靠左）
3. ✅ `simplify` 和 `expand_nested_no` 的结果一致
4. ✅ 所有测试通过

感谢 Oracle 的深入分析！
