# Simplify 函数正规排序功能整合总结

## 任务目标

将 `expand_nested_no` 的两种正规排序功能整合到 `simplify` 函数中：

1. **算符之间按标准顺序排列**：`NO(C, NO(A, B)) -> NO(A, NO(B, C)) + additional terms`
2. **NO 化成右嵌套**：`NO(NO(...), ...) -> NO(..., NO(..., NO(...)))`

## 实现方案

### 1. 增强右嵌套的排序检查

修改了 `_simplify_normal_ordered` 函数中处理右嵌套 `NO(B, NO(A, C))` 的逻辑：

- **检查外层顺序**：比较 `B` 和 `A` 的顺序
- **检查内层顺序**：比较 `A` 和 `C` 的顺序
- **递归简化**：确保内层的 NO 也被正确排序

关键修改：
```python
# 情况 1: B 和 A 需要交换（order_BA < 0），即 A 应该排在 B 前面
if order_BA < 0:
    expanded = _expand_nested_no_right(B, A, C, expand_derivatives)
    return simplify(expanded, expand_derivatives, preserve_nested_structure)

# 情况 2: B 和 A 的顺序正确，但 A 和 C 需要交换（order_AC < 0）
if order_AC < 0:
    simplified_right = simplify(right, expand_derivatives, preserve_nested_structure)
    if simplified_right != right:
        return simplify(NO(B, simplified_right), expand_derivatives, preserve_nested_structure)

# 情况 3: 所有顺序都正确，但仍需递归简化内层
simplified_right = simplify(right, expand_derivatives, preserve_nested_structure)
if simplified_right != right:
    return simplify(NO(B, simplified_right), expand_derivatives, preserve_nested_structure)
```

### 2. 修复左嵌套转换的递归问题

修改了左嵌套转换逻辑，确保转换后的结果被递归简化：

```python
# 没有非零 OPE，但仍需要转换为右嵌套形式
# NO(NO(A,B), C) -> NO(A, NO(B,C))（带符号因子）
parity_AB = (get_operator_parity(A) + get_operator_parity(B)) % 2
parity_C = get_operator_parity(C)
sign = (-1) ** (parity_AB * parity_C)
# 递归简化转换后的结果，因为 B 或 C 可能仍然是嵌套的 NO
result = sign * NO(A, NO(B, C))
return simplify(result, expand_derivatives, preserve_nested_structure)
```

这个修复解决了三重嵌套的问题：`NO(NO(NO(A,B),C),D)` 现在能正确转换为 `NO(A,NO(B,NO(C,D)))`。

### 3. 修正 compare_operators 返回值的理解

发现并修正了对 `compare_operators` 返回值语义的理解：
- `> 0`: 顺序正确（第一个算符应该在前）
- `= 0`: 两个算符相等
- `< 0`: 需要交换（第一个算符应该在后）

## 测试结果

### 新增测试

创建了 `tests/test_simplify_canonical_ordering.py`，包含9个测试用例：

1. ✅ `test_right_nested_ordering_simple` - 右嵌套简单排序
2. ✅ `test_right_nested_ordering_with_ope` - 右嵌套排序（有非零 OPE）
3. ✅ `test_left_nested_to_right_nested` - 左嵌套转换为右嵌套
4. ✅ `test_left_nested_to_right_nested_with_ope` - 左嵌套转换（有非零 OPE）
5. ✅ `test_triple_nested_ordering` - 三重嵌套排序
6. ✅ `test_preserve_nested_structure_flag` - preserve_nested_structure 参数测试
7. ✅ `test_ordering_consistency` - 排序一致性测试
8. ✅ `test_inner_layer_needs_reordering` - 内层需要重排序
9. ✅ `test_both_layers_need_reordering` - 内外层都需要重排序

**所有9个测试全部通过！**

### 完整测试套件

运行完整测试套件的结果：
- **292 passed** ✅
- **1 failed** ⚠️ (与本次修改无关，原分支中也失败)

失败的测试：
- `test_betagamma_stress_tensor_definition_lambda_3_2` - beta-gamma 系统的符号问题（原分支中也失败）

### 修复的测试

修改了 `tests/test_simplify.py` 中的 `test_expand_nested_no` 测试，以反映新的默认行为（左嵌套自动转换为右嵌套）。

## 功能验证

### 示例 1: 右嵌套排序

```python
A = BasisOperator("A", conformal_weight=1)
B = BasisOperator("B", conformal_weight=1)
C = BasisOperator("C", conformal_weight=1)

expr = NO(C, NO(A, B))
result = simplify(expr)
# 结果: NO(A, NO(B, C))
```

### 示例 2: 左嵌套转换

```python
expr = NO(NO(A, B), C)
result = simplify(expr)
# 结果: NO(A, NO(B, C))
```

### 示例 3: 三重嵌套

```python
D = BasisOperator("D", conformal_weight=1)
expr = NO(NO(NO(A, B), C), D)
result = simplify(expr)
# 结果: NO(A, NO(B, NO(C, D)))
```

### 示例 4: 内层排序

```python
expr = NO(A, NO(C, B))  # 内层 C > B（需要重排）
result = simplify(expr)
# 结果: NO(A, NO(B, C))
```

## 文档更新

更新了 `simplify` 函数的 docstring，明确说明了两种正规排序功能：

```python
def simplify(
    expr: Any, expand_derivatives: bool = True, preserve_nested_structure: bool = False
) -> Any:
    """
    化简 OPE 表达式为规范形式

    将表达式化简为排序的 NO product 的线性组合。执行以下操作：
    1. 展开嵌套的 NO 乘积（默认将左嵌套转换为右嵌套）
    2. 在 NO 内部按算符顺序排列（实现两种正规排序）
       - 算符之间按标准顺序排列：NO(C, NO(A, B)) -> NO(A, NO(B, C)) + additional terms
       - NO 化成右嵌套：NO(NO(...), ...) -> NO(..., NO(..., NO(...)))
    3. 合并同类项
    4. 标准化导数表示
    5. （可选）应用莱布尼茨法则展开正规序的导数
    ...
    """
```

## 关键改进

1. **完整的递归简化**：确保所有嵌套层级都被正确处理
2. **正确的顺序检查**：修正了 `compare_operators` 返回值的理解
3. **三重嵌套支持**：修复了三重嵌套无法完全展开的问题
4. **向后兼容**：通过 `preserve_nested_structure` 参数保持旧行为的可选性

## 总结

成功将 `expand_nested_no` 的两种正规排序功能整合到 `simplify` 函数中，实现了：

✅ 算符之间按标准顺序排列
✅ NO 化成右嵌套
✅ 递归处理多层嵌套
✅ 完整的测试覆盖
✅ 向后兼容性

所有新功能测试通过，完整测试套件保持高通过率（292/293），唯一失败的测试与本次修改无关。
