# Bug 修复报告：OPE q=0 情况的正确处理

## 问题描述

在 `src/pyope/api.py` 的 `_ope_composite_left` 函数中，当 `max_AC = max_BC = 0`（即两个 OPE 都没有奇异部分）时，代码错误地返回了空的 `OPEData({})`，而不是正确的 q=0 正规序乘积。

### 原始代码（第 584-586 行）

```python
# 如果两个 OPE 都是零，直接返回零
if max_AC == 0 and max_BC == 0:
    return OPEData({})  # ❌ 错误
```

## 理论依据

根据 Thielemans 论文 "An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings"：

### OPE 的完整定义（eq. 2.3.3）

$$T_1(z)T_2(w) = \sum_{n \leq h(T_1,T_2)} \frac{[T_1T_2]_n(w)}{(z-w)^n}$$

其中：
- **$n \geq 1$**: 奇异部分（singular part）
- **$n = 0$**: 正规序乘积 $[T_1T_2]_0 = :T_1T_2:$（eq. 2.3.7）
- **$n < 0$**: 正则部分（regular part）

### 关键公式

**eq. (3.3.18)** 适用于 $q \geq 1$（奇异部分）：
$$[[AB]_0 C]_q = \sum_{l \geq 0} \frac{1}{l!}[\partial^l A [BC]_{l+q}]_0 + (-1)^{|A||B|} \sum_{l \geq 0} \frac{1}{l!}[\partial^l B [AC]_{l+q}]_0 + (-1)^{|A||B|} \sum_{l=1}^{q-1} [B[AC]_{q-l}]_l$$

**eq. (3.3.19)** 适用于 $q = 0$（正规序）：
$$[[AB]_0 C]_0 = [A[BC]_0]_0 + \sum_{l > 0} \frac{1}{l!}[\partial^l A [BC]_l]_0 + (-1)^{|A||B|} \sum_{l > 0} \frac{1}{l!}[\partial^l B [AC]_l]_0$$

### 关键理解

当 $\max_{AC} = \max_{BC} = 0$ 时：
- 意味着 $[AC]_q = 0$ 和 $[BC]_q = 0$ 对所有 $q \geq 1$（**没有奇异部分**）
- 但根据 eq. (3.3.19)，此时：
  $$[[AB]_0 C]_0 = [A[BC]_0]_0 = [A \cdot NO(B,C)]_0 = NO(A, NO(B,C))$$
- 这**不是零**！而是一个正规序乘积

## 修复方案

### 修复后的代码（第 584-587 行）

```python
# 特殊情况：当 max_AC = max_BC = 0 时，根据 Thielemans eq. (3.3.19)
# [[AB]_0 C]_0 = [A[BC]_0]_0 = NO(A, NO(B,C))
# OPE 没有奇异部分（q >= 1），但有 q=0 的正规序乘积
if max_AC == 0 and max_BC == 0:
    return OPEData({0: NormalOrderedOperator(left, right)})
```

## 测试验证

创建了 `tests/test_q0_case_fix.py`，包含 5 个测试用例：

1. ✅ `test_no_singular_part_returns_normal_order` - 验证没有奇异部分时返回正规序
2. ✅ `test_bracket_q0_with_no_singular_part` - 验证 bracket(NO(A,B), C, 0) 的行为
3. ✅ `test_mixed_case_one_ope_has_poles` - 验证混合情况（一个有奇异部分，一个没有）
4. ✅ `test_both_opes_have_poles` - 验证两个都有奇异部分的情况（确保不破坏原有功能）
5. ✅ `test_theoretical_consistency_eq_3_3_19` - 验证理论一致性

### 测试结果

```
tests/test_q0_case_fix.py::TestQ0CaseFix::test_no_singular_part_returns_normal_order PASSED
tests/test_q0_case_fix.py::TestQ0CaseFix::test_bracket_q0_with_no_singular_part PASSED
tests/test_q0_case_fix.py::TestQ0CaseFix::test_mixed_case_one_ope_has_poles PASSED
tests/test_q0_case_fix.py::TestQ0CaseFix::test_both_opes_have_poles PASSED
tests/test_q0_case_fix.py::TestQ0CaseFix::test_theoretical_consistency_eq_3_3_19 PASSED

5 passed in 0.02s
```

### 完整测试套件结果

- **117 个测试通过** ✅
- **1 个测试失败**（与本次修复无关，是之前就存在的问题）

## 修复前后对比

### 修复前

```python
A = BasisOperator('A', conformal_weight=1)
B = BasisOperator('B', conformal_weight=1)
C = BasisOperator('C', conformal_weight=1)

OPE[A, C] = OPE.make([])  # 没有奇异部分
OPE[B, C] = OPE.make([])  # 没有奇异部分

result = OPE(NO(A, B), C)
print(result)  # 输出: Zero ❌
```

### 修复后

```python
A = BasisOperator('A', conformal_weight=1)
B = BasisOperator('B', conformal_weight=1)
C = BasisOperator('C', conformal_weight=1)

OPE[A, C] = OPE.make([])  # 没有奇异部分
OPE[B, C] = OPE.make([])  # 没有奇异部分

result = OPE(NO(A, B), C)
print(result)  # 输出: (NO(NO(A,B),C))/(z-w)^0 ✅
print(result.pole(0))  # 输出: NO(NO(A,B),C) ✅
```

## 影响范围

这个修复影响以下场景：
1. 当计算 `OPE(NO(A,B), C)` 且 `OPE(A,C)` 和 `OPE(B,C)` 都没有奇异部分时
2. 使用 `bracket(NO(A,B), C, 0)` 提取 q=0 系数时

修复后的行为符合 VOA 理论和 Thielemans 论文的定义。

## 结论

这是一个重要的 bug 修复，确保了代码与 VOA 理论的一致性。修复后：
- ✅ 正确处理了 q=0 的正规序乘积
- ✅ 符合 Thielemans eq. (3.3.19) 的定义
- ✅ 没有破坏任何现有功能
- ✅ 增加了 5 个新的测试用例

## 参考文献

- Thielemans, K. (1994). "An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings"
  - Section 2.3: Operator Product Expansions
  - eq. (2.3.3): OPE 的完整定义
  - eq. (2.3.7): 正规序乘积的定义
  - eq. (3.3.18): 复合算符 OPE 的奇异部分
  - eq. (3.3.19): 复合算符 OPE 的 q=0 部分
