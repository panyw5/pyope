# U(1)_k 测试修复总结

## 修复日期
2024-01-XX

## 修复的测试文件
`tests/test_u1_k.py`

---

## 修复 1: `test_j_derivative_ope` - 极点系数矛盾

### 问题描述
测试中存在逻辑矛盾：
```python
assert result.max_pole == 3  # 声称最高极点是 3
assert result.pole(3) == 0   # 声称 3 阶极点系数为 0 ❌
```

这在数学上是不可能的：如果 `max_pole == 3`，那么 `pole(3)` 必须非零。

### 数学推导
从基本 OPE `J(z)J(w) ~ k/(z-w)²` 出发：

```
J(z)(∂J)(w) ~ ∂_w[k/(z-w)²]
            = k · ∂_w[(z-w)^(-2)]
            = k · (-2)(z-w)^(-3) · (-1)
            = 2k/(z-w)³
```

因此：
- `max_pole = 3`
- `pole(3) = 2k`

### 修复方案
```python
def test_j_derivative_ope(self):
    """
    测试 OPE[J, ∂J]

    从 J(z)J(w) ~ k/(z-w)² 推导：
    J(z)(∂J)(w) ~ ∂_w[k/(z-w)²] = 2k/(z-w)³

    预期：max_pole = 3, pole(3) = 2k
    """
    J = BasisOperator("J", conformal_weight=1)
    Bosonic(J)

    k = Symbol("k", positive=True)
    OPE[J, J] = MakeOPE([k * One, 0])

    result = OPE(J, d(J))

    assert result.max_pole == 3
    assert result.pole(3) == 2 * k * One  # ✅ 修正

    print(f"✓ OPE[J, ∂J]: max_pole = {result.max_pole}")
    print(f"  - 3-pole: {result.pole(3)}")
```

### 一致性验证
修复后的结果与 `test_j_second_derivative_ope` 一致：
- `∂⁰J`: `OPE[J, J]` → `max_pole=2`, `pole(2)=k`
- `∂¹J`: `OPE[J, ∂J]` → `max_pole=3`, `pole(3)=2k` ✅
- `∂²J`: `OPE[J, ∂²J]` → `max_pole=4`, `pole(4)=6k` ✅

系数遵循导数规则：`∂_w[(z-w)^(-n)] = n(z-w)^(-(n+1))`

---

## 修复 2: `test_primary_field_under_t` - 错误的 OPE 注册

### 问题描述
测试试图为复合算符 `T` 注册 OPE：
```python
T = NO(J, J) / (2 * k)  # T 是复合表达式，不是 BasisOperator
phi = BasisOperator("φ", conformal_weight=2)

# 错误：不能为复合算符注册 OPE
OPE[T, phi] = MakeOPE([h * phi, d(phi)])  # ❌
```

### 为什么错误

1. **设计限制**：`MakeOPE` 只能用于注册 `BasisOperator` 之间的 OPE
2. **T 不是基本算符**：`T = NO(J,J)/(2k)` 是从 `J` 派生的复合表达式
3. **无法推导**：系统无法从 `OPE[J,J]` 推导出任意 `φ` 的行为
4. **测试失败**：注册被忽略，返回 `max_pole=0`（只有正规项）

### 数学背景

在 U(1)_k 理论中：
- **Sugawara 应力张量**：`T = NO(J,J)/(2k)`
- **Primary field 定义**：如果 `φ` 是 weight `h` 的 primary，则
  ```
  T(z)φ(w) ~ h·φ/(z-w)² + ∂φ/(z-w)
  ```

但是：
- 对于**任意** `φ`，这个关系不能从 `OPE[J,J]` 自动推导
- 只有特定的场（如 `J` 本身）才能从基本 OPE 推导出 primary 行为

### Oracle 建议

**正确的测试方法**：
- 测试 `J` 在 Sugawara 张量下的 primary 行为
- 这可以完全从 `OPE[J,J]` 推导出来
- 不需要注册额外的 OPE

### 修复方案

```python
def test_primary_field_under_t(self):
    """
    测试 Sugawara 张量使 J 成为 primary field

    验证 J 在 Sugawara 应力张量 T = NO(J,J)/(2k) 作用下是 weight 1 的 primary field

    预期：T(z)J(w) ~ J/(z-w)² + ∂J/(z-w)

    注意：这个测试从基本的 OPE[J,J] 推导出 OPE[T,J]，
    不需要（也不能）直接注册 OPE[T,J]，因为 T 是复合算符而非 BasisOperator。
    """
    J = BasisOperator("J", conformal_weight=1)
    Bosonic(J)

    k = Symbol("k", positive=True)

    OPE[J, J] = MakeOPE([k * One, 0])
    T = NO(J, J) / (2 * k)

    # 计算 OPE[T, J]（从 OPE[J,J] 推导）
    result = OPE(T, J)

    # J 是 weight 1 的 primary field
    assert result.max_pole == 2
    assert result.pole(2) == J
    assert result.pole(1) == d(J)

    print("✓ Sugawara 张量使 J 成为 weight 1 的 primary field")
    print(f"  - 2-pole: {result.pole(2)}")
    print(f"  - 1-pole: {result.pole(1)}")
```

### 关键改进

1. **移除了 `phi`**：不再测试任意场
2. **测试 `J` 本身**：这是可以从 `OPE[J,J]` 推导的
3. **移除了错误的注册**：不再尝试 `OPE[T, phi] = ...`
4. **依赖推导**：完全依赖 OPE 引擎从基本 OPE 计算复合算符的 OPE

### 为什么这样是正确的

- **数学上**：在任何 U(1)_k 理论中，current `J` 在 Sugawara 张量下都是 weight 1 的 primary
- **实现上**：这个结果可以从 `OPE[J,J]` 和正规序规则机械地推导出来
- **测试上**：验证了 OPE 引擎正确实现了 Sugawara 构造和 Jacobi 恒等式

---

## 测试结果

修复后，所有 16 个测试全部通过：

```
tests/test_u1_k.py::TestU1KDefinition::test_u1_current_declaration PASSED
tests/test_u1_k.py::TestU1KDefinition::test_u1_ope PASSED
tests/test_u1_k.py::TestU1KDefinition::test_sugawara_construction PASSED
tests/test_u1_k.py::TestU1KDefinition::test_u1_central_charge PASSED
tests/test_u1_k.py::TestU1KComputations::test_t_j_ope PASSED
tests/test_u1_k.py::TestU1KComputations::test_j_derivative_ope PASSED ✅
tests/test_u1_k.py::TestU1KComputations::test_j_second_derivative_ope PASSED
tests/test_u1_k.py::TestU1KComputations::test_normal_order_j_j PASSED
tests/test_u1_k.py::TestU1KComputations::test_t_no_j_j_ope PASSED
tests/test_u1_k.py::TestU1KProperties::test_kac_moody_algebra PASSED
tests/test_u1_k.py::TestU1KProperties::test_primary_field_under_t PASSED ✅
tests/test_u1_k.py::TestU1KProperties::test_virasoro_from_sugawara PASSED
tests/test_u1_k.py::TestU1KNumerical::test_u1_k_equals_1 PASSED
tests/test_u1_k.py::TestU1KNumerical::test_u1_k_equals_2 PASSED
tests/test_u1_k.py::TestU1KVertexOperators::test_vertex_operator_definition PASSED
tests/test_u1_k.py::TestU1KVertexOperators::test_j_vertex_operator_ope_concept PASSED

======================== 16 passed in 0.04s ========================
```

---

## 经验教训

### 1. 极点定义的一致性
- `max_pole` 的定义是**最高非零极点的阶数**
- 如果 `max_pole == n`，则 `pole(n) ≠ 0` 必须成立
- 测试断言必须与定义一致

### 2. OPE 注册的限制
- `MakeOPE` 只能用于 `BasisOperator` 之间
- 复合算符的 OPE 必须通过计算推导
- 不能为 `NO(A,B)` 或 `A/c` 这样的表达式注册 OPE

### 3. Primary field 的测试策略
- 测试可以从基本 OPE 推导的性质
- 对于 Sugawara 张量，测试 current 本身是最自然的
- 任意场的 primary 行为需要额外的结构（如 vertex operator 实现）

### 4. 与 Oracle 协作
- Oracle 提供了深入的数学背景和 CFT 理论
- 帮助理解测试的物理意义
- 给出了符合理论和实现的正确方案

---

## 参考资料

1. **Di Francesco et al., "Conformal Field Theory"**
   - Chapter 15: Kac-Moody algebras
   - Sugawara construction

2. **Thielemans, "An Algorithmic Approach to Operator Product Expansions"**
   - OPE 计算规则
   - 导数公式

3. **Oracle 分析**
   - Session: ses_42c77e9a1ffes1LV0mNuZL7kll
   - 提供了详细的数学推导和实现建议
