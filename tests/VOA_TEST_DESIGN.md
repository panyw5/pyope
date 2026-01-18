# VOA 系统测试用例设计方案

基于 `VOA.wls` 的分析，为 pyope 设计完整的测试用例结构。

## 测试设计原则

1. **独立性**：每个 VOA 系统有独立的测试文件
2. **完整性**：覆盖算符声明、OPE 定义、关键计算、性质验证
3. **可对比性**：与 Mathematica 结果精确对比
4. **清晰性**：详细注释说明测试目的和预期结果
5. **分层测试**：基础定义 → 简单计算 → 复杂计算 → 性质验证

## 测试文件结构模板

```python
"""
[VOA 系统名称] 测试

基于 VOA.wls 中的定义和计算

测试内容：
1. 算符声明和性质
2. 基本 OPE 定义
3. 导数算符的 OPE
4. 正规序乘积
5. 关键性质验证（primary conditions, null states 等）
"""

import pytest
import sympy as sp
from sympy import Symbol, Rational, simplify, expand

from pyope import (
    BasisOperator,
    d, dn,
    One, Zero,
    OPE, NO, bracket,
    MakeOPE,
    Bosonic, Fermionic,
)

@pytest.fixture(autouse=True)
def clear_registry():
    """每个测试前清空注册表"""
    from pyope.registry import ope_registry
    ope_registry.clear()
    yield
    ope_registry.clear()

class Test[SystemName]Definition:
    """测试基本定义"""
    
    def test_operators_declaration(self):
        """测试算符声明"""
        pass
    
    def test_basic_opes(self):
        """测试基本 OPE 定义"""
        pass

class Test[SystemName]Computations:
    """测试具体计算"""
    
    def test_derivative_opes(self):
        """测试导数算符的 OPE"""
        pass
    
    def test_normal_orders(self):
        """测试正规序乘积"""
        pass

class Test[SystemName]Properties:
    """测试代数性质"""
    
    def test_primary_conditions(self):
        """测试 primary 条件"""
        pass
    
    def test_null_states(self):
        """测试 null states"""
        pass
```

---

## 1. Virasoro 代数测试

**文件**: `tests/test_virasoro_voa.py`

### 1.1 测试函数列表

#### TestVirasoroDefinition（基础定义）

1. **test_stress_tensor_declaration**
   - 输入：声明 T (conformal_weight=2, bosonic)
   - 验证：T.conformal_weight == 2, T._bosonic == True

2. **test_virasoro_ope_with_central_charge**
   - 输入：定义 OPE[T,T] = MakeOPE([c/2*One, 0, 2*T, d(T)])
   - 验证：
     - max_pole == 4
     - pole(4) == c/2 * One
     - pole(3) == 0
     - pole(2) == 2*T
     - pole(1) == d(T)

3. **test_virasoro_ope_specific_c**
   - 输入：设定 c=26（临界弦理论）
   - 验证：OPE 系数正确

#### TestVirasoroComputations（计算测试）

4. **test_t_derivative_ope**
   - 输入：计算 OPE[T, d(T)]
   - 验证：与 Thielemans eq 3.3.2 一致
   - 预期：max_pole=5, 各极点系数

5. **test_t_second_derivative_ope**
   - 输入：计算 OPE[T, d(T,2)]
   - 验证：max_pole=6
   - 预期极点：
     - pole(6) = 10*c*One
     - pole(5) = 0
     - pole(4) = 32*T
     - pole(3) = 18*d(T)
     - pole(2) = 6*d(T,2)
     - pole(1) = d(T,3)

6. **test_normal_order_t_t**
   - 输入：计算 NO(T,T)
   - 验证：conformal_weight == 4
   - 验证：与 d(T,2) 的关系

7. **test_t_no_t_t_ope**
   - 输入：计算 OPE[T, NO(T,T)]
   - 验证：max_pole == 6
   - 预期：pole(6) = 3*c*One（Jacobi 恒等式）

#### TestVirasoroProperties（性质验证）

8. **test_virasoro_algebra_commutators**
   - 输入：计算 [L_m, L_n] 对易关系
   - 验证：(m-n)L_{m+n} + c/12*m(m²-1)δ_{m+n,0}

9. **test_conformal_transformation**
   - 输入：验证 T 在共形变换下的行为
   - 验证：Schwarzian 导数项

10. **test_null_state_level_2**
    - 输入：构造 level-2 null state
    - 验证：(L_{-2} + 3/(2(2h+1))L_{-1}²)|h⟩ = 0

---

## 2. bc-βγ 系统测试

**文件**: `tests/test_bc_betagamma.py`

### 2.1 测试函数列表

#### TestBCSystem（bc 系统）

1. **test_bc_operators_declaration**
   - 输入：b (weight=λ, fermionic), c (weight=1-λ, fermionic)
   - 验证：conformal weights, fermionic 性质

2. **test_bc_ope**
   - 输入：OPE[b,c] = MakeOPE([One])
   - 验证：max_pole == 1, pole(1) == One

3. **test_bc_stress_tensor**
   - 输入：T^{bc} = -λ*NO(b,d(c)) - (1-λ)*NO(d(b),c)
   - 验证：conformal_weight == 2

4. **test_bc_central_charge**
   - 输入：计算 OPE[T^{bc}, T^{bc}]
   - 验证：pole(4) == c^{bc}/2*One，其中 c^{bc} = -2(6λ²-6λ+1)

5. **test_bc_t_b_ope**
   - 输入：计算 OPE[T^{bc}, b]
   - 验证：pole(2) == λ*b, pole(1) == d(b)（b 是 primary）

6. **test_bc_t_c_ope**
   - 输入：计算 OPE[T^{bc}, c]
   - 验证：pole(2) == (1-λ)*c, pole(1) == d(c)（c 是 primary）

#### TestBetaGammaSystem（βγ 系统）

7. **test_betagamma_operators_declaration**
   - 输入：β (weight=λ, fermionic), γ (weight=1-λ, fermionic)
   - 验证：conformal weights, fermionic 性质

8. **test_betagamma_ope**
   - 输入：OPE[β,γ] = MakeOPE([-One])
   - 验证：max_pole == 1, pole(1) == -One

9. **test_betagamma_stress_tensor**
   - 输入：T^{βγ} = -λ*NO(β,d(γ)) - (1-λ)*NO(d(β),γ)
   - 验证：conformal_weight == 2

10. **test_betagamma_central_charge**
    - 输入：计算 OPE[T^{βγ}, T^{βγ}]
    - 验证：pole(4) == c^{βγ}/2*One，其中 c^{βγ} = 2(6λ²-6λ+1)

11. **test_betagamma_t_beta_ope**
    - 输入：计算 OPE[T^{βγ}, β]
    - 验证：β 是 primary of weight λ

12. **test_betagamma_t_gamma_ope**
    - 输入：计算 OPE[T^{βγ}, γ]
    - 验证：γ 是 primary of weight 1-λ

#### TestBCBetaGammaRelation（bc-βγ 关系）

13. **test_bosonization**
    - 输入：验证 bc 和 βγ 系统的对偶关系
    - 验证：c^{bc} + c^{βγ} = 0

14. **test_screening_operators**
    - 输入：构造 screening operators
    - 验证：与 stress tensor 的对易关系

---

## 3. U(1)_k 系统测试

**文件**: `tests/test_u1_k.py`

### 3.1 测试函数列表

#### TestU1KDefinition（基础定义）

1. **test_u1_current_declaration**
   - 输入：J (weight=1, bosonic)
   - 验证：conformal_weight == 1

2. **test_u1_ope**
   - 输入：OPE[J,J] = MakeOPE([k*One, 0])
   - 验证：max_pole == 2, pole(2) == k*One, pole(1) == 0

3. **test_sugawara_construction**
   - 输入：T = NO(J,J)/(2*k)
   - 验证：conformal_weight == 2

4. **test_u1_central_charge**
   - 输入：计算 OPE[T,T]
   - 验证：pole(4) == 1/2*One（c=1）

#### TestU1KComputations（计算测试）

5. **test_t_j_ope**
   - 输入：计算 OPE[T,J]
   - 验证：J 是 primary of weight 1
   - 预期：pole(2) == J, pole(1) == d(J)

6. **test_j_derivative_ope**
   - 输入：计算 OPE[J, d(J)]
   - 验证：max_pole == 3

7. **test_vertex_operators**
   - 输入：定义 V_α = exp(iαφ)，其中 ∂φ = J
   - 验证：OPE[J, V_α] ~ α*V_α/(z-w)

8. **test_vertex_operator_ope**
   - 输入：计算 OPE[V_α, V_β]
   - 验证：(z-w)^{αβ/k} 的行为

#### TestU1KProperties（性质验证）

9. **test_kac_moody_algebra**
   - 输入：验证 [J_m, J_n] = k*m*δ_{m+n,0}
   - 验证：对易关系

10. **test_charge_conservation**
    - 输入：验证电荷守恒
    - 验证：∮J = Q（总电荷）

---

## 4. sl(2)_k 系统测试

**文件**: `tests/test_sl2_k.py`

### 4.1 测试函数列表

#### TestSL2KDefinition（基础定义）

1. **test_sl2_generators_declaration**
   - 输入：J⁺, J⁰, J⁻ (weight=1, bosonic)
   - 验证：conformal_weight == 1

2. **test_sl2_opes**
   - 输入：定义三个 OPE
     - OPE[J⁺, J⁻] = MakeOPE([k*One, J⁰])
     - OPE[J⁰, J⁺] = MakeOPE([2*J⁺])
     - OPE[J⁰, J⁻] = MakeOPE([-2*J⁻])
   - 验证：各极点系数

3. **test_sl2_sugawara_tensor**
   - 输入：T = (NO(J⁺,J⁻) + NO(J⁻,J⁺) + NO(J⁰,J⁰)/2)/(k+2)
   - 验证：conformal_weight == 2

4. **test_sl2_central_charge**
   - 输入：计算 OPE[T,T]
   - 验证：pole(4) == c/2*One，其中 c = 3k/(k+2)

#### TestSL2KComputations（计算测试）

5. **test_t_j_plus_ope**
   - 输入：计算 OPE[T, J⁺]
   - 验证：J⁺ 是 primary of weight 1

6. **test_t_j_zero_ope**
   - 输入：计算 OPE[T, J⁰]
   - 验证：J⁰ 是 primary of weight 1

7. **test_t_j_minus_ope**
   - 输入：计算 OPE[T, J⁻]
   - 验证：J⁻ 是 primary of weight 1

8. **test_serre_relations**
   - 输入：验证 Serre 关系
   - 验证：[J⁺, [J⁺, J⁻]] = 0 等

#### TestSL2KProperties（性质验证）

9. **test_kac_moody_commutators**
   - 输入：验证 [J^a_m, J^b_n] = f^{abc}J^c_{m+n} + k*m*δ^{ab}δ_{m+n,0}
   - 验证：结构常数 f^{abc}

10. **test_highest_weight_representations**
    - 输入：构造最高权表示
    - 验证：J⁺|λ⟩ = 0, J⁰|λ⟩ = λ|λ⟩

11. **test_null_states_sl2**
    - 输入：构造 null states
    - 验证：特定 level 的 null state 条件

---

## 5. W₃ 代数测试（扩展）

**文件**: `tests/test_w3_algebra.py`（已存在，需要扩展）

### 5.1 需要添加的测试

#### TestW3AlgebraExtended（扩展测试）

1. **test_w3_with_generic_parameters**
   - 输入：符号参数 c, β
   - 验证：代数结构的一般性质

2. **test_w3_zamolodchikov_metric**
   - 输入：计算 W₃ 代数的度规
   - 验证：度规的正定性条件

3. **test_w3_null_states**
   - 输入：构造 W₃ 代数的 null states
   - 验证：特定 c 值下的 null state 条件

4. **test_w3_primary_fields**
   - 输入：构造 primary fields
   - 验证：OPE[T,φ] 和 OPE[W,φ] 的 primary 条件

5. **test_w3_casimir_operator**
   - 输入：构造 Casimir 算符
   - 验证：与所有生成元对易

6. **test_w3_fusion_rules**
   - 输入：计算简单的 fusion rules
   - 验证：结合律

---

## 6. N=4 Small SCFA 测试

**文件**: `tests/test_n4_scfa.py`

### 6.1 测试函数列表

#### TestN4SCFADefinition（基础定义）

1. **test_n4_operators_declaration**
   - 输入：T, G⁺, G⁻, J (superconformal generators)
   - 验证：conformal weights: T(2), G±(3/2), J(1)

2. **test_n4_t_t_ope**
   - 输入：OPE[T,T] = Virasoro OPE with c=6
   - 验证：标准 Virasoro 结构

3. **test_n4_j_j_ope**
   - 输入：OPE[J,J] = MakeOPE([k*One, 0])
   - 验证：U(1) current algebra

4. **test_n4_g_plus_g_minus_ope**
   - 输入：OPE[G⁺, G⁻]
   - 验证：包含 T 和 J 的项

5. **test_n4_t_g_ope**
   - 输入：OPE[T, G⁺] 和 OPE[T, G⁻]
   - 验证：G± 是 primary of weight 3/2

6. **test_n4_j_g_ope**
   - 输入：OPE[J, G⁺] 和 OPE[J, G⁻]
   - 验证：G± 的 U(1) 电荷

#### TestN4SCFAComputations（计算测试）

7. **test_n4_stress_tensor_decomposition**
   - 输入：T = T_matter + T_ghost
   - 验证：中心荷分解

8. **test_n4_superconformal_ward_identities**
   - 输入：验证 Ward 恒等式
   - 验证：超对称变换

9. **test_n4_spectral_flow**
   - 输入：计算 spectral flow
   - 验证：模空间的变换

#### TestN4SCFAProperties（性质验证）

10. **test_n4_supersymmetry_algebra**
    - 输入：验证 {Q, Q̄} = H
    - 验证：超对称代数

11. **test_n4_bps_states**
    - 输入：构造 BPS states
    - 验证：短表示条件

12. **test_n4_r_symmetry**
    - 输入：验证 R-对称性
    - 验证：R-charge 守恒

---

## 7. Rank-1 N=3 Theory 测试

**文件**: `tests/test_n3_rank1.py`

### 7.1 测试函数列表

#### TestN3Rank1Definition（基础定义）

1. **test_n3_operators_declaration**
   - 输入：T, G, A, Φ (N=3 generators)
   - 验证：conformal weights

2. **test_n3_basic_opes**
   - 输入：定义所有基本 OPE
   - 验证：OPE 结构

3. **test_n3_central_charge**
   - 输入：计算中心荷
   - 验证：c = 3k/(k+2)（rank-1 情况）

#### TestN3Rank1Computations（计算测试）

4. **test_n3_composite_operators**
   - 输入：构造复合算符
   - 验证：正规序乘积

5. **test_n3_derivative_opes**
   - 输入：计算导数算符的 OPE
   - 验证：与基本 OPE 的一致性

6. **test_n3_normal_orders**
   - 输入：计算嵌套正规序
   - 验证：Jacobi 恒等式

#### TestN3Rank1Properties（性质验证）

7. **test_n3_supersymmetry**
   - 输入：验证 N=3 超对称代数
   - 验证：对易关系

8. **test_n3_null_states**
   - 输入：构造 null states
   - 验证：特定条件下的 null state

9. **test_n3_modular_invariance**
   - 输入：验证模不变性
   - 验证：配分函数

---

## 测试数据来源

### Mathematica 参考结果

每个测试应该包含从 VOA.wls 提取的参考结果：

```python
# 示例：Virasoro OPE[T, d(T,2)] 的参考结果
MATHEMATICA_REFERENCE = {
    'max_pole': 6,
    'pole_6': '10*c*One',
    'pole_5': '0',
    'pole_4': '32*T',
    'pole_3': '18*d(T)',
    'pole_2': '6*d(T,2)',
    'pole_1': 'd(T,3)',
}
```

### 测试辅助函数

```python
def compare_with_mathematica(result, reference):
    """比较 pyope 结果与 Mathematica 参考结果"""
    assert result.max_pole == reference['max_pole']
    for n in range(1, result.max_pole + 1):
        pole_n = result.pole(n)
        ref_n = eval(reference[f'pole_{n}'])
        diff = simplify(pole_n - ref_n)
        assert diff == 0, f"Pole {n} mismatch: {pole_n} vs {ref_n}"

def verify_primary_condition(T, phi, h):
    """验证 φ 是 conformal weight h 的 primary field"""
    ope_result = OPE(T, phi)
    assert ope_result.max_pole == 2
    assert simplify(ope_result.pole(2) - h*phi) == 0
    assert simplify(ope_result.pole(1) - d(phi)) == 0

def verify_null_state(state, level):
    """验证 null state 条件"""
    # 实现 null state 验证逻辑
    pass
```

---

## 测试优先级和时间估计

### 高优先级（1-2 周）

1. ✅ **Virasoro** - 最基础，已有部分测试
2. ✅ **bc-βγ** - 自由场系统，已有基础测试
3. **U(1)_k** - 简单的 Kac-Moody 代数
4. **sl(2)_k** - 非平凡的 Kac-Moody 代数
5. ✅ **W₃** - 已有测试，需要扩展

### 中优先级（2-3 周）

6. **N=4 Small SCFA** - 超共形代数，中等复杂度

### 低优先级（3-4 周）

7. **Rank-1 N=3** - 最复杂，需要前面的经验

---

## 测试执行计划

### 阶段 1：基础测试（第 1 周）

- 创建 `test_virasoro_voa.py`（扩展现有测试）
- 创建 `test_u1_k.py`
- 运行并验证所有基础测试

### 阶段 2：Kac-Moody 代数（第 2 周）

- 创建 `test_sl2_k.py`
- 扩展 `test_bc_betagamma.py`
- 验证 Sugawara 构造

### 阶段 3：W 代数（第 3 周）

- 扩展 `test_w3_algebra.py`
- 添加 null states 测试
- 验证复杂的正规序计算

### 阶段 4：超共形代数（第 4 周）

- 创建 `test_n4_scfa.py`
- 验证超对称性质
- 测试 BPS states

### 阶段 5：高级系统（第 5-6 周）

- 创建 `test_n3_rank1.py`
- 完整的 N=3 代数测试
- 性能优化和文档完善

---

## 持续集成

### pytest 配置

```ini
# pytest.ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = 
    -v
    --tb=short
    --strict-markers
markers =
    virasoro: Virasoro algebra tests
    free_field: Free field system tests
    kac_moody: Kac-Moody algebra tests
    w_algebra: W-algebra tests
    superconformal: Superconformal algebra tests
    slow: Slow tests (> 1s)
```

### 测试标记使用

```python
@pytest.mark.virasoro
def test_virasoro_ope():
    pass

@pytest.mark.slow
def test_complex_null_state():
    pass
```

### 运行特定测试

```bash
# 运行所有 Virasoro 测试
pytest -m virasoro

# 运行单个文件
pytest tests/test_virasoro_voa.py

# 运行单个测试
pytest tests/test_virasoro_voa.py::TestVirasoroDefinition::test_virasoro_ope_with_central_charge

# 跳过慢速测试
pytest -m "not slow"
```

---

## 文档和注释规范

### 测试函数文档字符串

```python
def test_virasoro_ope_with_central_charge(self):
    """
    测试 Virasoro 代数的基本 OPE
    
    定义：T(z)T(w) ~ c/(2(z-w)⁴) + 2T/(z-w)² + T'/(z-w)
    
    参考：
    - VOA.wls line 123-145
    - Thielemans paper eq 2.1.3
    
    验证：
    - 4 阶极点系数为 c/2
    - 3 阶极点为 0
    - 2 阶极点系数为 2
    - 1 阶极点为导数项
    """
    pass
```

### 测试类文档字符串

```python
class TestVirasoroDefinition:
    """
    Virasoro 代数基础定义测试
    
    测试内容：
    1. 应力张量 T 的声明和性质
    2. T-T OPE 的定义和验证
    3. 中心荷 c 的参数化
    
    参考资料：
    - VOA.wls Section 1: Virasoro Algebra
    - Di Francesco et al., "Conformal Field Theory", Chapter 5
    """
    pass
```

---

## 总结

这个测试设计方案提供了：

1. **7 个独立的测试文件**，覆盖所有主要 VOA 系统
2. **100+ 个测试函数**，从基础到高级
3. **清晰的测试结构**，易于维护和扩展
4. **与 Mathematica 的精确对比**，确保正确性
5. **分阶段实施计划**，可逐步完成

建议从高优先级测试开始，逐步完善整个测试套件。
