# PyOPE 测试框架文档

本文档描述 PyOPE 测试框架的设计、使用方法和扩展指南。

## 概述

PyOPE 测试框架基于 pytest，用于验证 Python 实现与 Mathematica 参考实现（OPEdefs.m）的一致性。测试覆盖：

- **SL(2) Kac-Moody 代数**：嵌套正规序、导数算符、Sugawara 构造
- **W3 代数**：基本 OPE、复合算符、数值验证
- **其他代数**：可扩展支持更多 VOA 结构

## 文件结构

```
tests/
├── conftest.py                    # pytest 配置和共享 fixtures
├── algebras/
│   ├── test_sl2_nested_no.py      # SL(2) 嵌套 NO 测试（22个用例）
│   └── test_w3_algebra_ref.py     # W3 代数测试（8个用例）
├── utils/                         # 测试工具模块
│   ├── __init__.py
│   └── comparison.py              # 表达式比较工具
├── w3_algebra_test.wls            # Mathematica 参考实现（W3）
├── test_sl2_mathematica.wls       # Mathematica 参考实现（SL2）
└── TEST_FRAMEWORK.md              # 本文档

pytest.ini                         # pytest 配置文件
```

## 核心组件

### 1. Fixtures (`conftest.py`)

#### `sl2_algebra` Fixture

定义 SL(2) Kac-Moody 代数的三个生成元及其 OPE 规则。

**使用示例**：
```python
def test_example(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    Jm = sl2_algebra['Jminus']
    k = sl2_algebra['k']
    
    # 计算表达式
    result = simplify(NO(Jp, J0))
    
    # 测试完成后清理
    sl2_algebra['clear']()
```

**OPE 规则**（level k）：
- `J0(z) J0(w) = k/2 / (z-w)^2`
- `J0(z) J+(w) = J+(w) / (z-w)`
- `J0(z) J-(w) = -J-(w) / (z-w)`
- `J+(z) J-(w) = k / (z-w)^2 + 2*J0(w) / (z-w)`

#### `w3_algebra` Fixture

定义 W3 代数的两个生成元及其 OPE 规则。

**使用示例**：
```python
def test_example(w3_algebra):
    T = w3_algebra['T']
    W = w3_algebra['W']
    c = w3_algebra['c']
    beta = w3_algebra['beta']
    Lambda = w3_algebra['Lambda']
    
    # 计算 OPE
    result = OPE(T, W)
    
    # 测试完成后清理
    w3_algebra['clear']()
```

**OPE 规则**：
- `T(z) T(w) = c/2 / (z-w)^4 + 2*T(w) / (z-w)^2 + ∂T(w) / (z-w)`
- `T(z) W(w) = 3*W(w) / (z-w)^2 + ∂W(w) / (z-w)`
- `W(z) W(w) = c / (z-w)^6 + ... (见代码)`

**辅助算符**：
- `Λ = NO(T,T) - (3/10)*∂²T`

### 2. 比较工具 (`tests/utils/comparison.py`)

#### `canonicalize(expr, *, expand=True)`

将表达式规范化为标准形式。

**功能**：
- 使用 `simplify()` 展开和化简
- 合并同类项
- 标准化导数表示

**示例**：
```python
from tests.utils.comparison import canonicalize

T = BasicOperator("T")
expr = NO(T, T) + 2*T - T
result = canonicalize(expr)  # 返回 NO(T, T) + T
```

#### `assert_voa_equal(actual, expected, *, msg=None, canonicalize_first=True)`

断言两个 VOA 表达式相等（符号比较）。

**功能**：
- 自动规范化表达式
- 处理项顺序差异
- 提供清晰的错误消息

**示例**：
```python
from tests.utils.comparison import assert_voa_equal

actual = simplify(NO(T, T))
expected = NO(T, T)
assert_voa_equal(actual, expected)
```

#### `assert_voa_numeric_equal(actual, expected, *, subs, tol=1e-12, msg=None)`

断言两个 VOA 表达式在数值代入后相等（数值比较）。

**功能**：
- 符号参数数值代入
- 容差比较（默认 1e-12）
- 支持 OPEData 逐极点比较

**示例**：
```python
from tests.utils.comparison import assert_voa_numeric_equal
import sympy as sp

c = sp.Symbol('c')
actual = c * T
expected = 100 * T
assert_voa_numeric_equal(actual, expected, subs={c: 100})
```

#### `compare_expressions(expr1, expr2, *, canonicalize_first=True)`

比较两个表达式是否相等（返回布尔值）。

**示例**：
```python
from tests.utils.comparison import compare_expressions

expr1 = T + 2*T
expr2 = 3*T
assert compare_expressions(expr1, expr2)  # True
```

## 测试用例组织

### SL(2) 测试 (`algebras/test_sl2_nested_no.py`)

**22个测试用例**，分为5个类别：

1. **简单嵌套 NO**（6个测试：1.1-1.6）
   - `test_1_1__NO_NO_Jplus_Jzero_Jminus`
   - `test_1_2__NO_Jplus_NO_Jzero_Jminus`
   - 等等...

2. **三重嵌套 NO**（4个测试：2.1-2.4）
   - `test_2_1__NO_NO_Jminus_NO_Jplus_Jzero`
   - 等等...

3. **带导数的嵌套 NO**（5个测试：3.1-3.5）
   - `test_3_1__NO_NO_dJplus_Jzero_Jminus`
   - 等等...

4. **Sugawara 构造**（3个测试：4.1-4.3）
   - `test_4_1__NO_NO_Jplus_Jminus_NO_Jplus_Jminus`
   - 等等...

5. **复杂嵌套结构**（4个测试：5.1-5.4）
   - `test_5_1__NO_NO_NO_Jplus_Jminus_Jzero_Jplus`
   - 等等...

### W3 测试 (`algebras/test_w3_algebra_ref.py`)

**8个测试用例**，分为4个类别：

1. **基本 OPE**（3个测试）
   - `test_1__T_T_OPE`
   - `test_2__T_W_OPE`
   - `test_3__W_W_OPE`

2. **导数算符 OPE**（2个测试）
   - `test_4__T_with_NO_TT_minus_derivative`
   - `test_5__T_with_W_double_derivative`

3. **复合正规序 OPE**（2个测试）
   - `test_6__NO_dT_NO_ddW_T`
   - `test_7__NO_dW_NO_ddW_T`

4. **数值验证**（1个测试）
   - `test_8__numeric_verification_c100_beta_0p1`

## Pytest Markers

### 内置 Markers

- `@pytest.mark.mathematica_ref`：标记为参考 Mathematica 实现的测试
- `@pytest.mark.slow`：标记为慢速测试
- `@pytest.mark.requires_derivative`：标记为需要导数支持的测试
- `@pytest.mark.requires_no`：标记为需要正规序支持的测试

### 使用示例

```python
@pytest.mark.mathematica_ref
@pytest.mark.slow
class TestComplexCalculations:
    def test_heavy_computation(self, sl2_algebra):
        # 复杂计算...
        pass
```

### 选择性运行

```bash
# 跳过慢速测试
pytest -m "not slow"

# 只运行 Mathematica 参考测试
pytest -m "mathematica_ref"

# 跳过需要导数支持的测试
pytest -m "not requires_derivative"
```

## 运行测试

### 基本用法

```bash
# 运行所有测试
pytest

# 运行特定文件
pytest tests/algebras/test_sl2_nested_no.py

# 运行特定测试类
pytest tests/algebras/test_sl2_nested_no.py::TestSimpleNestedNO

# 运行特定测试
pytest tests/algebras/test_sl2_nested_no.py::TestSimpleNestedNO::test_1_1__NO_NO_Jplus_Jzero_Jminus

# 详细输出
pytest -v

# 显示所有输出（包括 print）
pytest -s

# 失败时进入调试器
pytest --pdb
```

### 高级用法

```bash
# 并行运行（需要 pytest-xdist）
pytest -n auto

# 生成覆盖率报告（需要 pytest-cov）
pytest --cov=pyope --cov-report=html

# 只运行失败的测试
pytest --lf

# 运行失败的测试，然后运行其他测试
pytest --ff

# 在第一个失败时停止
pytest -x

# 在 N 个失败后停止
pytest --maxfail=3
```

## 添加新测试

### 步骤 1：创建 Mathematica 参考实现

在 `tests/` 目录下创建 `.wls` 文件：

```mathematica
#!/usr/bin/env wolframscript
(* 新代数测试 *)

Get["/path/to/OPEdefs.m"];

(* 定义算符 *)
Bosonic[A, B, C]

(* 定义 OPE *)
OPE[A, B] = MakeOPE[...];

(* 测试用例 *)
expr1 = NO[A, NO[B, C]];
result1 = Expand[expr1];
Print["Test 1: ", result1];
```

### 步骤 2：在 `conftest.py` 中添加 Fixture

```python
@pytest.fixture
def new_algebra():
    """新代数 fixture"""
    # 定义算符
    A = BasicOperator("A")
    B = BasicOperator("B")
    C = BasicOperator("C")
    
    Bosonic(A, B, C)
    
    # 定义 OPE
    OPE[A, B] = MakeOPE([...])
    
    def clear():
        # 清理逻辑
        pass
    
    return {
        'A': A,
        'B': B,
        'C': C,
        'clear': clear,
    }
```

### 步骤 3：创建测试文件

```python
# tests/algebras/test_new_algebra.py

import pytest
from pyope import NO, simplify
from tests.utils.comparison import assert_voa_equal

@pytest.mark.mathematica_ref
class TestNewAlgebra:
    def test_1__basic_no(self, new_algebra):
        """Test 1: 基本正规序"""
        A = new_algebra['A']
        B = new_algebra['B']
        C = new_algebra['C']
        
        expr = NO(A, NO(B, C))
        result = simplify(expr)
        
        # 验证结果
        assert result is not None
```

### 步骤 4：运行并验证

```bash
# 运行新测试
pytest tests/algebras/test_new_algebra.py -v

# 与 Mathematica 结果对比
wolframscript tests/new_algebra_test.wls > mathematica_output.txt
# 手动对比或编写自动对比脚本
```

## 调试技巧

### 1. 打印中间结果

```python
def test_debug(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    
    expr = NO(Jp, J0)
    print(f"Expression: {expr}")
    
    result = simplify(expr)
    print(f"Simplified: {result}")
    
    assert result is not None
```

运行：`pytest -s tests/algebras/test_sl2_nested_no.py::test_debug`

### 2. 使用 pytest 调试器

```python
def test_debug(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    
    expr = NO(Jp, J0)
    
    # 设置断点
    import pdb; pdb.set_trace()
    
    result = simplify(expr)
    assert result is not None
```

运行：`pytest --pdb`

### 3. 比较失败时的详细信息

```python
from tests.utils.comparison import assert_voa_equal

def test_comparison(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    
    actual = 2 * Jp
    expected = 3 * Jp
    
    # 会显示详细的差异
    assert_voa_equal(actual, expected, msg="Custom error message")
```

### 4. 检查 OPE 定义

```python
def test_check_ope(sl2_algebra):
    from pyope.registry import ope_registry
    
    Jp = sl2_algebra['Jplus']
    Jm = sl2_algebra['Jminus']
    
    # 检查 OPE 是否已定义
    key = (Jp, Jm)
    if key in ope_registry._opes:
        print(f"OPE defined: {ope_registry._opes[key]}")
    else:
        print(f"OPE not defined for {key}")
```

## 与 Mathematica 对比策略

### 策略 1：结构比较

对于简单表达式，直接比较结构：

```python
# Python
result = simplify(NO(Jp, J0))

# Mathematica
(* result = Expand[NO[Jplus, Jzero]] *)

# 比较：手动检查项和系数
```

### 策略 2：数值验证

对于复杂表达式，使用数值代入：

```python
# Python
result = OPE(W, W)
assert_voa_numeric_equal(
    result.pole(6),
    100,
    subs={c: 100, beta: sp.Rational(1, 10)}
)

# Mathematica
(* OPE[W, W] /. {c -> 100, β -> 1/10} *)
```

### 策略 3：自动化对比（高级）

创建脚本自动运行 Mathematica 并解析输出：

```python
import subprocess
import re

def run_mathematica_test(wls_file):
    """运行 Mathematica 测试并解析输出"""
    result = subprocess.run(
        ['wolframscript', wls_file],
        capture_output=True,
        text=True
    )
    
    # 解析输出
    output = result.stdout
    # ... 解析逻辑 ...
    
    return parsed_results

def test_with_mathematica_comparison(sl2_algebra):
    # 运行 Mathematica
    mathematica_results = run_mathematica_test('tests/test_sl2_mathematica.wls')
    
    # 运行 Python
    python_result = simplify(NO(Jp, J0))
    
    # 比较
    assert compare_with_mathematica(python_result, mathematica_results['test_1_1'])
```

## 常见问题

### Q1: 测试失败，但结果看起来正确？

**A**: 可能是项顺序或表示形式不同。尝试：
1. 使用 `canonicalize()` 规范化
2. 使用数值比较 `assert_voa_numeric_equal()`
3. 检查是否有符号差异（如 `-1` vs `(-1)`）

### Q2: 如何跳过某些测试？

**A**: 使用 `@pytest.mark.skip` 或条件跳过：

```python
@pytest.mark.skipif(not supports_derivative(), reason="Derivative not supported")
def test_with_derivative(sl2_algebra):
    # ...
```

### Q3: 测试运行很慢？

**A**: 
1. 使用 `-m "not slow"` 跳过慢速测试
2. 使用 `pytest-xdist` 并行运行：`pytest -n auto`
3. 只运行特定测试文件或类

### Q4: 如何验证 Jacobi 恒等式？

**A**: 参考现有的 Jacobi 测试：

```python
def test_jacobi(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    Jm = sl2_algebra['Jminus']
    
    # 计算 Jacobi 恒等式的三项
    term1 = simplify(NO(Jp, NO(J0, Jm)))
    term2 = simplify(NO(J0, NO(Jp, Jm)))
    term3 = simplify(NO(NO(Jp, J0), Jm))
    
    # 验证恒等式
    # {A{BC}_p}_q = (-1)^{|A||B|} {B{AC}_q}_p + Σ ...
```

## 扩展指南

### 添加新的比较函数

在 `tests/utils/comparison.py` 中添加：

```python
def assert_ope_pole_equal(actual_ope, expected_ope, pole, **kwargs):
    """比较 OPE 的特定极点"""
    actual_coeff = actual_ope.pole(pole)
    expected_coeff = expected_ope.pole(pole)
    assert_voa_equal(actual_coeff, expected_coeff, **kwargs)
```

### 添加新的 Fixture

在 `conftest.py` 中添加：

```python
@pytest.fixture
def virasoro_algebra():
    """Virasoro 代数 fixture"""
    # 实现...
    pass
```

### 添加新的 Marker

在 `pytest.ini` 中添加：

```ini
markers =
    integration: 标记为集成测试
```

在 `conftest.py` 中配置：

```python
def pytest_configure(config):
    config.addinivalue_line("markers", "integration: integration tests")
```

## 参考资料

- [pytest 文档](https://docs.pytest.org/)
- [Mathematica OPEdefs.m](../OPEdefs/OPEdefs.m)
- [PyOPE README](../README.md)
- [Thielemans 论文](../papers/)

## 贡献指南

1. 添加新测试时，确保：
   - 有对应的 Mathematica 参考实现
   - 测试名称清晰描述测试内容
   - 包含详细的文档字符串
   - 使用适当的 markers

2. 修改比较工具时，确保：
   - 向后兼容
   - 添加单元测试
   - 更新文档

3. 提交前运行：
   ```bash
   pytest tests/
   ```

## 更新日志

- **2024-01-XX**: 初始版本
  - 创建 SL(2) 和 W3 测试
  - 实现比较工具
  - 编写测试框架文档
