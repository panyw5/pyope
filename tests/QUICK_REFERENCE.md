# PyOPE 测试快速参考

## 快速开始

### 运行新创建的测试

```bash
# 运行所有新测试（32个）
pytest tests/test_sl2_nested_no.py tests/test_w3_algebra_ref.py -v

# 只运行 SL(2) 测试（22个）
pytest tests/test_sl2_nested_no.py -v

# 只运行 W3 测试（10个）
pytest tests/test_w3_algebra_ref.py -v
```

### 测试结果
✅ **32/32 测试通过 (100%)**
- SL(2) Kac-Moody: 22/22 ✅
- W3 代数: 10/10 ✅

## 文件清单

| 文件 | 说明 | 行数 |
|------|------|------|
| `tests/utils/comparison.py` | 表达式比较工具 | 430 |
| `tests/conftest.py` | Pytest fixtures | 260 |
| `tests/test_sl2_nested_no.py` | SL(2) 测试 | 422 |
| `tests/test_w3_algebra_ref.py` | W3 测试 | 330 |
| `tests/TEST_FRAMEWORK.md` | 详细文档 | 600+ |
| `tests/TESTING_SUMMARY.md` | 实施总结 | 400+ |
| `pytest.ini` | Pytest 配置 | 40 |

## 核心 API

### 比较工具

```python
from tests.utils.comparison import (
    assert_voa_equal,           # 符号比较
    assert_voa_numeric_equal,   # 数值比较
    canonicalize,               # 规范化
)

# 符号比较
assert_voa_equal(actual, expected)

# 数值比较
assert_voa_numeric_equal(actual, expected, subs={c: 100, beta: 0.1})
```

### Fixtures

```python
# SL(2) Kac-Moody 代数
def test_example(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    Jm = sl2_algebra['Jminus']
    k = sl2_algebra['k']

# W3 代数
def test_example(w3_algebra):
    T = w3_algebra['T']
    W = w3_algebra['W']
    c = w3_algebra['c']
    beta = w3_algebra['beta']
    Lambda = w3_algebra['Lambda']
```

## 测试示例

### SL(2) 嵌套 NO 测试

```python
@pytest.mark.mathematica_ref
def test_nested_no(sl2_algebra):
    Jp = sl2_algebra['Jplus']
    J0 = sl2_algebra['Jzero']
    Jm = sl2_algebra['Jminus']
    
    # 计算嵌套正规序
    expr = NO(NO(Jp, J0), Jm)
    result = simplify(expr)
    
    # 验证结果
    assert result is not None
```

### W3 OPE 测试

```python
@pytest.mark.mathematica_ref
def test_w3_ope(w3_algebra):
    T = w3_algebra['T']
    c = w3_algebra['c']
    
    # 计算 OPE
    result = OPE(T, T)
    
    # 验证极点
    assert result.max_pole == 4
    assert_voa_equal(result.pole(4), sp.Rational(1, 2) * c)
    assert_voa_equal(result.pole(2), 2 * T)
```

## Pytest Markers

```bash
# 跳过慢速测试
pytest -m "not slow"

# 只运行 Mathematica 参考测试
pytest -m "mathematica_ref"

# 跳过需要导数的测试
pytest -m "not requires_derivative"
```

## 常用命令

```bash
# 详细输出
pytest -v

# 显示 print 输出
pytest -s

# 失败时进入调试器
pytest --pdb

# 只运行失败的测试
pytest --lf

# 在第一个失败时停止
pytest -x

# 并行运行（需要 pytest-xdist）
pytest -n auto
```

## 对应关系

### Mathematica → PyOPE

| Mathematica | PyOPE | 说明 |
|-------------|-------|------|
| `Bosonic[A, B]` | `Bosonic(A, B)` | 声明玻色算符 |
| `NO[A, B]` | `NO(A, B)` | 正规序 |
| `OPE[A, B]` | `OPE(A, B)` | 计算 OPE |
| `Derivative[1][A]` | `d(A)` | 一阶导数 |
| `A''` | `dn(2, A)` | 二阶导数 |
| `Expand[expr]` | `simplify(expr)` | 展开/化简 |

### 测试编号对应

**SL(2) 测试**：
- 1.1-1.6: 简单嵌套 NO
- 2.1-2.4: 三重嵌套 NO
- 3.1-3.5: 带导数的嵌套 NO
- 4.1-4.3: Sugawara 构造
- 5.1-5.4: 复杂嵌套结构

**W3 测试**：
- Test 1-3: 基本 OPE (T-T, T-W, W-W)
- Test 4-5: 导数算符 OPE
- Test 6-7: 复合正规序
- Test 8: 数值验证

## 参考文档

- **详细文档**: `tests/TEST_FRAMEWORK.md`
- **实施总结**: `tests/TESTING_SUMMARY.md`
- **Mathematica 参考**:
  - `tests/test_sl2_mathematica.wls`
  - `tests/w3_algebra_test.wls`

## 问题排查

### 测试失败？

1. 检查表达式规范化：`canonicalize(expr)`
2. 使用数值比较：`assert_voa_numeric_equal()`
3. 查看详细输出：`pytest -v -s`
4. 进入调试器：`pytest --pdb`

### 添加新测试？

1. 创建 Mathematica 参考 (`.wls`)
2. 在 `conftest.py` 添加 fixture
3. 创建测试文件
4. 使用 `@pytest.mark.mathematica_ref`
5. 运行并验证

详见 `tests/TEST_FRAMEWORK.md` 的"添加新测试"章节。

## 贡献

提交前请确保：
```bash
# 运行所有测试
pytest tests/

# 运行新测试
pytest tests/test_sl2_nested_no.py tests/test_w3_algebra_ref.py -v
```

所有测试应该通过！✅
