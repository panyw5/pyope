# pyope 测试套件索引

本文档提供 pyope 库所有测试文件的完整索引和说明。

---

## 测试文件结构

```
tests/
├── 文档
│   ├── JACOBI_IDENTITY_TEST.md          # Jacobi 恒等式测试文档
│   ├── JACOBI_TEST_REPORT.md            # Jacobi 测试详细报告
│   ├── JACOBI_TEST_SUMMARY.md           # Jacobi 测试快速总结
│   ├── W3_ALGEBRA_TESTS.md              # W3 代数测试文档
│   └── TEST_INDEX.md                    # 本文档
│
├── Mathematica 参考测试
│   ├── ref_jacobi_virasoro.wls          # Jacobi 恒等式参考
│   └── w3_algebra_test.wls              # W3 代数参考
│
└── Python 测试文件
    ├── test_api.py                      # API 接口测试
    ├── test_constants.py                # 常量定义测试
    ├── test_ope_data.py                 # OPE 数据结构测试
    ├── test_operator.py                 # 算符基础测试
    ├── test_local_operator.py           # 局域算符测试
    ├── test_thielemans_eqs.py           # Thielemans 方程测试
    ├── test_thielemans_eqs_extended.py  # Thielemans 方程扩展测试
    ├── test_conformal_weight.py         # 共形权重测试
    ├── test_derivative_ope_fix.py       # 导数 OPE 修复测试
    ├── test_composite_left_ope.py       # 复合左 OPE 测试
    ├── test_simplify.py                 # 简化函数测试
    ├── test_advanced_ope.py             # 高级 OPE 测试
    ├── test_ope_examples_comprehensive.py # OPE 示例综合测试
    ├── test_opedefs_comparison.py       # OPEdefs 对比测试
    ├── test_voa_manual_examples.py      # VOA 手册示例测试
    ├── test_jacobi_virasoro.py          # Jacobi 恒等式测试
    └── test_w3_algebra.py               # W3 代数测试
```

---

## 测试分类

### 1. 基础功能测试

#### test_constants.py
**功能**: 测试常量算符（One, Zero）的基本性质
**覆盖**:
- One 算符的定义和性质
- Zero 算符的定义和性质
- 常量算符的 OPE 规则

**运行**: `pytest tests/test_constants.py -v`

#### test_operator.py
**功能**: 测试算符的基本功能
**覆盖**:
- BasisOperator 创建
- 算符属性（名称、权重）
- 算符相等性判断

**运行**: `pytest tests/test_operator.py -v`

#### test_ope_data.py
**功能**: 测试 OPEData 数据结构
**覆盖**:
- OPEData 创建和访问
- 极点（pole）操作
- OPEData 算术运算

**运行**: `pytest tests/test_ope_data.py -v`

#### test_api.py
**功能**: 测试核心 API 接口
**覆盖**:
- OPE 定义和查询
- bracket 函数
- 基本 OPE 计算

**运行**: `pytest tests/test_api.py -v`

---

### 2. 局域算符测试

#### test_local_operator.py
**功能**: 测试局域算符的性质
**覆盖**:
- 局域算符定义
- Parity（统计性）
- 共形权重

**运行**: `pytest tests/test_local_operator.py -v`

---

### 3. 导数和复合算符测试

#### test_derivative_ope_fix.py
**功能**: 测试导数算符的 OPE 计算修复
**覆盖**:
- 导数算符 OPE 规则
- 左右导数的正确处理
- Bug 修复验证

**运行**: `pytest tests/test_derivative_ope_fix.py -v`

**参考**: Git commit `59cfecb` - 修复左右两侧导数算符 OPE 计算的关键 bug

#### test_composite_left_ope.py
**功能**: 测试复合算符的左 OPE 计算
**覆盖**:
- 复合算符定义
- 左 OPE 计算规则
- 嵌套算符处理

**运行**: `pytest tests/test_composite_left_ope.py -v`

---

### 4. 共形权重测试

#### test_conformal_weight.py
**功能**: 测试共形权重的计算和验证
**覆盖**:
- 基础算符的共形权重
- 导数算符的权重增加
- 复合算符的权重计算
- OPE 中的权重守恒

**运行**: `pytest tests/test_conformal_weight.py -v`

**关键测试**:
- 导数算符: `h(∂A) = h(A) + 1`
- OPE 权重守恒: `h(A) + h(B) = h(C) + n`

---

### 5. 简化函数测试

#### test_simplify.py
**功能**: 测试表达式简化功能
**覆盖**:
- `simplify()` 函数
- 符号表达式简化
- OPE 结果简化

**运行**: `pytest tests/test_simplify.py -v`

**参考**: Git commit `1bbc5c9` - 添加可选的 simplify() 函数

---

### 6. Thielemans 方程测试

#### test_thielemans_eqs.py
**功能**: 测试 Thielemans 论文中的基本方程
**覆盖**:
- 基本 OPE 规则
- Thielemans 论文中的示例

**运行**: `pytest tests/test_thielemans_eqs.py -v`

#### test_thielemans_eqs_extended.py
**功能**: 扩展的 Thielemans 方程测试
**覆盖**:
- 更复杂的 OPE 计算
- 多重嵌套 OPE
- 边界情况

**运行**: `pytest tests/test_thielemans_eqs_extended.py -v`

---

### 7. Jacobi 恒等式测试

#### test_jacobi_virasoro.py ⭐
**功能**: 测试 Jacobi 恒等式在 Virasoro 代数上的应用
**覆盖**:
- Virasoro OPE 结构验证
- `check_jacobi_identity(T, T, T)` 计算
- `verify_jacobi_identity(T, T, T)` 验证
- 与 Mathematica 参考对比
- 简单流代数测试

**运行**:
```bash
# pytest 方式
pytest tests/test_jacobi_virasoro.py -v

# 独立运行
python tests/test_jacobi_virasoro.py
```

**参考测试**: `ref_jacobi_virasoro.wls`

**文档**:
- `JACOBI_IDENTITY_TEST.md` - 测试说明
- `JACOBI_TEST_REPORT.md` - 详细报告
- `JACOBI_TEST_SUMMARY.md` - 快速总结

**状态**: ✅ 所有测试通过（6/6）

**参考**:
- OPEdefs.m 第 1601-1637 行
- Git commit `6703739` - 修复 Jacobi 恒等式实现

---

### 8. W3 代数测试

#### test_w3_algebra.py
**功能**: 测试 W3 代数的 OPE 计算
**覆盖**:
- W3 代数定义
- W 算符的 OPE
- W3 代数的 Jacobi 恒等式

**运行**: `pytest tests/test_w3_algebra.py -v`

**参考测试**: `w3_algebra_test.wls`

**文档**: `W3_ALGEBRA_TESTS.md`

---

### 9. 综合测试

#### test_advanced_ope.py
**功能**: 高级 OPE 计算测试
**覆盖**:
- 复杂 OPE 结构
- 多重嵌套 OPE
- 高阶极点

**运行**: `pytest tests/test_advanced_ope.py -v`

#### test_ope_examples_comprehensive.py
**功能**: OPE 示例的综合测试
**覆盖**:
- OPEdefs.m 中的示例
- ope-examples.nb 中的计算
- 各种代数结构

**运行**: `pytest tests/test_ope_examples_comprehensive.py -v`

#### test_opedefs_comparison.py
**功能**: 与 OPEdefs.m 的系统对比测试
**覆盖**:
- 逐函数对比
- 数值精度验证
- 边界情况对比

**运行**: `pytest tests/test_opedefs_comparison.py -v`

#### test_voa_manual_examples.py
**功能**: VOA 手册中的示例测试
**覆盖**:
- voa-manual.md 中的示例
- 理论验证
- 教学示例

**运行**: `pytest tests/test_voa_manual_examples.py -v`

---

## Mathematica 参考测试

### ref_jacobi_virasoro.wls
**功能**: Jacobi 恒等式的 Mathematica 参考实现
**用途**: 为 `test_jacobi_virasoro.py` 提供参考结果

**运行**: `wolframscript tests/ref_jacobi_virasoro.wls`

**输出**: 5×5 全零矩阵（Jacobi 恒等式成立）

### w3_algebra_test.wls
**功能**: W3 代数的 Mathematica 参考实现
**用途**: 为 `test_w3_algebra.py` 提供参考结果

**运行**: `wolframscript tests/w3_algebra_test.wls`

---

## 运行所有测试

### 运行全部测试
```bash
cd /Users/lelouch/pyope
pytest tests/ -v
```

### 运行特定类别
```bash
# 基础功能测试
pytest tests/test_constants.py tests/test_operator.py tests/test_ope_data.py tests/test_api.py -v

# Jacobi 恒等式测试
pytest tests/test_jacobi_virasoro.py -v

# W3 代数测试
pytest tests/test_w3_algebra.py -v

# 综合测试
pytest tests/test_advanced_ope.py tests/test_ope_examples_comprehensive.py -v
```

### 运行并生成覆盖率报告
```bash
pytest tests/ --cov=src/pyope --cov-report=html
```

---

## 测试统计

### 测试文件数量
- Python 测试文件: 18 个
- Mathematica 参考文件: 2 个
- 文档文件: 5 个

### 测试覆盖的代数结构
- ✅ Virasoro 代数
- ✅ 简单流代数（current algebra）
- ✅ W3 代数
- ⚠️ 待测试: 超对称代数、费米算符

### 测试覆盖的功能模块
- ✅ 算符定义和操作
- ✅ OPE 计算
- ✅ 导数算符
- ✅ 复合算符
- ✅ 共形权重
- ✅ Jacobi 恒等式
- ✅ 表达式简化
- ✅ 与 OPEdefs.m 对比

---

## 测试质量指标

### 代码覆盖率
- 目标: > 80%
- 当前: 待测量

### 测试通过率
- Jacobi 恒等式: 100% (6/6)
- 其他测试: 待统计

### 与 Mathematica 一致性
- Jacobi 恒等式: ✅ 100% 一致
- W3 代数: 待验证
- 其他: 待验证

---

## 测试开发指南

### 创建新测试的步骤

1. **确定测试对象**
   - 选择要测试的功能或代数结构
   - 查找 OPEdefs.m 中的对应实现

2. **创建 Mathematica 参考测试**
   - 文件命名: `ref_<feature>.wls`
   - 使用 OPEdefs.m 计算参考结果
   - 输出结果到 JSON 或直接打印

3. **创建 Python 测试**
   - 文件命名: `test_<feature>.py`
   - 使用 pytest 框架
   - 与 Mathematica 参考对比

4. **编写测试文档**
   - 说明测试目的
   - 列出测试覆盖
   - 提供运行示例

### 测试命名规范

- **测试文件**: `test_<feature>.py`
- **参考文件**: `ref_<feature>.wls`
- **测试类**: `TestFeatureName`
- **测试方法**: `test_specific_case`

### 测试结构模板

```python
"""
<Feature> 测试

测试内容：
1. ...
2. ...

参考：
- OPEdefs.m 第 X-Y 行
- ref_<feature>.wls
"""

import pytest
from pyope import ...

class TestFeature:
    """测试 <Feature> 功能"""

    def setup_method(self):
        """设置测试环境"""
        # 初始化算符、参数等
        pass

    def test_basic_case(self):
        """测试基本情况"""
        # 测试代码
        pass

    def test_comparison_with_mathematica(self):
        """与 Mathematica 参考对比"""
        # 对比代码
        pass
```

---

## 持续集成

### GitHub Actions（待配置）
```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
      - name: Install dependencies
        run: pip install -e .[dev]
      - name: Run tests
        run: pytest tests/ -v --cov=src/pyope
```

---

## 参考资料

### 代码参考
- **OPEdefs.m**: Kris Thielemans 的 Mathematica 包
- **voa-manual.md**: VOA 实现手册
- **ope-examples.nb**: OPE 计算示例

### 理论参考
- Di Francesco et al., "Conformal Field Theory" (1997)
- V. Kac, "Vertex Algebras for Beginners" (1998)
- Thielemans, K., "A Mathematica package for computing operator product expansions" (1991)

---

## 未来工作

### 待添加的测试
- [ ] 超对称代数（N=1, N=2）
- [ ] 费米算符的 Jacobi 恒等式
- [ ] 混合玻色-费米算符
- [ ] W_N 代数（N > 3）
- [ ] 性能基准测试

### 待改进的测试
- [ ] 增加代码覆盖率
- [ ] 添加性能测试
- [ ] 添加压力测试（大规模计算）
- [ ] 改进测试文档

---

## 联系与贡献

如需添加新测试或报告测试问题，请：
1. 查看现有测试是否覆盖
2. 参考测试开发指南
3. 创建测试文件和文档
4. 提交 Pull Request

---

**最后更新**: 2026-01-07
**维护者**: Claude (Sonnet 4.5)
**版本**: 1.0
