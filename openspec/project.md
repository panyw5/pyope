# Project Context

## Purpose

本项目旨在将 Mathematica 代码包 `OPEdefs.m` 重构为 Python 库 `pyope`，用于进行**顶点算符代数 (Vertex Operator Algebra, VOA)** 的符号计算。

### 核心功能

- **算符四则运算与导数**: 实现 VOA 元素的加、减、数乘和求导运算
- **算符积展开 (OPE)**: 计算算符的奇异部分展开
- **Bracket 与正规序**: 计算 $\{AB\}_n$ 和正规序乘积 $(AB)$
- **LaTeX 渲染**: 在 Jupyter Notebook 中自动渲染数学公式
- **对称性处理**: 自动处理玻色子和费米子的对易/反对易关系

### 目标用户

- 理论物理研究人员（共形场论、弦论）
- 数学研究人员（顶点算符代数）
- 需要进行 VOA 符号计算的教育工作者

## Tech Stack

### 核心语言与框架
- **Python 3.10+**: 主要开发语言
- **SymPy**: 符号计算核心库
- **NumPy**: 数值计算支持
- **IPython.display**: Jupyter Notebook LaTeX 渲染

### 构建工具
- **pyproject.toml**: 现代 Python 项目配置（PEP 517/518）
- **setuptools**: 构建后端
- **pytest**: 测试框架
- **pytest-cov**: 代码覆盖率

### 开发工具
- **black**: 代码格式化
- **mypy**: 类型检查
- **ruff**: 快速 Linter

### 参考实现
- **Mathematica**: `OPEdefs.m` (Wolfram Language)
- **Wolfram Script**: 用于验证计算正确性

### 项目结构
```
src/pyope/          # 源代码包
├── operators.py       # 算符类定义
├── local_operator.py  # 局域算符基类
├── ope_data.py        # OPE 数据结构
├── api.py             # 公共 API (OPE, NO, bracket, d, dn)
├── registry.py        # OPE 注册表
├── constants.py       # 常数算符 (One, Zero)
└── utils.py           # 工具函数

tests/             # 测试文件
demo/              # Jupyter Notebook 演示
papers/            # VOA 参考文献
OPEdefs/           # Mathematica 原始代码（参考）
```

## Project Conventions

### Code Style

#### 命名规范
- **类名**: `PascalCase` (例: `LocalOperator`, `OPEData`)
- **函数名**: `snake_case` (例: `normal_ordered`, `bracket`)
- **常量**: `UPPER_SNAKE_CASE` (例: `One`, `Zero`, `ope_registry`)
- **私有方法**: `_leading_underscore` (例: `_ope_composite_left`)

#### 代码风格
- **精简高效**: 代码风格始终定位为精简高效、毫无冗余
- **注释原则**: 非必要不形成 - 只在逻辑不显而易见时添加注释
- **类型注解**: 所有公共 API 必须包含类型提示
- **中文注释**: 允许使用中文注释解释数学概念

#### 格式化工具
- **Black**: 默认配置（line-length=88）
- **Ruff**: 快速 linting

### Architecture Patterns

#### 算符类层次结构
```
Operator (抽象基类)
├── BasicOperator          # 基本算符 a(z)
├── DerivativeOperator     # 导数算符 ∂A
├── NormalOrderedOperator  # 正规序算符 (AB)
└── ConstantOperator       # 常数算符 (One, Zero)

LocalOperator (Operator 子类)
├── OperatorSum            # 线性组合 A + 2B - C
└── OperatorProduct        # 算符乘积形式
```

#### 设计模式
- **注册表模式** (`registry.py`): 管理基本算符的 OPE 定义
- **访问者模式**: 用于 OPE 计算的递归遍历
- **表达式树**: 算符组成树状结构，便于简化

#### 关键设计决策
1. **不可变性**: 算符对象一旦创建不应被修改（返回新对象）
2. **自动简化**: 所有运算后自动调用 `simplify()` 进行合并同类项
3. **延迟计算**: OPE 计算在需要时才执行（按需求值）
4. **类型安全**: 使用 `@abstractmethod` 确保接口完整性

### Testing Strategy

#### 测试分类
1. **单元测试** (`tests/test_*.py`): 测试单个类/函数
2. **集成测试** (`tests/test_*_examples.py`): 测试完整计算流程
3. **对比验证** (`*_reference.py`): 与 Mathematica 结果对比

#### 测试文件组织
- `test_operator.py`: 算符类基础测试
- `test_local_operator.py`: 局域算符运算测试
- `test_ope_data.py`: OPE 数据结构测试
- `test_api.py`: 公共 API 测试
- `test_thielemans_eqs.py`: Thielemans 论文方程验证
- `test_voa_manual_examples.py`: VOA manual 示例测试

#### 测试驱动开发 (TDD)
- 每个功能实现前先编写测试
- 测试用例参考:
  - `OPEdefs/` 中的 Mathematica 示例
  - `papers/` 中的论文公式
  - `voa` SKILL 中的示例

#### 对比验证策略
每个测试都需要两个文件：
1. `test_*.py`: Python (pyope) 实现
2. `test_*.wls`: Wolfram Script 参考实现
3. 对比两个结果的等价性

#### 覆盖率目标
- 核心模块: ≥ 90%
- API 层: 100%
- 整体: ≥ 80%

### Git Workflow

#### 分支策略
- `main`: 稳定发布版本
- `feature/*`: 功能开发分支
- `fix/*`: Bug 修复分支

#### Commit 规范
使用 Conventional Commits:
- `feat: 添加新功能`
- `fix: 修复 bug`
- `refactor: 代码重构`
- `test: 添加测试`
- `docs: 文档更新`
- `chore: 构建/工具更新`

#### 分支命名
- `feature/operator-multiplication`: 添加算符乘法功能
- `fix/normal-ordering-simplify`: 修复正规序简化 bug
- `refactor/ope-registry`: 重构 OPE 注册表

#### PR 流程
1. 创建 feature 分支
2. 编写代码和测试
3. 运行 `pytest` 确保测试通过
4. 创建 PR 并描述变更
5. Code review
6. 合并到 main

## Domain Context

### 顶点算符代数 (VOA) 基础

#### 核心概念

**局域算符** (Local Operator): $A(z)$
- 具有 **conformal weight** (共形权重) $h(A)$
- 具有 **parity** (宇称) $|A|$ (玻色子/费米子)
- 满足 **定域性** (locality): $(z-w)^N [A(z), B(w)] = 0$ 对于足够大的 $N$

**算符积展开** (Operator Product Expansion):
$$
A(z)B(w) = \sum_{n=0}^{\infty} \frac{\{AB\}_n(w)}{(z-w)^n}
$$
其中 $\{AB\}_n$ 称为 **bracket**。

**正规序乘积** (Normal Ordered Product):
$$
(AB)(z) = \{AB\}_0(z) = \oint_z \frac{dw}{2\pi i} \frac{A(w)B(z)}{w-z}
$$

**导数算符**:
$$
\partial A(z) = \sum_n a_{n-h-1} z^{-n-1}
$$

#### 关键恒等式

**Jacobi 恒等式** (VOA 的核心):
$$
\{A\{BC\}_p\}_q = (-1)^{|A||B|} \{B\{AC\}_q\}_p + \sum_{\ell \ge 1} \binom{q-1}{\ell-1} \{\{AB\}_{\ell}C\}_{p+q-\ell}
$$

**Thielemans 方程** (导数规则):
- eq 3.3.1: $[\partial A, B]_q = -(q-1)[AB]_{q-1}$
- eq 3.3.2: $[A, \partial B]_q = (q-1)[AB]_{q-1} + \partial[AB]_q$

#### 对易关系

**玻色子** (Bosonic): $|A| = 0$
$$
[A(z), B(w)] = \sum_{n \ge 0} \frac{\{AB\}_n(w)}{(z-w)^n} - \sum_{n \ge 0} \frac{(-1)^{|A||B|}\{BA\}_n(w)}{(w-z)^n}
$$

**费米子** (Fermionic): $|A| = 1$
满足反对易关系

#### 典型 VOA 示例

**Virasoro 代数**:
- 能量动量张量 $T(z)$ with $h(T) = 2$
- OPE: $T(z)T(w) = \frac{c/2}{(z-w)^4} + \frac{2T(w)}{(z-w)^2} + \frac{\partial T(w)}{z-w}$
- 其中 $c$ 是中心荷

**Kac-Moody 代数**:
- 流算符 $J^a(z)$ with $h(J) = 1$
- OPE: $J^a(z)J^b(w) = \frac{k \delta^{ab}}{(z-w)^2} + \frac{f^{ab}_c J^c(w)}{z-w}$

**$W_3$ 代数**:
- 主场 $T(z)$ with $h(T) = 2$
- $W(z)$ with $h(W) = 3$
- 非线性 OPE 结构

### 符号约定

- 算符: 大写字母 $A, B, C, \ldots$
- 基本算符: 小写字母 $a, b, c, \ldots$
- 导数: $\partial = \frac{d}{dz}$
- 模态: $A_n$ (第 $n$ 个模)
- Bracket: $\{AB\}_n$
- 正规序: $(AB)$ 或 $:AB:$

### 参考文献优先级

1. **[Thielemans] An Algorithmic Approach to ...**: 主要参考
2. **OPEdefs.m**: Mathematica 实现参考
3. **voa-manual.md**: VOA 计算手册
4. **Frenkel-Lepowsky-Meurman**: VOA 理论基础

## Important Constraints

### 技术约束

1. **与 Mathematica 结果一致**: 所有计算必须与 `OPEdefs.m` 结果完全一致
2. **符号计算精度**: 使用 SymPy 的精确符号计算，避免浮点数误差
3. **不可变对象**: 算符对象一旦创建不应修改（函数式编程风格）
4. **类型安全**: 所有公共 API 必须有类型注解

### 性能约束

1. **OPE 计算复杂度**: 复合算符的 OPE 计算是指数级的，需要优化
2. **内存管理**: 表达式树可能很深，需要考虑内存使用
3. **简化算法**: `simplify()` 函数需要高效（避免重复简化）

### 兼容性约束

1. **Python 版本**: 支持 Python 3.10+
2. **SymPy 版本**: ≥ 1.12
3. **Jupyter Notebook**: 6.x 和 7.x

### 业务约束

1. **数学正确性**: 任何变更都不能破坏现有测试（226/228 通过）
2. **API 稳定性**: 公共 API 一旦发布不应轻易更改
3. **文档同步**: 代码变更时必须同步更新文档

## External Dependencies

### Python 包依赖

#### 运行时依赖
- **sympy** (≥1.12): 符号计算核心
  - 用于符号表达式、求导、简化
  - 符号参数（如 $c$, $k$）表示
- **numpy** (≥1.24): 数值计算支持
- **IPython** (≥7.0): Jupyter Notebook LaTeX 渲染

#### 开发依赖
- **pytest** (≥7.0): 测试框架
- **pytest-cov**: 覆盖率报告
- **black**: 代码格式化
- **mypy**: 类型检查
- **ruff**: 快速 linter

### 外部系统

#### Wolfram Language (Mathematica)
- **用途**: 验证计算正确性的参考实现
- **调用方式**: `wolframscript` 或 `wolfram-engine`
- **关键文件**: `OPEdefs/OPEdefs.m`

#### 参考文档

- **`voa` SKILL**: Claude Code 技能包，包含 VOA 计算手册
  - `voa-manual.md`: 实现指南
  - `voa-examples.md`: 示例代码
- **`papers/` 目录**: VOA 理论文献
  - Thielemans 论文（主要参考）
  - Frenkel-Lepowsky-Meurman 书籍

### 项目内部依赖

#### 源代码依赖关系
```
api.py (顶层 API)
├── operators.py
├── local_operator.py
├── ope_data.py
├── registry.py
├── constants.py
└── utils.py
```

#### 测试依赖
- 所有测试依赖 `src/pyope/` 模块
- 部分测试依赖 `OPEdefs/` 中的 Mathematica 示例

## 开发状态总览

### 已完成功能 (2026-01-23)

- ✅ **阶段 0**: 项目基础设施
- ✅ **阶段 1**: 核心算符类实现
- ✅ **阶段 2**: 导数运算实现
- ✅ **阶段 3**: OPE 数据结构与注册表
- ✅ **阶段 4**: 核心 OPE 计算（包括高级规则）
- ✅ **阶段 5**: Thielemans 方程验证
- ✅ **阶段 7**: LaTeX 渲染与可视化

### 测试状态

- **总测试数**: 228
- **通过**: 226 (99.1%)
- **失败**: 2 (与当前功能无关)

### 待完成功能

- 🔄 **阶段 9**: Null States 计算（进行中）
- ⏳ **阶段 6**: 模态运算（可选）
- ⏳ **阶段 8**: 完整文档与示例

### 已知限制

1. **嵌套正规序性能**: 深度嵌套的正规序计算较慢
2. **自动简化**: 某些复杂表达式可能无法完全简化
3. **模态运算**: 尚未实现（阶段 6）
