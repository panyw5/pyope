本项目的目标是开发 `python` 版本的 `OPEdefs.m`，用于计算 vertex operator algebras (VOA) 中的算符积展开 (OPE)。

# 基本功能初步设想 (待完善)

实现算符的四则运算以及导数

- 父类: `Operator`
- VOA 元素类: `LocalOperator`
  - 具有相同 conformal weight 和 parity 的 local operator 能够做加法 `A+B`，减法 `A-B`，结果还是 `LocalOperator`
  - 能够做数乘 $3 * A$ ，结果还是 `LocalOperator`
  - 数乘与加减法能够自动做基本的合并同类项，比如 `A + (-1)*A = Zero`
  - 数乘与加减法自动完成**分配律**，输出 `LocalOperator` 类的实例
  - 能够求导，结果还是 `LocalOperator` 类的实例
  - 输出结果可以在 `Jupyter Notebook` 中以 `LaTex` 渲染，也可以手动调用 `.latex()` 方法得到 `LaTex` 代码，用 `.render()` 输出渲染的结果
- 针对不同的应用场景设置一些 `LocalOperator` 子类，添加更多的方便的性质与方法
  - 特殊的 `LocalOperator`：`One`，`Zero` (`0 * A = Zero`)
  - 比如 `BasicLocalOperator` 代表用户指定的基本的局域算符 $a(z)$，具有特定的 conformal weight $h(a)$ 以及 parity (文献中记为 $|a|$，标记玻色或费米子)，这些性质也应该能够传递到复合算符上
  - 算符的各种复合，比如 normal order product $(AB)$，bracket $\{AB\}$ 等
- 求导：`d(A)` 或者 `dn(n, A)`  (结果还是 `LocalOperator` 类): 计算 VOA 元素 A 的导数 $∂A(z)$ 或者 $∂^n A(z)$。
  应当自动实现如下性质
  - **可加性**：$∂(A+B) = ∂A + ∂B$
  - **线性性**：$∂(cA) = c∂A$
  - **莱布尼茨律**: 比如 $\partial (NO(AB)) = NO(\partial A B) + NO(A \partial B)$，$\{A \partial B\}_{n + 1} = n\{AB\}_n + \partial (\{AB\}_{n + 1})$，$\{\partial A B\}_{n + 1} = -n \{AB\}_n$

# 关键功能构想 (待完善)

- 根据顶点算符代数的解析形式计算设计 `OPE`、`NO`、`bracket`、`d`、`dn` 等函数
  - 参考 `OPEdefs.m` 的实现
  - 充分利用你的 `voa` SKILL

## Operator Product Expansion 与 Bracket

算符积展开:

$$
a(z)b(w) = \sum_{n} \frac{\{ab\}_n(w)}{(z - w)^n}
$$

其中 $\{ab\}_n$ 称为 $a, b$ 的 bracket。

> 注意，在 `VOA-manual.md` 中，用的符号是 $[ab]_n$，其实等价于我用的记号 $\{ab\}_n$

- `OPE[A, B] = MakeOPE(...)` 用于指定 `BasicLocalOperator` `A, B` 的 OPE
  - `MakeOPE` 返回 `OPEData` 类的实例，包含各个极点的信息
- `OPE(A, B)`: 计算两个 VOA 元素 A 和 B 的算符积展开 $A(z)B(w)$，结果是 `OPEData`，包含各个极点的信息
  - `OPEData` 也应该可以进行加减数乘，自动实现分配律和合并同类项，结果还是 `OPEData` 类的实例
  - 从 `OPE(A, B)` 可以提取 `.pole(n)` 得到 $(z-w)^{-n}$ 的系数 $\{AB\}_n$ (`LocalOperator` 类)，可以直接进行后续的运算
  - `OPE(A,B)` 的输出可以自动在 Jupyter notebook 中以 `LaTex` 渲染，也可以手动调用 `.latex()` 方法得到 `LaTex` 代码，用 `.render()` 输出渲染的结果
- `bracket(A, B, n)` (`LocalOperator` 类): 计算 A 和 B 的 $\{AB\}_n(z)$

## 正规乘积 (Normal ordered product)

正规乘积 $NO(AB)(z) = \{AB\}_0(z)$，由 bracket 给出

- `NO(A, B)` (`LocalOperator` 类): 计算 A 和 B 的正则序列化 (normal ordering) $(AB)(z)$。

## 模态

`A.mode(n)` (`Operator`): 获取 `LocalOperator` `A` (代表 $A(z)$) 的第 n 个模态 (mode) $A_n$。

- 模态的加减和数乘依然是模态
- 对易子: `commutator(A.mode(m), B.mode(n))` 计算两个模态的对易子
  $$
  [a_m, b_n] = \sum_{\ell \ge 1} 
  \begin{pmatrix}
  m + h_a - 1\\
  \ell - 1  
  \end{pmatrix} (\{ab\}_{\ell})_{m + n}
  $$
- 若 local operator `a` 的 conformal weight 是 $h(a)$，$a_{n}$ 的 conformal weight 为 $-n$，特别地 $a_{-h(a)}$ 的 conformal weight 是 $h(a)$
- 模态 $a_n$ 对 `LocalOperator` $b(z)$ 的作用: $(a_{-h_a + m} b)(z) = \{ab\}_{m}(z)$， 或等价地 $(a_{m} b)(z) = \{ab\}_{m + h_a}(z)$

## Jacobi 恒等式

$$
\{A\{ BC\}_p\}_q
= (-1)^{|A||B|} \{B \{AC\}_q\}_p
+ \sum_{\ell \ge 1}

\begin{pmatrix}
  q - 1\\
  \ell - 1
\end{pmatrix}  \{\{AB\}_{\ell}C\}_{p + q - \ell}
$$

其中 $|A|, |B|$ 是 $A, B$ 的 parity。

---

# 项目结构规划

## 目标目录结构

```
pyope/
├── pyproject.toml           # 项目配置（PEP 517/518 标准）
├── README.md                # 项目说明文档
├── CLAUDE.md                # 开发指南
├── .gitignore
├── src/
│   └── pyope/              # 源代码包（采用 src layout）
│       ├── __init__.py          # 包初始化，导出公共 API
│       ├── operators.py         # 算符类定义
│       ├── local_operator.py    # 局域算符基类
│       ├── ope_data.py          # OPE 数据结构
│       ├── api.py               # 公共 API (OPE, NO, bracket)
│       ├── registry.py          # OPE 注册表
│       ├── constants.py         # 常数算符 (One, Zero, Delta)
│       └── utils.py             # 工具函数
├── tests/                   # 测试文件
│   ├── __init__.py
│   ├── test_operator.py
│   ├── test_ope_data.py
│   ├── test_api.py
│   ├── test_constants.py
│   ├── test_advanced_ope.py
│   ├── test_local_operator.py
│   └── test_thielemans_eqs.py
├── demo/                    # 演示 Jupyter Notebooks
│   ├── pyope_demo.ipynb
│   └── latex_display_example.ipynb
├── papers/                  # 参考文献
└── OPEdefs/                # Mathematica 原始代码（参考）
```

## 当前状态分析

### 存在的问题

1. **源代码缺失**：`src/` 目录下只有分析文档，没有实际的 Python 源代码
2. **包配置缺失**：缺少 `pyproject.toml` 等现代 Python 项目配置文件
3. **目录结构混乱**：之前的 `pyope/pyope/` 源代码已被删除但未重建
4. **测试文件孤立**：`tests/` 目录中的测试文件引用不存在的 `pyope` 模块

### 需要完成的工作

- 重建 `src/pyope/` 源代码目录
- 创建 `pyproject.toml` 项目配置
- 逐步实现核心功能模块
- 确保测试文件能够正常运行

---

# 分阶段开发计划

## 阶段 0：项目基础设施搭建 ✅

**目标**：建立项目基本结构和开发环境

### 任务清单

- [X] 创建 `pyproject.toml` 配置文件
  - 定义项目元数据（名称、版本、作者等）
  - 配置依赖项（sympy, numpy, IPython 等）
  - 设置构建系统（setuptools 或 hatchling）
  - 配置开发工具（pytest, black, mypy 等）
- [X] 创建 `src/pyope/` 目录结构
- [X] 创建基础的 `__init__.py` 文件
- [X] 验证包可以正常安装（`pip install -e .`）
- [X] 确保现有测试文件能够被发现（即使暂时失败）

### 参考资料

- Python Packaging User Guide
- 现有的 `tests/` 目录中的测试文件

---

## 阶段 1：核心算符类实现 ✅

**目标**：实现基础的算符类层次结构

### 1.1 基础算符类 (`operators.py`) ✅

**实现内容**：

- [X] `Operator` 基类
  - 定义算符的基本接口
  - 实现 `__repr__` 和 `__str__` 方法
  - 实现 `__eq__` 和 `__hash__` 方法
- [X] `BasisOperator` 类（基本算符）
  - 存储算符名称、conformal weight、parity
  - 实现 `is_bosonic` 和 `is_fermionic` 属性
  - 实现 LaTeX 渲染方法
- [X] `DerivativeOperator` 类（导数算符）
  - 表示 $\partial^n A(z)$ 形式的导数
  - 实现导数的嵌套和简化
- [X] `NormalOrderedOperator` 类（正规序算符）
  - 表示 $NO(AB)$ 形式的正规序乘积
  - 实现基本的代数性质

**测试**：

- 运行 `tests/test_operator.py` 中的相关测试
- 创建对应的 Mathematica 测试脚本进行对比验证

**参考资料**：

- `OPEdefs/OPEdefs.m` 中的算符定义
- `src/OPEdefs_Analysis.md` 中的分析
- `voa` SKILL 中的 `voa-manual.md`

### 1.2 局域算符类 (`local_operator.py`) ✅

**实现内容**：

- [X] `LocalOperator` 基类
  - 表示 VOA 中的局域算符 $A(z)$
  - 实现加法、减法、数乘运算
  - 实现自动合并同类项
  - 实现分配律
- [X] `OperatorSum` 类（算符和）
  - 表示多个算符的线性组合
  - 实现自动简化和合并
- [X] `OperatorProduct` 类（算符积）
  - 表示算符的乘积形式
  - 为后续 OPE 计算做准备

**测试**：

- 运行 `tests/test_local_operator.py` 中的相关测试
- 验证加减法、数乘的正确性
- 验证合并同类项功能

**参考资料**：

- `OPEdefs/OPEdefs.m` 中的 `LocalOperator` 相关代码

### 1.3 常数算符 (`constants.py`) ✅

**实现内容**：

- [X] `ConstantOperator` 基类
- [X] `One` 单位算符
- [X] `Zero` 零算符
- [X] `Delta` 函数（如果需要）

**测试**：

- [X] 运行 `tests/test_constants.py`
- [X] 验证常数算符的特殊性质

---

## 阶段 2：导数运算实现 ✅

**目标**：实现算符的导数运算及其代数性质

### 2.1 导数函数 (`api.py` 部分) ✅

**实现内容**：

- [X] `d(A)` 函数：计算 $\partial A(z)$
- [X] `dn(n, A)` 函数：计算 $\partial^n A(z)$
- [X] 实现导数的可加性：$\partial(A+B) = \partial A + \partial B$
- [X] 实现导数的线性性：$\partial(cA) = c\partial A$
- [X] 实现莱布尼茨律（针对不同类型的算符积）

**测试**：

- 创建导数运算的测试用例
- 与 Mathematica 版本对比验证

**参考资料**：

- `OPEdefs/OPEdefs.m` 中的导数实现
- `voa-manual.md` 中的导数性质

---

## 阶段 3：OPE 数据结构与注册表 ✅

**目标**：实现 OPE 数据的存储和管理机制

### 3.1 OPE 数据结构 (`ope_data.py`) ✅

**实现内容**：

- [X] `OPEData` 类
  - 存储 OPE 的极点信息：$\{AB\}_n$ 对应的系数
  - 实现 `.pole(n)` 方法提取特定极点的系数
  - 实现加法、减法、数乘运算
  - 实现自动合并同类项
  - 实现 LaTeX 渲染
- [X] OPE 的代数运算
  - 分配律的自动展开
  - 极点的自动合并

**测试**：

- [X] 运行 `tests/test_ope_data.py`
- [X] 验证 OPE 数据的存储和提取

### 3.2 OPE 注册表 (`registry.py`) ✅

**实现内容**：

- [X] `OPERegistry` 类
  - 存储用户定义的基本算符的 OPE
  - 提供注册和查询接口
  - 实现 `Bosonic` 和 `Fermionic` 辅助函数
- [X] 全局注册表 `ope_registry`
- [X] OPE 定义的装饰器或函数接口

**测试**：

- 验证 OPE 的注册和查询功能
- 测试多个 OPE 的管理

**参考资料**：

- `OPEdefs/OPEdefs.m` 中的 OPE 存储机制

---

## 阶段 4：核心 OPE 计算实现 ✅

**目标**：实现 OPE、bracket、NO 等核心计算功能

### 4.1 基本 OPE 计算 (`api.py`) ✅

**实现内容**：

- [X] `OPE(A, B)` 函数
  - 对于基本算符：从注册表查询
  - 对于复合算符：递归计算
  - 实现导数的 OPE 计算规则
  - 实现正规序算符的 OPE 计算规则
- [X] `bracket(A, B, n)` 函数
  - 计算 $\{AB\}_n(z)$
  - 从 OPE 中提取对应极点
- [X] `NO(A, B)` 函数
  - 计算正规序乘积 $(AB)(z) = \{AB\}_0(z)$

**测试**：

- [X] 运行 `tests/test_api.py`
- [X] 运行 `tests/test_composite_left_ope.py` (8/8 通过)
- [X] 运行 `tests/test_voa_manual_examples.py` (8/8 通过)
- [X] 与 Mathematica 版本对比验证简单算例

**参考资料**：

- `OPEdefs/OPEdefs.m` 中的 OPE 计算逻辑
- `voa-manual.md` Section 3.3 Implementation
- `papers/` 中的 Thielemans 论文

### 4.2 高级 OPE 计算规则 ✅

**实现内容**：

- [X] 导数算符的 OPE 计算
  - $\partial A(z) \cdot B(w)$ 的计算规则
  - $A(z) \cdot \partial B(w)$ 的计算规则
- [X] 正规序算符的 OPE 计算
  - $NO(AB)(z) \cdot C(w)$ 的计算规则（左侧和右侧）
  - 完整的 Jacobi 恒等式实现
- [X] 复合算符的递归展开
- [X] 自动应用莱布尼茨律

**测试**：

- [X] 创建复杂算例进行测试
- [X] 验证 Jacobi 恒等式（部分）
- [X] 测试嵌套正规序算符

---

## 阶段 5：Thielemans 方程验证 ✅

**目标**：验证实现的正确性，通过 Thielemans 论文中的方程

### 5.1 实现论文中的测试用例 ✅

**实现内容**：

- [X] 实现 Virasoro 代数的 OPE
- [X] 实现 $W_3$ 代数的 OPE
- [X] 验证论文中的关键方程
  - [X] eq 3.3.1: [∂A, B]_q = -(q-1)[A, B]_{q-1}
  - [X] eq 3.3.2: [A, ∂B]_q = (q-1)[A, B]_{q-1} + ∂[A, B]_q
  - [X] eq 3.3.4: Jacobi 恒等式（复合算符 OPE）

**测试**：

- ✅ 运行 `tests/test_thielemans_eqs.py` - 全部通过
- ✅ 运行 `tests/test_thielemans_eqs_extended.py` - 7/7 通过
- ✅ 验证 Virasoro 代数的导数规则
- ✅ 验证 Kac-Moody 代数的导数规则
- ✅ 验证复合算符的 OPE 计算

**参考资料**：

- `papers/[Thielemans] An Algorithmic Approach...`
- `OPEdefs/` 中的示例 notebook

---

## 阶段 6：模态运算（可选）

**目标**：实现算符的模态表示和对易子计算

### 6.1 模态类实现

**实现内容**：

- [ ] `Mode` 类表示 $A_n$
- [ ] `.mode(n)` 方法从 `LocalOperator` 提取模态
- [ ] 模态的加减数乘运算
- [ ] `commutator(A.mode(m), B.mode(n))` 对易子计算
- [ ] 模态对局域算符的作用

**测试**：

- 创建模态运算的测试用例
- 验证对易关系

**参考资料**：

- `voa-manual.md` 中的模态定义
- README.md 中的模态公式

---

## 阶段 7：LaTeX 渲染与可视化 ✅

**目标**：实现美观的 LaTeX 输出和 Jupyter Notebook 集成

### 7.1 LaTeX 渲染 ✅

**实现内容**：

- [X] 为所有算符类实现 `_latex()` 方法
  - [X] `BasisOperator` - 基本算符
  - [X] `DerivativeOperator` - 导数算符（$\partial A$, $\partial^n A$）
  - [X] `NormalOrderedOperator` - 正规序算符（$:AB:$）
  - [X] `ConstantOperator` - 常数算符（One → 1, Zero → 0）
- [X] 为 `OPEData` 实现 `_repr_latex_()` 方法用于 Jupyter 自动渲染
- [X] OPE 结果自动以 LaTeX 格式显示

**测试**：

- ✅ 创建了 `demo/latex_rendering_demo.ipynb` 演示渲染效果
- ✅ 验证了各种复杂表达式的渲染
- ✅ 所有算符和 OPE 结果都能正确渲染

**功能特性**：

- ✅ 导数算符：$\partial T$, $\partial^2 T$
- ✅ 正规序：$(TT)$, $(T\partial T)$
- ✅ 嵌套正规序：$((AB)C)$
- ✅ OPE 分式：$\frac{c/2}{(z-w)^4} + \frac{2T}{(z-w)^2} + \frac{\partial T}{z-w}$
- ✅ 符号参数：支持 sympy 符号的自动渲染

---

## 阶段 8：文档与示例

**目标**：完善文档和使用示例

### 8.1 API 文档

**实现内容**：

- [ ] 为所有公共类和函数添加 docstring
- [ ] 生成 API 参考文档（使用 Sphinx 或 MkDocs）

### 8.2 使用示例

**实现内容**：

- [ ] 完善 `demo/pyope_demo.ipynb`
  - 基本用法示例
  - Virasoro 代数示例
  - $W_3$ 代数示例
- [ ] 创建更多示例 notebook

### 8.3 README 完善

**实现内容**：

- [ ] 添加安装说明
- [ ] 添加快速开始指南
- [ ] 添加使用示例
- [ ] 添加贡献指南

---

## 开发原则

### 测试驱动开发

- 每个功能实现后立即编写测试
- 与 Mathematica 版本对比验证结果
- 使用 `pytest` 进行自动化测试

### 渐进式开发

- 先实现核心功能，再添加高级特性
- 每个阶段完成后进行充分测试
- 确保每个阶段的代码质量

### 代码质量

- 遵循 PEP 8 代码风格
- 使用类型注解（Type Hints）
- 编写清晰的中文注释
- 保持代码简洁和可读性

### 协作开发

- 对于大型文件（如 `.nb` 文件）的分析，使用 `gemini` 协助
- 对于复杂的数学逻辑，参考 `voa` SKILL
- 充分利用 `OPEdefs.m` 作为参考实现

---

## 当前优先级

**已完成**：

- ✅ 阶段 0：项目基础设施搭建
- ✅ 阶段 1：核心算符类实现
- ✅ 阶段 2：导数运算实现
- ✅ 阶段 3：OPE 数据结构与注册表
- ✅ 阶段 4：核心 OPE 计算实现（包括高级规则）

**测试状态**：

- ✅ 137/137 测试全部通过 (重大更新！)
  - `tests/test_composite_left_ope.py`: 8/8 ✅
  - `tests/test_voa_manual_examples.py`: 8/8 ✅
  - `tests/test_ope_examples_comprehensive.py`: 17/17 ✅
  - `tests/test_simplify.py`: 14/14 ✅
  - `tests/test_w3_algebra.py`: 10/10 ✅
  - `tests/test_advanced_ope.py`: 17/17 ✅ (修复对易测试)
  - 以及其他所有测试模块

**重大修复 (2026-01-07)**：

- ✅ **修复 Jacobi 恒等式实现** - 在 `_ope_composite_left` 和 `_ope_composite_right` 中添加了缺失的第三项
- ✅ **OPE[T, NO(T,T)] 极点修正** - 现在正确返回 max_pole=6（之前错误地返回 4）
- ✅ **Sugawara 张量计算修正** - OPE[TSugawara, TSugawara] 现在正确包含 4 阶极点 2*One
- ✅ **算符对易关系修复** - 修复了 `test_bosonic_commutation` 和 `test_fermionic_anticommutation` 测试
- ✅ **W₃ 代数验证** - 所有 W₃ 代数测试与 Mathematica 结果一致

**最新功能**：

- ✅ OPE 对称性自动处理 - `OPE(B,A)` 自动使用对易公式从 `OPE(A,B)` 计算
- ✅ 算符排序机制 - 按声明顺序自动排列
- ✅ 可选的化简功能 - `simplify(expr)` 函数
- ✅ W₃ 代数完整支持 - 包括 T-T, T-W, W-W OPE 和辅助算符 Λ
- ✅ 正确的复合算符 OPE - 完整实现 Jacobi 恒等式的三项求和
- ✅ Thielemans 方程验证 - 验证了论文中的关键导数规则和 Jacobi 恒等式
- ✅ **LaTeX 自动渲染** - OPE 结果和算符在 Jupyter Notebook 中自动以美观的数学公式显示

**下一步**：

- ✅ 阶段 5：Thielemans 方程验证（已完成）
- ✅ 阶段 7：LaTeX 渲染与可视化（已完成）
- 🔄 阶段 9：Null States 计算（进行中 - 阶段 1 完成）
- 阶段 6：模态运算（可选）
- 阶段 8：文档与示例（进行中）

建议按照阶段顺序逐步推进，确保每个阶段完成并测试通过后再进入下一阶段。
