# How to Write Spec

> 让 AI agent 在你的代码库中**可靠、自主地**产出高质量代码的关键：写好 spec。

本文档以 pyope 项目的实际 spec 为蓝本，讲解 Trellis 框架下 spec 的**设计理念、注入机制、写作方法和常见陷阱**。

**所有示例均取自本项目**，标注格式为 `📎 source: .trellis/spec/<path>`，你可以直接打开对照阅读。

---

## 目录

1. [什么是 Spec](#1-什么是-spec)
2. [Spec 的注入链路](#2-spec-的注入链路)
3. [目录组织](#3-目录组织)
4. [写 index.md — 导航入口](#4-写-indexmd--导航入口)
5. [写 guideline 文件 — 规范正文](#5-写-guideline-文件--规范正文)
6. [有效写作的 7 个原则](#6-有效写作的-7-个原则)
7. [反面模式](#7-反面模式)
8. [完整示例拆解](#8-完整示例拆解)

---

## 1. 什么是 Spec

Spec（specification）是一组**写给 AI agent 看**的项目编码规范文档。

它和传统 coding style guide 的区别：

| 对比维度 | 传统 Style Guide | Trellis Spec |
|---------|-----------------|-------------|
| **读者** | 人类开发者 | AI agent（LLM） |
| **加载方式** | 人主动阅读 | Hook **自动注入** agent 上下文 |
| **写作风格** | 散文体、可省略前提 | **结构化**、每条规则必须自包含 |
| **执行保障** | 靠 code review | Agent 在生成代码时实时遵循 |
| **更新频率** | 低频 | 每次发现新 pattern 都应更新 |

**核心思想**：Spec 不是文档，而是 **agent 的运行时配置**。Agent 按 spec 写代码，就像程序按配置文件运行。

---

## 2. Spec 的注入链路

理解 spec 怎么流入 agent，决定你怎么写 spec。

### 2.1 会话启动时

```
session-start.py (SessionStart hook)
      │
      ├─ 遍历 .trellis/spec/*/index.md
      ├─ 把每个 index.md 的内容注入为 <guidelines>
      └─ → 主对话 agent 看到所有 index.md 概览
```

**你写的 index.md 直接出现在每次会话开始的上下文中。** 所以 index.md 必须简短 — 它是导航，不是正文。

### 2.2 任务执行时

```
task.py init-context  →  生成 implement.jsonl / check.jsonl / debug.jsonl
      │                         │
      │  记录需要注入的 spec 文件路径          │
      ↓                         ↓
inject-subagent-context.py (PreToolUse hook)
      │
      ├─ 读取 JSONL 中的所有文件路径
      ├─ 读取每个文件的完整内容
      └─ → 包装进 implement/check/debug agent 的 prompt
```

**Agent 只能看到 JSONL 中列出的文件。** 如果一个 guideline 文件没有被列入 JSONL，agent 根本不知道它的存在。

### 2.3 关键推论

| 推论 | 对写作的影响 |
|------|-----------|
| index.md **每次会话都注入** | 保持简短（< 60 行），只放导航和最重要的约束 |
| guideline 文件**按需注入** | 可以写得详细，包含完整代码示例 |
| Agent **不会跟链接** | index.md 中的链接是给人看的；agent 依赖 JSONL 指定要读哪些文件 |
| Agent 读到的是**纯文本** | Markdown 渲染无效；用代码块和缩进来组织，不要依赖排版 |

---

## 3. 目录组织

### 3.1 推荐结构

```
.trellis/spec/
├── <topic-1>/           # 按主题分目录
│   ├── index.md         # 该主题的导航入口（必须）
│   ├── guideline-a.md   # 具体规范 A
│   └── guideline-b.md   # 具体规范 B
├── <topic-2>/
│   ├── index.md
│   └── ...
└── guides/              # 通用思维指南（非特定主题）
    ├── index.md
    └── ...
```

### 3.2 主题拆分原则

按**使用场景**分，不按技术模块分。以本项目为例：

> 📎 source: `.trellis/spec/` 目录结构

```
# 本项目的实际组织（按使用场景）      # 反面：按模块（不推荐）
spec/                                 spec/
├── core/       ← 写库代码时要读        ├── operators/
├── testing/    ← 写测试时要读           ├── api/
├── domain/     ← 理解数学背景时要读     ├── simplify/
└── guides/     ← 跨层开发时要读         └── ...（太碎了）
```

**为什么**：Agent 被分配的任务是"写一个测试"或"实现一个 API"，而不是"改 operators.py"。按场景分，JSONL 配置更自然。

### 3.3 文件命名

- `index.md` — 固定名，session-start.py 依赖它
- guideline 文件：`<topic>-<aspect>.md`，如 `api-conventions.md`、`test-conventions.md`
- 不要用数字编号（`01-xxx.md`）— 文件名本身就是描述

---

## 4. 写 index.md — 导航入口

index.md 有**三个职责**：

1. **30 秒概览** — 让 agent 快速理解这个主题的范围
2. **文件索引** — 列出所有 guideline 文件及其用途
3. **读取顺序指引** — 告诉 agent 先读什么、什么情况下读什么

### 4.1 模板

```markdown
# <主题名称>

> 一句话描述该主题的范围和目的。

---

## Overview

<2-3 句话说明该主题的上下文、关键约束。>

**Key constraint**: <最重要的一条规则，粗体高亮。>

---

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Name](./file.md) | 一句话描述 | Filled / Draft / TODO |

---

## Pre-Development Checklist

Before <doing what>:

1. **Always**: [File A](./a.md) — 为什么
2. **If <condition>**: [File B](./b.md) — 为什么
3. **Before <milestone>**: [File C](./c.md) — 为什么

---

## Quick Reference: <某种快速查阅表>

<一个表格或代码块，让 agent 不用打开子文件就能回答常见问题>
```

### 4.2 实际项目示例

> 📎 source: `.trellis/spec/core/index.md`（62 行，完整文件可直接查看）

来自 `core/index.md`：

```markdown
# Core Library Guidelines

> Architecture, API conventions, and coding patterns for pyope.

## Overview

pyope is a Python library for computing OPE in VOA.
It is a port of Thielemans' OPEdefs.m to Python/SymPy.

**Key constraint**: All symbolic algebra is built on SymPy.
Operators are sympy.Symbol subclasses.

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Architecture](./architecture.md) | Module organization, type hierarchy, dependency flow | Filled |
| [API Conventions](./api-conventions.md) | OPE/NO/bracket/MakeOPE usage patterns | Filled |
| ...

## Pre-Development Checklist

1. **Always**: Architecture — understand module boundaries
2. **Always**: API Conventions — know the public API surface
3. **If touching SymPy expressions**: SymPy Patterns
4. **Before PR/review**: Quality Guidelines

## Quick Reference: Module Map

src/pyope/
├── operators.py   # Type system
├── api.py         # Core API: OPE(), NO(), bracket()
...
```

### 4.3 常见错误

| 错误 | 后果 |
|------|------|
| index.md 写了 100+ 行详细内容 | 每次会话都注入，浪费 token 预算 |
| 没有 Pre-Development Checklist | Agent 不知道在哪种情况下读哪个文件 |
| 没有 Quick Reference | Agent 不得不打开子文件查基础信息 |
| 列了文件但 Status 是 TODO | Agent 打开后发现是空模板，浪费一次工具调用 |

---

## 5. 写 guideline 文件 — 规范正文

每个 guideline 文件是**一份具体的编码规范**，被 JSONL 按需注入到 agent 上下文中。

### 5.1 骨架结构

一个好的 guideline 文件通常包含这些节：

```
1. 标题 + 一句话描述（和 index.md 的描述一致）
2. 核心概念 / 架构图 / 类型关系
3. 具体规则，每条配代码示例
4. Forbidden Patterns（不能做的）
5. Required Patterns（必须做的）
6. Checklist（快速自查）
```

不是每个文件都需要所有节。**有代码示例的规则** > 没有代码示例的规则。

### 5.2 写一条规则

每条规则遵循 **Claim → Code → Rationale** 结构：

> 📎 source: `.trellis/spec/core/quality-guidelines.md` 第 26-39 行（F1 规则）

```markdown
### F1: Direct operator multiplication

（Claim: 断言）
FORBIDDEN: Do not multiply operators directly.

（Code: 代码示例，对比 BAD/GOOD）
​```python
# FORBIDDEN
result = T * W

# CORRECT: use NO() for normal-ordered products
result = NO(T, W)
​```

（Rationale: 为什么）
**Why**: In VOA, operator "multiplication" is not well-defined
without specifying ordering. NO(A, B) makes the normal ordering explicit.
```

**为什么要 Rationale？** LLM 不是无脑执行指令。给出理由，agent 在**边缘场景**也能做出正确判断（而不是"我猜这里应该也用 NO"）。

### 5.3 Forbidden vs Required

区分**禁止模式**和**要求模式**，是 spec 最核心的写作技巧。

**Forbidden（禁止模式）** — 告诉 agent "**不要**这样做"：

> 📎 source: `.trellis/spec/core/quality-guidelines.md` 第 40-52 行（F2 规则）

```markdown
### F2: Hand-written expected values in computation tests

​```python
# FORBIDDEN
assert result.pole(4) == "something I think is right"

# CORRECT: derive from OPEdefs.m
# See testing/mathematica-reference.md for the proper workflow
​```

**Why**: VOA computations are subtle. Expected values must come from
Thielemans' OPEdefs.m or published literature.
```

**Required（要求模式）** — 告诉 agent "**必须**这样做"：

> 📎 source: `.trellis/spec/core/quality-guidelines.md` 第 90-96 行（R1 规则）

```markdown
### R1: Register parity before OPE definition

​```python
T = BasisOperator("T", conformal_weight=2)
Bosonic(T)                              # ← REQUIRED
OPE[T, T] = MakeOPE([...])
​```
```

**经验法则**：如果一个错误曾经发生过，就把它加到 Forbidden 中。如果一个步骤容易被遗忘，就加到 Required 中。

### 5.4 Checklist — 快速自查表

每个 guideline 文件末尾放一个 checklist，供 agent（和 check agent）快速扫描：

> 📎 source: `.trellis/spec/core/quality-guidelines.md` 第 122-134 行（Review Checklist）

```markdown
## Review Checklist

Before merging any change:

- [ ] `black --check` passes
- [ ] `ruff check` passes
- [ ] No new forbidden patterns introduced
- [ ] If new computation: Mathematica reference documented
- [ ] If touching operator types: guard preserved
```

Check agent 会逐项验证这些 checklist。写得越具体，check agent 越可靠。

---

## 6. 有效写作的 7 个原则

### 原则 1：面向 LLM，不是面向人

人读文档会跳读、联想上下文。LLM 逐 token 处理，容易丢失远距离上下文。

> 📎 source: `.trellis/spec/core/api-conventions.md` 第 65-70 行（bracket 定义）

```markdown
# 差：假设读者已知上下文
Use the standard bracket convention.

# 好：显式写出
bracket(A, B, n) returns the n-th POLE of OPE(A, B).
WARNING: this is the n-th POLE, NOT the n-th MODE.
```

### 原则 2：代码 > 描述

> 📎 source: `.trellis/spec/core/api-conventions.md` 第 29-42 行（MakeOPE pole ordering）

```markdown
# 差：纯文字
Pole list should go from highest to lowest order.

# 好：代码示例
# Pole ordering: highest-to-lowest
MakeOPE([c/2*One, Zero, 2*T, d(T)])
#        ↑4th pole  ↑3rd   ↑2nd  ↑1st pole
```

Agent 能从代码示例**直接模仿**。纯文字描述需要额外"翻译"步骤，容易出错。

### 原则 3：一条规则 = 一个节

不要把多条规则混在一个段落里。每条规则用独立标题，方便 agent 定位和引用。

> 📎 source: `.trellis/spec/core/quality-guidelines.md`（R1、F5、R2 三条规则拆分示范）

```markdown
# 差：糅在一起
Operators must be registered with Bosonic() or Fermionic() before
defining OPE. Also use Rational instead of float. And don't forget
to clear the registry in tests.

# 好：拆成独立规则
### R1: Register parity before OPE definition
...
### F5: Using Python float for mathematical constants
...
### R2: Test isolation via registry clear
...
```

### 原则 4：用编号标识规则

给 Forbidden/Required 规则编号（F1、F2、R1、R2...），好处：

- Check agent 可以在报告中引用 "违反了 F3"
- Checklist 可以引用 "确认 R1-R4 已遵守"
- 团队沟通时有明确指代

### 原则 5：写**边界条件**，不只写 happy path

> 📎 source: `.trellis/spec/testing/mathematica-reference.md`（三级参考源 + 豁免清单）

```markdown
### Acceptable Reference Sources

| Tier | Source | When to use |
|------|--------|-------------|
| Tier 1 | OPEdefs.m computation | **Preferred** for all OPE/bracket/NO tests |
| Tier 2 | Published literature | When OPEdefs.m cannot compute it |
| Tier 3 | VOA axiom derivation | For structural properties only |

### What Does NOT Need Reference

- Operator declaration (weight, parity)
- Registry setup and teardown
- Error handling (exception tests)
```

告诉 agent **什么时候不适用**这条规则，和告诉它什么时候适用同样重要。

### 原则 6：保持可维护的长度

| 文件类型 | 建议行数 | 原因 |
|---------|---------|------|
| index.md | 30–60 行 | 每次会话都注入，短=省 token |
| guideline 文件 | 50–150 行 | Agent 注意力有限，太长规则被忽略 |
| thinking guide | 40–100 行 | 提供思维框架，不是详细清单 |

如果一个文件超过 150 行，考虑拆分。

### 原则 7：从 bug 中学习

Spec 不是一次性写完的，而是**逐步生长**的：

```
发现 bug → 定位根因 → 写成 Forbidden/Required 规则 → 加入 spec
```

**例子**：第一次发现 `Matrix.rank()` 在大有理数矩阵上会 hang，于是写入 `sympy-patterns.md`：

> 📎 source: `.trellis/spec/core/sympy-patterns.md`（Performance Pitfalls 章节）

```markdown
### Performance Pitfall: Matrix.rank() on Rational matrices

sympy.Matrix.rank() can HANG on large rational matrices.

​```python
# BAD: may hang indefinitely
rank = sp.Matrix(data).rank()

# GOOD: convert to numpy for rank computation
import numpy as np
rank = np.linalg.matrix_rank(np.array(data, dtype=float))
​```
```

这是 spec 最有价值的内容 — **从实战中提取的陷阱和解法**。

---

## 7. 反面模式

### Anti-pattern 1：空模板

> 📎 反面教材: 本项目 `.trellis/spec/backend/` 曾有 5 个空模板文件（`database-guidelines.md` 等），已在本次修复中删除。

```markdown
# Database Guidelines

> To be filled by the team.

## Query Patterns
<!-- TODO -->

## Migrations
<!-- TODO -->
```

**问题**：Agent 打开后发现是空的，浪费 token 和工具调用。**如果还没有内容，不要创建文件。** 等有内容再创建，然后加入 index.md。

### Anti-pattern 2：重定向地狱

> 📎 反面教材: 本项目 `.trellis/spec/backend/index.md` 曾只是一个重定向页（指向 `core/`），而 `init-context` 默认 JSONL 指向它。Agent 只看到"去别处看"，却无法访问实际内容。已在本次修复中将 JSONL 改为直接指向 `core/index.md`。

```markdown
# Backend Guidelines

This project is not a web backend.
Go read core/index.md instead.
```

**问题**：Agent 无法"跟链接"。当 JSONL 指向这个文件时，agent 只看到"去别处看"，但别处的文件并没有注入。

**修复**：让 JSONL 直接指向实际内容文件。在 init-context 的默认配置中使用真正的 spec 路径。

### Anti-pattern 3：百科全书式

```markdown
# SymPy Complete Reference

## Chapter 1: What is SymPy
SymPy is a Python library for symbolic mathematics...
(500 lines of SymPy tutorial)

## Chapter 2: Expression Trees
...
```

**问题**：Agent 不需要学 SymPy。Agent 需要知道**在这个项目中**使用 SymPy 的特定规则和陷阱。

**修复**：只写和项目直接相关的 SymPy pattern，不写通用教程。

### Anti-pattern 4：只有 BAD，没有 GOOD

```markdown
### Don't use float
### Don't import across layers
### Don't assume cache state
```

**问题**：Agent 知道不该做什么，但不知道**应该做什么**。

**修复**：每个 FORBIDDEN 都配一个 CORRECT 示例。

### Anti-pattern 5：规则和项目脱节

> 📎 反面教材: 本项目是数学计算库，Trellis 模板自带的 `database-guidelines.md`、`logging-guidelines.md` 等 web 后端规范完全不适用，已删除。

```markdown
# API Error Response Format

All API responses must follow RFC 7807 Problem Details format...
```

**问题**：这个项目是数学计算库，不是 REST API。保留不适用的模板误导 agent。

**修复**：删除不适用的 spec 文件。只保留和项目相关的规范。

---

## 8. 完整示例拆解

下面以 pyope 的 `testing/` 目录为例，展示一个从 index.md 到 guideline 文件的完整 spec 集合。

> 📎 所有文件均可在 `.trellis/spec/testing/` 和 `.trellis/spec/core/` 目录下直接查看原文件。

### 8.1 目录

```
testing/
├── index.md                  # 40 行：导航 + 速查表
├── test-conventions.md       # ~100 行：测试结构、命名、fixture
├── mathematica-reference.md  # ~80 行：参考数据工作流
└── algebra-fixtures.md       # ~90 行：标准 VOA fixture 复用
```

### 8.2 index.md 拆解

> 📎 source: `.trellis/spec/testing/index.md`（52 行，完整文件）

```markdown
# Testing Guidelines                    ← ① 标题
> How to write and organize tests...    ← ② 一句话描述

## Overview                              ← ③ 关键约束（粗体高亮）
pyope testing has a **unique constraint**: computational results
must be validated against OPEdefs.m or published literature.
Hand-written expected values are **forbidden**.

## Guidelines Index                      ← ④ 文件索引表
| Guide | Description | Status |
|-------|-------------|--------|
| [Test Conventions] | Test structure, fixtures, naming | Filled |
| [Mathematica Reference] | OPEdefs.m reference workflow | Filled |
| [Algebra Fixtures] | Standard VOA fixtures and reuse | Filled |

## Pre-Development Checklist             ← ⑤ 条件化读取顺序
1. **Always**: Test Conventions
2. **If testing OPE/NO/bracket**: Mathematica Reference — mandatory
3. **If defining a VOA algebra**: Algebra Fixtures

## Quick Reference: Test Categories      ← ⑥ 速查表
| Category | Reference source | Example |
|----------|-----------------|---------|
| Computation | OPEdefs.m (mandatory) | test_t_derivative_ope |
| Identity | Mathematical theorem | test_jacobi_identity |
| ...

## The Golden Rule                       ← ⑦ 最重要的一条规则
> Computation tests MUST cite their reference.
```

**观察**：
- 40 行做到了全部 7 个要素
- "The Golden Rule" 让 agent 即使不打开子文件也记住最重要的约束
- Quick Reference 让 agent 在写测试之前快速判断属于哪个分类

### 8.3 guideline 文件拆解

> 📎 source: `.trellis/spec/core/quality-guidelines.md`（134 行，完整文件）

以 `quality-guidelines.md` 为例：

```
第 1-5 行    ：标题 + 描述
第 7-20 行   ：Code Style（工具清单 + 命令）
第 24-85 行  ：Forbidden Patterns（F1-F5，每条含 BAD/GOOD 代码 + Why）
第 89-119 行 ：Required Patterns（R1-R4，每条含代码示例）
第 123-134 行：Review Checklist（8 条 checkbox）
```

**观察**：
- 134 行总长，在推荐范围内
- Forbidden（占最大篇幅）— 因为防错比教对更重要
- 每条规则都有编号（F1-F5, R1-R4）— check agent 可以引用
- Checklist 结尾 — check agent 的执行清单

---

## 附录：Spec 写作 Checklist

写完 spec 后，用这个清单自检：

### index.md

- [ ] 控制在 60 行以内
- [ ] 包含 Overview + **Key constraint**
- [ ] 包含 Guidelines Index 表（含 Status 列）
- [ ] 包含 Pre-Development Checklist（条件化读取顺序）
- [ ] 包含至少一个 Quick Reference 速查表
- [ ] 没有 TODO/空占位符条目

### guideline 文件

- [ ] 控制在 150 行以内
- [ ] 每条规则有独立标题和编号
- [ ] 每条 Forbidden 都有 CORRECT 对照代码
- [ ] 每条规则都有 **Why**（解释原因）
- [ ] 末尾有 Review Checklist
- [ ] 所有代码示例可以直接复制使用

### 整体

- [ ] 没有空模板文件（未填写的文件不应存在）
- [ ] 没有重定向地狱（JSONL 应直接指向实际内容）
- [ ] 每个 spec 文件都和项目相关（删除了不适用的模板）
- [ ] Spec 内容来自实际 bug/经验（不是凭空想象的规则）
