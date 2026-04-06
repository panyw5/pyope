---
name: voa
description: 进行形顶点算符代数 (vertex operator algebra) 的解析、符号化的计算。当需要计算顶点算符代数的 OPE、normal ordered product, Zhu's C2 代数，Zhu's 代数，associate variety 等对象时使用，当需要为其它顶点算符代数程序包书写测试 (test) 时使用。
---


# `OPEdefs` manual

阅读 `OPEdefs-manual.md` (相对于本 SKILL.md 的路径是 [OPEdefs-manual.md](manual/OPEdefs-manual.md))

# vertex operator algebra preliminaries

阅读 `VOA-manual.md` (相对于本 SKILL.md 的路径是 [VOA-manual.md](manual/VOA-manual.md))


# `wls` Coding convention
1. **CRITICAL**: 多行表达式必须用括号 `()` 括起来作为一个整体
   否则 Wolfram Script 可能会误解表达式的结构，将第二行以及之后的行当成新的表达式，导致语法错误或计算错误
2. 函数要用中括号 `[arg1, arg2, ...]` 包裹 arguments
3. **注意**: 变量名、argument 名字不能加下划线
4. 下划线代表 `pattern`，小心使用

# `pyope` local operator API

如果任务涉及 `pyope` 里的局域算符枚举、标准化、坐标展开或 `C2` 相关计算，优先按当前接口约定工作：

1. `LocalOperatorBasis` 是当前主接口，负责局域算符基底的枚举、canonicalize、coordinates 等操作。
2. `LocalOperatorCanonicalizer` 现在只是 `LocalOperatorBasis` 的 **backward-compatible alias**。
3. 新代码、新测试、新文档都应优先写成 `LocalOperatorBasis(...)`，不要再把 `LocalOperatorCanonicalizer` 当成独立主类介绍。
4. 如果需要提兼容性，可以说明：旧代码里出现 `LocalOperatorCanonicalizer` 时，可将其理解为 `LocalOperatorBasis` 的旧名/兼容别名。
