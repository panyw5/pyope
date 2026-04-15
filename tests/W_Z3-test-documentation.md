# `W_Z3` 测试文档

本文档说明 `W_Z3` 相关测试在验证什么、覆盖了哪些函数/类/数据构造，以及这些测试所要求的 expected behavior。

相关文件：

- `tests/algebras/test_w_z3_algebra.py`
- `tests/algebras/test_w_z3_algebra_ref.py`
- `tests/utils/w_z3_algebra.py`
- `tests/W_Z3-algebra.md`

## 测试目标

这组测试围绕 $\mathcal{W}_{\mathbb{Z}_3}$ 代数做三类验证：

1. 从 `bc beta-gamma` free-field realization 出发，检查 `pyope` 是否能恢复文档给出的 OPE。
2. 只保留抽象生成元和 OPE 表时，检查 Jacobi 残差是否与文档整理出的 null expressions 一致。
3. 检查文档列出的选定 weight-8 null states 在 free-field realization 下是否全部化为零。

同时，测试还保留了一条 Wolfram 参考链路：若环境中安装了 `wolframscript`，就额外用 `OPEdefs.m` / `OPEdefs.wls` 做交叉验证；若没有安装，则相关测试必须自动跳过。

## 被测试的对象

### 1. 共享数据与解析工具

定义于 `tests/utils/w_z3_algebra.py`：

- `W_Z3_OPE_GROUND_TRUTH`
  - 文档中非平凡 OPE 的结构化转录。
  - expected behavior：每一项都能被解析成精确的 `pyope`/SymPy 表达式，并作为后续断言的 ground truth。

- `W_Z3_NULL_GROUND_TRUTH`
  - 记录被纳入测试的 weight-8 null-state 编号与特殊关系。
  - expected behavior：测试只检查这份白名单中的 relation，数量必须稳定，不能漏项。

- `W_Z3_JACOBI_GROUND_TRUTH`
  - 记录抽象 Jacobi 检查使用的生成元顺序与可直接比较的三元组。
  - expected behavior：支持直接检查的三元组应当给出零 Jacobi 矩阵；其余三元组产生的非零残差应落在文档列出的 null expressions 中。

- `_WLParser` / `parse_wl_expr`
  - 把 Wolfram 风格字符串解析为 `pyope` 表达式。
  - expected behavior：能正确识别 `NO[...]`、`Derivative[n][...]`、`One`、有理数和生成元名字，供 OPE 与 null-state 资料复用。

### 2. 生成元构造函数

- `make_z3_free_field_data(prefix="...")`
  - 构造 `b, c, beta, gamma` free fields，以及 `T, J, W, Wbar, G, Gbar, GW, GbarWbar` 的具体 realization。
  - expected behavior：
    - 正确声明权重与费米/玻色统计；
    - 正确注册 `b(z)c(w)` 与 `beta(z)gamma(w)` 的基础 OPE；
    - `GW` 与 `GbarWbar` 必须由一级极点 `bracket(..., 1)` 得到，而不是手写常量。

- `make_z3_abstract_data(prefix="...")`
  - 构造抽象强生成元，并把 `W_Z3_OPE_GROUND_TRUTH` 装入 OPE 注册表。
  - expected behavior：
    - 抽象对象只保留生成元和 OPE 关系，不依赖 free-field realization；
    - 后续 Jacobi 检查应当完全基于这张 OPE 表运行。

- `make_z3_realization_target_data(prefix="...")`
  - 构造一套单独的“目标生成元 + realization map”。
  - expected behavior：抽象表达式在 realization 时不会直接复用原生成元对象，从而避免对象混淆。

### 3. 文档转断言的辅助函数

- `expected_ope_expr_map(ops)`
  - 将文档 OPE 表解析为精确表达式。
  - expected behavior：返回的 pole map 能与 `OPE(...)` 结果逐极点比较。

- `load_selected_null_relation_sources()`
  - 从 `tests/W_Z3-algebra.md` 提取选定的 weight-8 null-state 源文本。
  - expected behavior：必须能找到所有 `basis_ids` 与 `Particular T4 Relation`；缺块时应抛错，避免静默漏测。

- `build_null_relations(ops)`
  - 解析并规范化 null-state 关系。
  - expected behavior：
    - 从 Markdown 文档读取关系；
    - 转成 `pyope` 表达式；
    - 经过 `simplify_with_wolfram` 规范化后，作为“应当等于零”的关系集合返回。

### 4. Jacobi 相关辅助函数

- `collect_python_jacobi_residuals(ops=None, include_zero=False)`
  - 收集抽象生成元上的 Jacobi 残差。
  - expected behavior：
    - 默认只收集非零残差；
    - 这些残差应该与文档整理出的 Jacobi null expressions 对齐。

- `expected_jacobi_null_expressions(ops)`
  - 解析完整的 Jacobi null-expression 列表。
  - expected behavior：生成出的表达式在 free-field realization 下应全部化零。

- `expected_jacobi_null_expressions_canonical(ops)`
  - 解析去重后的 canonical 列表。
  - expected behavior：适合与 `collect_python_jacobi_residuals` 的结果做集合比较。

- `canonicalize_up_to_sign(expr)`
  - 对表达式做“忽略整体符号”的规范化。
  - expected behavior：若两个残差只差一个负号，则应当被视为同一个 relation。

- `realize_exprs_with_free_fields(exprs, abstract_ops=None)`
  - 将抽象表达式递归替换到 free-field realization。
  - expected behavior：抽象 Jacobi/null-expression 经 realization 后，交给 Wolfram 简化时应能变成零。

## 测试文件与 expected behavior

### `tests/algebras/test_w_z3_algebra.py`

这是 `pyope` 主链路测试，直接检查 Python 侧的计算与零化行为。

#### `TestWZ3AlgebraComputations.test_pyope_recovers_documented_ope_from_free_field_ground_truth`

- 被测对象：`make_z3_free_field_data`、`expected_ope_expr_map`、`pyope.OPE`
- expected behavior：
  - 对文档列出的每个非平凡 OPE，`OPE(left, right)` 的每一个极点都应与 ground truth 完全一致；
  - 若某一极点在文档中没有给出，则视为零；
  - 即使 `result.max_pole` 和 expected pole 数不同，也不能漏检高阶零极点。

#### `TestWZ3AlgebraComputations.test_pyope_checks_abstract_jacobi_against_ground_truth`

- 被测对象：`make_z3_abstract_data`、`collect_python_jacobi_residuals`、`expected_jacobi_null_expressions_canonical`、`canonicalize_up_to_sign`
- expected behavior：
  - 只用抽象生成元与 OPE 表时，`pyope` 收集到的非零 Jacobi 残差集合，应与文档去重后的 canonical null-expression 集合完全相同；
  - 比较时允许相差整体负号。

#### `TestWZ3AlgebraComputations.test_pyope_realizes_abstract_jacobi_residuals_to_zero`

- 被测对象：`collect_python_jacobi_residuals`、`realize_exprs_with_free_fields`、`simplify_with_wolfram`
- expected behavior：
  - 抽象层面出现的 Jacobi 残差，不应是“真实矛盾”；
  - 它们在 free-field realization 下都应该被解释为零；
  - Wolfram 简化后的结果必须全部为 `0` 或 `Zero`。

#### `TestWZ3AlgebraComputations.test_pyope_checks_selected_weight8_null_states_against_ground_truth`

- 被测对象：`build_null_relations`、`W_Z3_NULL_GROUND_TRUTH`、`simplify_with_wolfram`
- expected behavior：
  - 选定的 `19` 个 basis relation 加上 `1` 个特殊关系，数量必须固定；
  - 每条关系在 free-field realization 下经 Wolfram 简化后都应为零；
  - 若有任何一条不为零，就说明文档转录、表达式构造或实现存在错误。

### `tests/algebras/test_w_z3_algebra_ref.py`

这是 Wolfram 参考链路测试，重点是“外部参考是否与 Python 测试口径一致”。

整个测试类带有：

- `@pytest.mark.mathematica_ref`
- `@pytest.mark.skipif(shutil.which("wolframscript") is None, ...)`

expected behavior：如果环境里没有 `wolframscript`，这些测试必须被跳过，而不是失败。

#### `TestWZ3WolframReference.test_free_field_ope_recovery_matches_ground_truth`

- 被测对象：free-field OPE 恢复链路
- 当前 expected behavior：
  - 设计目标是用 Wolfram 检查 free-field OPE 差值全部为零；
  - 但当前测试体内显式 `pytest.skip(...)`，因为“大型 free-field OPE diff 的 Wolfram round-trip 还不稳定”；
  - 因而此测试现在的预期结果是“主动跳过”，不是通过或失败。

#### `TestWZ3WolframReference.test_abstract_jacobi_null_expressions_vanish_after_realization`

- 被测对象：`expected_jacobi_null_expressions`、`realize_exprs_with_free_fields`、`simplify_with_wolfram`
- expected behavior：
  - 文档中完整列出的 Jacobi null expressions，在 free-field realization 下经 Wolfram 简化都应为零；
  - 这是对“文档中的 Jacobi 关系本身自洽”的参考验证。

#### `TestWZ3WolframReference.test_free_field_null_relations_vanish_in_wolfram`

- 被测对象：`build_null_relations`、`W_Z3_NULL_GROUND_TRUTH`、`simplify_with_wolfram`
- expected behavior：
  - 文档选取的 null relations 经 Wolfram 计算应全部化零；
  - 关系条数必须等于 `len(basis_ids) + 1`；
  - 这是对 Python 侧同名测试的外部参考版复核。

## 这组测试实际在保证什么

从整体上看，这组测试不是只在检查某个单独函数是否返回了某个值，而是在保证下面几件事同时成立：

- 文档 `tests/W_Z3-algebra.md` 中记录的 OPE、Jacobi 关系和 null states 被正确转录进测试代码；
- `pyope` 在具体 free-field realization 上算出来的 OPE 与文档一致；
- `pyope` 在抽象 OPE 表上得到的 Jacobi 残差，与文档整理出的 null-expression 结构一致；
- 这些抽象残差和 null states 在具体 realization 下确实全部化零；
- Wolfram 参考环境存在时，Python 侧结论能被外部实现再次核对；
- Wolfram 环境不存在时，测试流程仍然稳定，不会因为缺少外部工具而误报失败。

## 适合后续补充到文档里的内容

如果后续还要继续扩展这份文档，最自然的补充方向有两类：

1. 给每个 Jacobi/null relation 补上更细的来源位置，例如对应 `tests/W_Z3-algebra.md` 中的节名或 basis 编号。
2. 单独增加“当前已知限制”小节，记录为什么 `test_free_field_ope_recovery_matches_ground_truth` 现在需要主动跳过。
