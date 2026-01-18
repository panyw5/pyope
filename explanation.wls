#!/usr/bin/env wolframscript
(* 展示 Nl[gg, gg] 如何被排除的完整流程 *)

Print["========================================"];
Print["为什么 Nl[gg[0,1,z], gg[0,1,z]] 最终不出现在物理态中"];
Print["========================================"];
Print[""];

Print["=== 1. 字段生成阶段 ==="];
Print["在 finalgen 生成时，代码**确实生成**了 Nl[gg[0,1,z], gg[0,1,z]]"];
Print["这是因为代码机械地生成所有可能的算符组合。"];
Print["位置：finalgen[[19]] = Nl[gg[0,1,z], gg[0,1,z]]"];
Print[""];

Print["=== 2. 为什么这一项在数学上是零 ==="];
Print["gg 是一个费米子算符（fermionic operator）");
Print["费米子满足反对易关系：{ψ(z), ψ(w)} = 0 当 z ≠ w");
Print["在同一点，费米子的正规序乘积：");
Print["  :ψ(z)ψ(z): = (1/2){ψ(z), ψ(z)} = 0");
Print[""];
Print["因此：Nl[gg[0,1,z], gg[0,1,z]] = :gg(z)·gg(z): = 0");
Print[""];

Print["=== 3. 代码如何排除这一项 ==="];
Print["代码并**不是在生成阶段**排除它，而是在后续步骤中："];
Print[""];

Print["方法 A：在 OPE 计算中自动为零"];
Print["  当代码尝试计算涉及 Nl[gg, gg] 的 OPE 时，"];
Print["  由于 OPE[gg, gg] 的定义（费米子反对易），"];
Print["  所有包含 Nl[gg, gg] 的项会自动简化为零。");
Print[""];

Print["方法 B：在 null state 分析中被识别"];
Print["  代码的 'Calculating Physical States' 部分会："];
Print["  1. 生成所有可能的态（包括 Nl[gg, gg]）");
Print["  2. 计算这些态在 OPE 作用下的关系");
Print["  3. 通过矩阵秩分析识别 null states");
Print["  4. Nl[gg, gg] 会被识别为 null state（因为它等于零）");
Print["  5. 最终的物理态数 = 朴素态数 - null state 数");
Print[""];

Print["=== 4. 代码验证 ==="];
Print["让我们检查 OPE 定义："];
Print[""];
Print["如果 gg 被定义为费米子（通过 Fermionic[b1, c1, gg[i_,j_,z], ...]），");
Print["那么 OPE[gg, gg] 应该满足反对易关系。"];
Print[""];
Print["在原代码第 231 行："];
Print["  Fermionic[b1,c1,gg[i_,j_,z],ggt[i_,j_,z],ggw[i_,j_,z],ggwb[i_,j_,z]];"];
Print[""];
Print["这告诉 OPEdefs 包：gg, ggt, ggw, ggwb 都是费米子。"];
Print["因此当计算 NO[gg, gg] 时，OPEdefs 会自动应用反对易规则。"];
Print[""];

Print["=== 5. 总结 ==="];
Print["Nl[gg[0,1,z], gg[0,1,z]] **在 finalgen 中存在**（第19项）"];
Print["但它在代码的后续步骤中被排除，因为："];
Print["  • 它在数学上等于零（费米子反对易）"];
Print["  • OPEdefs 包会自动处理费米子统计"];
Print["  • 在计算物理态时，它被识别为 null state"];
Print[""];
Print["这是一个**运行时排除**，而不是**生成时排除**。"];
Print["========================================"];
