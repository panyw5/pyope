#!/usr/bin/env wolframscript
(* 测试 lev=3 时的字段生成过程 *)

(* 定义必要的函数和变量 *)
opfields={ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};
fields={ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};

(* 定义 gen 函数 *)
gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===ww[0,1,z]):=ww[gr-3/2,1,z];
gen[op_, gr_]/;(gr>=1&&IntegerQ[gr]&&op===jj0[0,1,z]):=jj0[gr-1,1,z];
gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===wwb[0,1,z]):=wwb[gr-3/2,1,z];
gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===tt[0,1,z]):=tt[gr-2,1,z];
gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===gg[0,1,z]):=gg[gr-3/2,1,z];
gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===ggt[0,1,z]):=ggt[gr-3/2,1,z];
gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===ggw[0,1,z]):=ggw[gr-2,1,z];
gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===ggwb[0,1,z]):=ggwb[gr-2,1,z];
gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===ww[0,1,z]):=Nothing;
gen[op_, gr_]/;(!IntegerQ[gr]&&op===jj0[0,1,z]):=Nothing;
gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===wwb[0,1,z]):=Nothing;
gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===tt[0,1,z]):=Nothing;
gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===gg[0,1,z]):=Nothing;
gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===ggt[0,1,z]):=Nothing;
gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===ggw[0,1,z]):=Nothing;
gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===ggwb[0,1,z]):=Nothing;

OperatorList[oplis_,grlis_]:=Table[gen[op, gr],{gr,grlis},{op,oplis}]

OperatorMult[a_,b_,right__]:=NO[a,b,right];
OperatorMult[a_,b_]:=NO[a,b];
OperatorMult[a_]:=NO[a];

GenOp[partlis_]:=Flatten@Table[Outer[OperatorMult,Sequence@@OperatorList[opfields,partition]],{partition,partlis}]

(* 设置 lev = 3 *)
lev = 3;

(* 生成分区 *)
Print["=== Integer Partitions for lev=3 ==="];
part = 1/2 IntegerPartitions[2 lev];
Print["part = ", part];
Print[""];

(* 生成所有算符组合 *)
Print["=== Generating operators ==="];
generated = Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
Print["Number of generated operators: ", Length[generated]];
Print[""];

(* 删除重复项 *)
Print["=== After DeleteDuplicates ==="];
deduplicated = DeleteDuplicates[generated];
Print["Number after deduplication: ", Length[deduplicated]];
Print[""];

(* 检查是否包含 NO[gg[0,1,z], gg[0,1,z]] *)
Print["=== Checking for NO[gg[0,1,z], gg[0,1,z]] ==="];
hasGGGG = MemberQ[deduplicated, NO[gg[0,1,z], gg[0,1,z]]];
Print["Contains NO[gg[0,1,z], gg[0,1,z]]: ", hasGGGG];
Print[""];

(* 查找所有包含 gg 的项 *)
Print["=== All terms containing gg ==="];
ggTerms = Select[deduplicated, !FreeQ[#, gg]&];
Print["Number of terms with gg: ", Length[ggTerms]];
Do[Print[ggTerms[[i]]], {i, Min[10, Length[ggTerms]]}];
Print[""];

(* 最终转换 *)
Print["=== Final generation (first 20 terms) ==="];
finalgen = Table[If[Head[deduplicated[[i]]]===NO, ReplaceAll[deduplicated[[i]], NO->Nl], Nl[deduplicated[[i]]]], {i, Length[deduplicated]}];
Print["Number in finalgen: ", Length[finalgen]];

(* 检查 finalgen 中是否有 Nl[gg[0,1,z], gg[0,1,z]] *)
hasNlGGGG = MemberQ[finalgen, Nl[gg[0,1,z], gg[0,1,z]]];
Print["Contains Nl[gg[0,1,z], gg[0,1,z]]: ", hasNlGGGG];
Print[""];

(* 显示前20个 finalgen 项 *)
Do[Print[finalgen[[i]]], {i, Min[20, Length[finalgen]]}];
