#!/usr/bin/env wolframscript
(* 详细追踪 finalgen 的处理过程 *)

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
part = 1/2 IntegerPartitions[2 lev];

Print["=== Step 1: GenOp[part] ==="];
step1 = GenOp[part];
Print["Length: ", Length[step1]];
ggInStep1 = Select[step1, !FreeQ[#, gg]&];
Print["Terms with gg: ", Length[ggInStep1]];
Print[""];

Print["=== Step 2: Select[..., !NumberQ[#]&] ==="];
step2 = Select[step1, !NumberQ[#]&];
Print["Length: ", Length[step2]];
ggInStep2 = Select[step2, !FreeQ[#, gg]&];
Print["Terms with gg: ", Length[ggInStep2]];
Print[""];

Print["=== Step 3: /.-1*x_:>x ==="];
step3 = step2/.-1*x_:>x;
Print["Length: ", Length[step3]];
ggInStep3 = Select[step3, !FreeQ[#, gg]&];
Print["Terms with gg: ", Length[ggInStep3]];
Print[""];

Print["=== Step 4: DeleteDuplicates ==="];
step4 = DeleteDuplicates[step3];
Print["Length: ", Length[step4]];
ggInStep4 = Select[step4, !FreeQ[#, gg]&];
Print["Terms with gg: ", Length[ggInStep4]];
Do[Print["  ", ggInStep4[[i]]], {i, Length[ggInStep4]}];
Print[""];

Print["=== Step 5: NO -> Nl conversion ==="];
step5 = Table[If[Head[step4[[i]]]===NO, ReplaceAll[step4[[i]], NO->Nl], Nl[step4[[i]]]], {i, Length[step4]}];
Print["Length: ", Length[step5]];
ggInStep5 = Select[step5, !FreeQ[#, gg]&];
Print["Terms with gg: ", Length[ggInStep5]];
Do[Print["  ", ggInStep5[[i]]], {i, Length[ggInStep5]}];
Print[""];

Print["=== 检查 Nl[gg[0,1,z], gg[0,1,z]] 是否在 step5 中 ==="];
Print["MemberQ result: ", MemberQ[step5, Nl[gg[0,1,z], gg[0,1,z]]]];
