#!/usr/bin/env wolframscript
(* 正确使用 % 符号重现原代码逻辑 *)

opfields={ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};
fields={ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};

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

lev = 3;
part = 1/2 IntegerPartitions[2 lev];

Print["=== 原代码的三行 ==="];
Print["Line 1: generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;"];
generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
Print["  generated length: ", Length[generated]];
Print["  Contains NO[gg[0,1,z], gg[0,1,z]]: ", MemberQ[generated, NO[gg[0,1,z], gg[0,1,z]]]];
Print[""];

Print["Line 2: DeleteDuplicates[%];"];
(* 这里 % 指的是上一行的输出，即 generated *)
deduplicated = DeleteDuplicates[generated];
Print["  deduplicated length: ", Length[deduplicated]];
Print["  Contains NO[gg[0,1,z], gg[0,1,z]]: ", MemberQ[deduplicated, NO[gg[0,1,z], gg[0,1,z]]]];
Print[""];

Print["Line 3: finalgen=Table[If[Head[%[[i]]]===NO,ReplaceAll[%[[i]],NO->Nl],Nl[%[[i]]]],{i,Length[%]}]"];
(* 这里的 % 指的是上一行的输出，即 deduplicated *)
finalgen=Table[If[Head[deduplicated[[i]]]===NO,ReplaceAll[deduplicated[[i]],NO->Nl],Nl[deduplicated[[i]]]],{i,Length[deduplicated]}];
Print["  finalgen length: ", Length[finalgen]];
Print["  Contains Nl[gg[0,1,z], gg[0,1,z]]: ", MemberQ[finalgen, Nl[gg[0,1,z], gg[0,1,z]]]];
Print[""];

Print["=== finalgen 中包含 gg 的所有项 ==="];
ggTerms = Select[finalgen, !FreeQ[#, gg]&];
Print["Number of terms with gg: ", Length[ggTerms]];
Do[Print[i, ": ", ggTerms[[i]]], {i, Length[ggTerms]}];
