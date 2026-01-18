#!/usr/bin/env wolframscript
(* 完整追踪到字符计算 *)

(* 加载 OPE 包 - 如果可用 *)
(* <<ope.math.m *)
(* <<OPEDefs` *)

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

dww[x_,y_,z_]:=ww;
djj0[x_,y_,z_]:=jj0;
dwwb[x_,y_,z_]:=wwb;
dtt[x_,y_,z_]:=tt;
dgg[x_,y_,z_]:=gg;
dggt[x_,y_,z_]:=ggt;
dggw[x_,y_,z_]:=ggw;
dggwb[x_,y_,z_]:=ggwb;

mnumber[x_]:=Association[{ww->3/2,jj0->0,wwb->-(3/2),tt->0,gg->-(1/2),ggt->1/2,ggw->1,ggwb->-1}]@x
rnumber[x_]:=Association[{ww->0,jj0->0,wwb->0,tt->0,gg->1/2,ggt->-(1/2),ggw->1/2,ggwb->-(1/2)}]@x
fnumber[x_]:=Association[{ww->0,jj0->0,wwb->0,tt->0,gg->1,ggt->1,ggw->1,ggwb->1}]@x

mtotal[x_]:=Module[{y=x,m=0,i=1},For[i=1,i<Length[fields]+1,i++,m=m+Exponent[y,fields[[i]]]mnumber[fields[[i]]]];m]
rtotal[x_]:=Module[{y=x,r=0,i=1},For[i=1,i<Length[fields]+1,i++,r=r+Exponent[y,fields[[i]]]rnumber[fields[[i]]]];r]
ftotal[x_]:=Module[{y=x,r=0,i=1},For[i=1,i<Length[fields]+1,i++,r=r+Exponent[y,fields[[i]]]fnumber[fields[[i]]]];r]

lev = 3;
part = 1/2 IntegerPartitions[2 lev];
generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
deduplicated = DeleteDuplicates[generated];
finalgen=Table[If[Head[deduplicated[[i]]]===NO,ReplaceAll[deduplicated[[i]],NO->Nl],Nl[deduplicated[[i]]]],{i,Length[deduplicated]}];

Print["=== Step 1: finalgen 生成 ==="];
Print["finalgen length: ", Length[finalgen]];
Print["Contains Nl[gg[0,1,z], gg[0,1,z]]: ", MemberQ[finalgen, Nl[gg[0,1,z], gg[0,1,z]]]];
Print["Position: ", Position[finalgen, Nl[gg[0,1,z], gg[0,1,z]]]];
Print[""];

Print["=== Step 2: 转换为 fulllist ==="];
fulllist = finalgen/.Nl->List;
fulllist = Times@@#&/@fulllist;
fulllist = fulllist/.{ww->dww,jj0->djj0,wwb->dwwb,tt->dtt,gg->dgg,ggt->dggt,ggw->dggw,ggwb->dggwb};
Print["fulllist length: ", Length[fulllist]];
Print["Contains gg^2: ", MemberQ[fulllist, gg^2]];
ggPos = Position[fulllist, gg^2];
If[Length[ggPos] > 0,
  Print["Position of gg^2: ", ggPos];
  Print["Corresponding finalgen term: ", finalgen[[ggPos[[1,1]]]]];
];
Print[""];

Print["=== Step 3: 计算量子数 ==="];
Print["Quantum numbers of gg^2:"];
Print["  m = ", mtotal[gg^2]];
Print["  r = ", rtotal[gg^2]];
Print["  f = ", ftotal[gg^2]];
Print[""];

Print["=== 关键观察 ==="];
Print["gg 是费米子 (fnumber[gg] = ", fnumber[gg], ")"];
Print["两个相同费米子的正规序乘积 Nl[gg, gg] 在物理上应该为 0"];
Print["这是因为费米子满足反对易关系：{ψ, ψ} = 0"];
Print["因此 :ψψ: = (1/2){ψ, ψ} = 0"];
Print[""];
Print["代码在 finalgen 生成阶段**保留**了这一项，"];
Print["但在后续的 OPE 计算或 null state 分析中会发现它为零。"];
