#!/usr/bin/env wolframscript
(* 追踪量子数分组过程 *)

(* 加载所有必要的定义 *)
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

(* 量子数函数 *)
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

(* 生成 finalgen *)
lev = 3;
part = 1/2 IntegerPartitions[2 lev];
generated = Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
DeleteDuplicates[generated];
finalgen = Table[If[Head[%[[i]]]===NO, ReplaceAll[%[[i]], NO->Nl], Nl[%[[i]]]], {i, Length[%]}];

Print["=== finalgen 生成完成 ==="];
Print["Length: ", Length[finalgen]];
Print["Contains Nl[gg[0,1,z], gg[0,1,z]]: ", MemberQ[finalgen, Nl[gg[0,1,z], gg[0,1,z]]]];
Print[""];

(* 按照原代码进行量子数分组 *)
Print["=== 量子数转换 ==="];
fulllist = finalgen/.Nl->List;
fulllist = Times@@#&/@fulllist;
fulllist = fulllist/.{ww->dww,jj0->djj0,wwb->dwwb,tt->dtt,gg->dgg,ggt->dggt,ggw->dggw,ggwb->dggwb};
Print["fulllist length: ", Length[fulllist]];
Print[""];

(* 检查 gg^2 *)
Print["=== 检查 gg^2 在 fulllist 中 ==="];
Print["MemberQ[fulllist, gg^2]: ", MemberQ[fulllist, gg^2]];
ggTermsInFulllist = Select[fulllist, !FreeQ[#, gg]&];
Print["Terms with gg in fulllist: ", Length[ggTermsInFulllist]];
Do[Print["  ", ggTermsInFulllist[[i]]], {i, Length[ggTermsInFulllist]}];
Print[""];

(* 计算 m 量子数 *)
Print["=== 计算 m 量子数 ==="];
mtset = mtotal[#]&/@fulllist;
Print["mtset: ", Union[mtset]];
Print["m quantum number of gg^2: ", mtotal[gg^2]];
Print[""];

(* 按 m 分组 *)
Print["=== 按 m 量子数分组 ==="];
mgroupstates = Table[Append[{},finalgen[[#]]]&/@Position[mtset,i]//Flatten,{i,Min[mtset],Max[mtset],1/2}]/.{}->Nothing;
Print["Number of m groups: ", Length[mgroupstates]];
Print[""];

(* 检查每个 m 组中是否有 gg 项 *)
Print["=== 检查每个 m 组 ==="];
Do[
  ggInGroup = Select[mgroupstates[[i]], !FreeQ[#, gg]&];
  If[Length[ggInGroup] > 0,
    Print["m group ", i, " (m=", Min[mtset] + (i-1)/2, "): ", Length[ggInGroup], " terms with gg"];
    Do[Print["    ", ggInGroup[[j]]], {j, Min[3, Length[ggInGroup]]}];
  ];
, {i, Length[mgroupstates]}];
