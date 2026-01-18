#!/usr/bin/env wolframscript
(* 完整测试：包括量子数分组步骤 *)

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

(* 设置 lev = 3 *)
lev = 3;
part = 1/2 IntegerPartitions[2 lev];
generated = Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
DeleteDuplicates[generated];
finalgen = Table[If[Head[%[[i]]]===NO, ReplaceAll[%[[i]], NO->Nl], Nl[%[[i]]]], {i, Length[%]}];

Print["=== finalgen 中包含 gg 的项 ==="];
ggTermsInFinalgen = Select[finalgen, !FreeQ[#, gg]&];
Print["Number of gg terms in finalgen: ", Length[ggTermsInFinalgen]];
Do[Print[i, ": ", ggTermsInFinalgen[[i]]], {i, Length[ggTermsInFinalgen]}];
Print[""];

(* 计算量子数 *)
Print["=== 计算 Nl[gg[0,1,z], gg[0,1,z]] 的量子数 ==="];
testTerm = Nl[gg[0,1,z], gg[0,1,z]];
testTermConverted = testTerm/.Nl->List/.{ww->dww,jj0->djj0,wwb->dwwb,tt->dtt,gg->dgg,ggt->dggt,ggw->dggw,ggwb->dggwb};
testTermProduct = Times@@testTermConverted;
Print["Converted term: ", testTermProduct];
Print["m quantum number: ", mtotal[testTermProduct]];
Print["r quantum number: ", rtotal[testTermProduct]];
Print["f quantum number (fermion parity): ", ftotal[testTermProduct]];
Print[""];

(* 检查费米子统计 *)
Print["=== 费米子统计检查 ==="];
Print["gg 是费米子，fnumber[gg] = ", fnumber[gg]];
Print["Nl[gg, gg] 的费米子数 = ", ftotal[testTermProduct]];
Print["费米子数为偶数 (2)，说明这是两个费米子的乘积"];
Print[""];

(* 检查 Nl 的性质 *)
Print["=== Nl 的反对称性 ==="];
Print["对于费米子，Nl[gg, gg] 在数学上应该等于多少？"];
Print["如果 Nl 代表正规序 (normal ordered product)，"];
Print["对于费米子算符 ψ，有 :ψψ: = 0 (费米子的反对称性)"];
Print["因此 Nl[gg[0,1,z], gg[0,1,z]] 在数学上应该为 0"];
