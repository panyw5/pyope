#!/usr/bin/env wolframscript
(* Direct copy from null_states.wls with lev=4 *)

(* Load Packages *)
SetDirectory["/Users/lelouch/pyope/OPEdefs"];
<<ope.math.m
<<OPEDefs`

(* Settings *)
fcal=Checkbox[False];
lev=4;

(* Definitions *)
opfields={ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};
fields={ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};

Bosonic[β1,γ1,ww[i_,j_,z],jj0[i_,j_,z],wwb[i_,j_,z],tt[i_,j_,z]];
Fermionic[b1,c1,gg[i_,j_,z],ggt[i_,j_,z],ggw[i_,j_,z],ggwb[i_,j_,z]];

OPE[β1,γ1]=MakeOPE[-One (z-w)^-1+Ord[z,w,0]];
OPE[b1,c1]=MakeOPE[One (z-w)^-1+Ord[z,w,0]];

Derivativer[x_]:=x';

wwww[i_,j_,z_]:=Nest[Derivativer,w,i];
j0j0[i_,j_,z_]:=Nest[Derivativer,j0,i];
wbwb[i_,j_,z_]:=Nest[Derivativer,wb,i];
tttt[i_,j_,z_]:=Nest[Derivativer,t,i];
gggg[i_,j_,z_]:=Nest[Derivativer,g,i];
gtgt[i_,j_,z_]:=Nest[Derivativer,gt,i];
gwgw[i_,j_,z_]:=Nest[Derivativer,gw,i];
gwbgwb[i_,j_,z_]:=Nest[Derivativer,gwb,i];

w=β1;
j0=2 NO[b1,c1]+3 NO[β1,γ1];
wb=NO[β1,NO[β1,NO[γ1,NO[γ1,γ1]]]]+2 NO[β1,NO[γ1,NO[γ1,NO[b1,c1]]]]-4 NO[β1,NO[Derivative[1][γ1],γ1]]-4/3 NO[γ1,NO[b1,Derivative[1][c1]]]+2/3 NO[γ1,NO[Derivative[1][b1],c1]]+2/3 NO[Derivative[1][β1],NO[γ1,γ1]]-8/3 NO[Derivative[1][γ1],NO[b1,c1]]+(10 γ1'')/9;
t=-2 NO[b1,Derivative[1][c1]]-3/2 NO[β1,Derivative[1][γ1]]-NO[Derivative[1][b1],c1]-1/2 NO[Derivative[1][β1],γ1];
g=NO[γ1,b1];
gt=2NO[β1',c1]+3NO[β1,c1'];
gw=b1;
gwb=8/3 NO[b1,NO[c1'',c1]]+3 NO[β1,NO[β1,NO[γ1,NO[γ1,Derivative[1][c1]]]]]-4 NO[β1,NO[γ1,NO[b1,NO[Derivative[1][c1],c1]]]]-4 NO[β1,NO[γ1,c1'']]-4 NO[β1,NO[Derivative[1][γ1],Derivative[1][c1]]]-2/3 NO[Derivative[1][b1],NO[Derivative[1][c1],c1]]+2 NO[Derivative[1][β1],NO[β1,NO[γ1,NO[γ1,c1]]]]-8/3 NO[Derivative[1][β1],NO[Derivative[1][γ1],c1]]+2/3 NO[β1'',NO[γ1,c1]]+(10 c1''')/9;

(* Functions *)
NO[x_]:=x;

pwrexpand[x_]:=Module[
{y=x,list={},i,j,p=Variables[x]},
For[i=1,i<Length[p]+1,i++,For[j=1,j<(Exponent[y,p[[i]]]+1),j++,list=Append[list,Head[p[[i]]]]]];
list
]

parorder[x_]:=Block[
{y=x,l1,r},
l1=Table[NO@@y[[i]],{i,Length[y]}];
r=Total[l1];
r
]

parordermod[x_]:=Block[
{y=x,r},
r=NO@@y;
r
]

(* Field Generation *)
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

(* Level *)
part=1/2IntegerPartitions[2 lev]

generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
DeleteDuplicates[%];
finalgen=Table[If[Head[%[[i]]]===NO,ReplaceAll[%[[i]],NO->Nl],Nl[%[[i]]]],{i,Length[%]}]

Print["Generated fields at level ", lev, ":"];
Print[finalgen];
Print["\nTotal count: ", Length[finalgen]];
