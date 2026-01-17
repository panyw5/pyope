#!/usr/bin/env wolframscript
(* Final comparison: original code vs expected results *)

(* Load Packages *)
SetDirectory["/Users/lelouch/pyope/OPEdefs"];
<<ope.math.m
<<OPEDefs`

(* Expected result from user *)
expected = {
  Nl[jj0[3, 1, z]], Nl[tt[2, 1, z]], Nl[ggw[2, 1, z]],
  Nl[ggwb[2, 1, z]], Nl[jj0[0, 1, z], jj0[2, 1, z]],
  Nl[jj0[0, 1, z], tt[1, 1, z]], Nl[jj0[0, 1, z], ggw[1, 1, z]],
  Nl[jj0[0, 1, z], ggwb[1, 1, z]], Nl[ww[0, 1, z], ww[1, 1, z]],
  Nl[ww[1, 1, z], wwb[0, 1, z]], Nl[ww[1, 1, z], gg[0, 1, z]],
  Nl[ww[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], wwb[1, 1, z]],
  Nl[wwb[0, 1, z], wwb[1, 1, z]], Nl[wwb[1, 1, z], gg[0, 1, z]],
  Nl[wwb[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], gg[1, 1, z]],
  Nl[wwb[0, 1, z], gg[1, 1, z]], Nl[gg[0, 1, z], gg[1, 1, z]],
  Nl[gg[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], ggt[1, 1, z]],
  Nl[wwb[0, 1, z], ggt[1, 1, z]], Nl[gg[0, 1, z], ggt[1, 1, z]],
  Nl[ggt[0, 1, z], ggt[1, 1, z]], Nl[jj0[1, 1, z], jj0[1, 1, z]],
  Nl[jj0[1, 1, z], tt[0, 1, z]], Nl[jj0[1, 1, z], ggw[0, 1, z]],
  Nl[jj0[1, 1, z], ggwb[0, 1, z]], Nl[tt[0, 1, z], tt[0, 1, z]],
  Nl[tt[0, 1, z], ggw[0, 1, z]], Nl[tt[0, 1, z], ggwb[0, 1, z]],
  Nl[ggw[0, 1, z], ggwb[0, 1, z]],
  Nl[jj0[0, 1, z], jj0[0, 1, z], jj0[1, 1, z]],
  Nl[jj0[0, 1, z], jj0[0, 1, z], tt[0, 1, z]],
  Nl[jj0[0, 1, z], jj0[0, 1, z], ggw[0, 1, z]],
  Nl[jj0[0, 1, z], jj0[0, 1, z], ggwb[0, 1, z]],
  Nl[ww[0, 1, z], ww[0, 1, z], jj0[0, 1, z]],
  Nl[ww[0, 1, z], jj0[0, 1, z], wwb[0, 1, z]],
  Nl[ww[0, 1, z], jj0[0, 1, z], gg[0, 1, z]],
  Nl[ww[0, 1, z], jj0[0, 1, z], ggt[0, 1, z]],
  Nl[jj0[0, 1, z], wwb[0, 1, z], wwb[0, 1, z]],
  Nl[jj0[0, 1, z], wwb[0, 1, z], gg[0, 1, z]],
  Nl[jj0[0, 1, z], wwb[0, 1, z], ggt[0, 1, z]],
  Nl[jj0[0, 1, z], gg[0, 1, z], ggt[0, 1, z]],
  Nl[jj0[0, 1, z], jj0[0, 1, z], jj0[0, 1, z], jj0[0, 1, z]]
};

(* Run the original code *)
fcal=Checkbox[False];
lev=4;

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

NO[x_]:=x;

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

part=1/2IntegerPartitions[2 lev];

generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;
generated=DeleteDuplicates[generated];
actual=Table[If[Head[generated[[i]]]===NO,ReplaceAll[generated[[i]],NO->Nl],Nl[generated[[i]]]],{i,Length[generated]}];

(* Compare *)
Print["Expected count: ", Length[expected]];
Print["Actual count: ", Length[actual]];
Print["\n=== COMPARISON ==="];

If[Length[expected] == Length[actual],
  Print["✓ Counts match!"],
  Print["✗ Counts differ"]
];

If[Sort[expected] === Sort[actual],
  Print["✓✓✓ RESULTS MATCH PERFECTLY! ✓✓✓"],
  Print["✗ Results differ"];
  Print["Missing: ", Complement[expected, actual]];
  Print["Extra: ", Complement[actual, expected]]
];
