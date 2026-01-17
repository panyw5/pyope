#!/usr/bin/env wolframscript

(* Settings *)
lev = 4;

(* Definitions *)
opfields = {ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};
fields = {ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};

(* Functions *)
NO[x_] := x;

(* Derivative functions *)
Derivativer[x_] := D[x, z];
wwww[i_,j_,z_] := Nest[Derivativer, w, i];
j0j0[i_,j_,z_] := Nest[Derivativer, j0, i];
wbwb[i_,j_,z_] := Nest[Derivativer, wb, i];
tttt[i_,j_,z_] := Nest[Derivativer, t, i];
gggg[i_,j_,z_] := Nest[Derivativer, g, i];
gtgt[i_,j_,z_] := Nest[Derivativer, gt, i];
gwgw[i_,j_,z_] := Nest[Derivativer, gw, i];
gwbgwb[i_,j_,z_] := Nest[Derivativer, gwb, i];

(* Field Generation *)
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === ww[0,1,z]) := ww[gr-3/2,1,z];
gen[op_, gr_] /; (gr >= 1 && IntegerQ[gr] && op === jj0[0,1,z]) := jj0[gr-1,1,z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === wwb[0,1,z]) := wwb[gr-3/2,1,z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === tt[0,1,z]) := tt[gr-2,1,z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === gg[0,1,z]) := gg[gr-3/2,1,z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === ggt[0,1,z]) := ggt[gr-3/2,1,z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === ggw[0,1,z]) := ggw[gr-2,1,z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === ggwb[0,1,z]) := ggwb[gr-2,1,z];
gen[op_, gr_] /; ((!OddQ[2 gr] || gr < 3/2) && op === ww[0,1,z]) := Nothing;
gen[op_, gr_] /; (!IntegerQ[gr] && op === jj0[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!OddQ[2 gr] || gr < 3/2) && op === wwb[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!IntegerQ[gr] || gr < 2) && op === tt[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!OddQ[2 gr] || gr < 3/2) && op === gg[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!OddQ[2 gr] || gr < 3/2) && op === ggt[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!IntegerQ[gr] || gr < 2) && op === ggw[0,1,z]) := Nothing;
gen[op_, gr_] /; ((!IntegerQ[gr] || gr < 2) && op === ggwb[0,1,z]) := Nothing;

OperatorList[oplis_, grlis_] := Table[gen[op, gr], {gr, grlis}, {op, oplis}];

OperatorMult[a_, b_, right__] := NO[a, b, right];
OperatorMult[a_, b_] := NO[a, b];
OperatorMult[a_] := NO[a];

GenOp[partlis_] := Flatten@Table[Outer[OperatorMult, Sequence @@ OperatorList[opfields, partition]], {partition, partlis}];

(* Level - generate all field combinations *)
part = 1/2 IntegerPartitions[2 lev];
generated = Select[GenOp[part], !NumberQ[#] &] /. -1*x_ :> x;
DeleteDuplicates[%];
finalgen = Table[If[Head[%[[i]]] === NO, ReplaceAll[%[[i]], NO -> Nl], Nl[%[[i]]]], {i, Length[%]}];

(* Output the result *)
Print["Generated fields at level ", lev, ":"];
Print[finalgen];