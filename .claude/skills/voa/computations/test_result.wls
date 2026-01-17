#!/usr/bin/env wolframscript

SetDirectory["/Users/lelouch/pyope/.claude/skills/voa/computations"];
Get["ope.math.m"];
Get["OPEDefs.m"];

lev = 4;

(* Definitions *)
opfields = {ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};
fields = {ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};

Bosonic[\[Beta]1,\[Gamma]1,ww[i_,j_,z],jj0[i_,j_,z],wwb[i_,j_,z],tt[i_,j_,z]];
Fermionic[b1,c1,gg[i_,j_,z],ggt[i_,j_,z],ggw[i_,j_,z],ggwb[i_,j_,z]];

OPE[\[Beta]1,\[Gamma]1] = MakeOPE[-One (z-w)^-1+Ord[z,w,0]];
OPE[b1,c1] = MakeOPE[One (z-w)^-1+Ord[z,w,0]];

w = \[Beta]1;
j0 = 2 NO[b1,c1]+3 NO[\[Beta]1,\[Gamma]1];
wb = NO[\[Beta]1,NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,\[Gamma]1]]]]+2 NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,NO[b1,c1]]]]-4 NO[\[Beta]1,NO[Derivative[1][\[Gamma]1],\[Gamma]1]]-4/3 NO[\[Gamma]1,NO[b1,Derivative[1][c1]]]+2/3 NO[\[Gamma]1,NO[Derivative[1][b1],c1]]+2/3 NO[Derivative[1][\[Beta]1],NO[\[Gamma]1,\[Gamma]1]]-8/3 NO[Derivative[1][\[Gamma]1],NO[b1,c1]]+(10 \[Gamma]1'')/9;
t = -2 NO[b1,Derivative[1][c1]]-3/2 NO[\[Beta]1,Derivative[1][\[Gamma]1]]-NO[Derivative[1][b1],c1]-1/2 NO[Derivative[1][\[Beta]1],\[Gamma]1];
g = NO[\[Gamma]1,b1];
gt = 2NO[\[Beta]1',c1]+3NO[\[Beta]1,c1'];
gw = b1;
gwb = 8/3 NO[b1,NO[c1'',c1]]+3 NO[\[Beta]1,NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,Derivative[1][c1]]]]]-4 NO[\[Beta]1,NO[\[Gamma]1,NO[b1,NO[Derivative[1][c1],c1]]]]-4 NO[\[Beta]1,NO[\[Gamma]1,c1'']]-4 NO[\[Beta]1,NO[Derivative[1][\[Gamma]1],Derivative[1][c1]]]-2/3 NO[Derivative[1][b1],NO[Derivative[1][c1],c1]]+2 NO[Derivative[1][\[Beta]1],NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,c1]]]]-8/3 NO[Derivative[1][\[Beta]1],NO[Derivative[1][\[Gamma]1],c1]]+2/3 NO[\[Beta]1'',NO[\[Gamma]1,c1]]+(10 c1''')/9;

(* Field generation rules *)
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

(* Operator functions *)
OperatorList[oplis_, grlis_] := Table[gen[op, gr], {gr, grlis}, {op, oplis}];

OperatorMult[a_, b_, right__] := NO[a, b, right];
OperatorMult[a_, b_] := NO[a, b];
OperatorMult[a_] := NO[a];

GenOp[partlis_] := Flatten@Table[Outer[OperatorMult, Sequence @@ OperatorList[opfields, partition]], {partition, partlis}];

(* Generate partitions *)
part = 1/2 IntegerPartitions[2 lev];
Print["part = ", part];

(* Generate operators step by step *)
generated = Select[GenOp[part], !NumberQ[#] &] /. -1*x_ :> x;
generatedNoDup = DeleteDuplicates[generated];
Print["Generated unique operators: Length = ", Length[generatedNoDup]];
Print[generatedNoDup];

(* Convert to final format *)
finalgen = Table[If[Head[generatedNoDup[[i]]] === NO, ReplaceAll[generatedNoDup[[i]], NO -> Nl], Nl[generatedNoDup[[i]]]], {i, Length[generatedNoDup]}];
Print["Final result (Nl notation):"];
Print[finalgen];

(* Compare with expected result *)
expected = {Nl[jj0[3, 1, z]], Nl[tt[2, 1, z]], Nl[ggw[2, 1, z]], Nl[ggwb[2, 1, z]], Nl[jj0[0, 1, z], jj0[2, 1, z]], Nl[jj0[0, 1, z], tt[1, 1, z]], Nl[jj0[0, 1, z], ggw[1, 1, z]], Nl[jj0[0, 1, z], ggwb[1, 1, z]], Nl[ww[0, 1, z], ww[1, 1, z]], Nl[ww[1, 1, z], wwb[0, 1, z]], Nl[ww[1, 1, z], gg[0, 1, z]], Nl[ww[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], wwb[1, 1, z]], Nl[wwb[0, 1, z], wwb[1, 1, z]], Nl[wwb[1, 1, z], gg[0, 1, z]], Nl[wwb[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], gg[1, 1, z]], Nl[wwb[0, 1, z], gg[1, 1, z]], Nl[gg[0, 1, z], gg[1, 1, z]], Nl[gg[1, 1, z], ggt[0, 1, z]], Nl[ww[0, 1, z], ggt[1, 1, z]], Nl[wwb[0, 1, z], ggt[1, 1, z]], Nl[gg[0, 1, z], ggt[1, 1, z]], Nl[ggt[0, 1, z], ggt[1, 1, z]], Nl[jj0[1, 1, z], jj0[1, 1, z]], Nl[jj0[1, 1, z], tt[0, 1, z]], Nl[jj0[1, 1, z], ggw[0, 1, z]], Nl[jj0[1, 1, z], ggwb[0, 1, z]], Nl[tt[0, 1, z], tt[0, 1, z]], Nl[tt[0, 1, z], ggw[0, 1, z]], Nl[tt[0, 1, z], ggwb[0, 1, z]], Nl[ggw[0, 1, z], ggwb[0, 1, z]], Nl[jj0[0, 1, z], jj0[0, 1, z], jj0[1, 1, z]], Nl[jj0[0, 1, z], jj0[0, 1, z], tt[0, 1, z]], Nl[jj0[0, 1, z], jj0[0, 1, z], ggw[0, 1, z]], Nl[jj0[0, 1, z], jj0[0, 1, z], ggwb[0, 1, z]], Nl[ww[0, 1, z], ww[0, 1, z], jj0[0, 1, z]], Nl[ww[0, 1, z], jj0[0, 1, z], wwb[0, 1, z]], Nl[ww[0, 1, z], jj0[0, 1, z], gg[0, 1, z]], Nl[ww[0, 1, z], jj0[0, 1, z], ggt[0, 1, z]], Nl[jj0[0, 1, z], wwb[0, 1, z], wwb[0, 1, z]], Nl[jj0[0, 1, z], wwb[0, 1, z], gg[0, 1, z]], Nl[jj0[0, 1, z], wwb[0, 1, z], ggt[0, 1, z]], Nl[jj0[0, 1, z], gg[0, 1, z], ggt[0, 1, z]], Nl[jj0[0, 1, z], jj0[0, 1, z], jj0[0, 1, z], jj0[0, 1, z]]};

Print["Expected result: Length = ", Length[expected]];
Print[expected];

(* Check if the result matches *)
Print["Checking comparison:"];
Print["Length match: ", Length[finalgen] == Length[expected]];

(* Check for elements in expected but not in finalgen *)
missing = Complement[expected, finalgen];
extra = Complement[finalgen, expected];
Print["Elements in expected but not in result: ", Length[missing]];
Print[missing];
Print["Elements in result but not in expected: ", Length[extra]];
Print[extra];