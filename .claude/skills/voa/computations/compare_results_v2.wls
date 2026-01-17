#!/usr/bin/env wolframscript
(* Compare expected vs actual results - with order normalization *)

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

(* Actual result from our script *)
$IterationLimit = 100000;
lev = 4;
opfields = {ww[0,1,z], jj0[0,1,z], wwb[0,1,z], tt[0,1,z], gg[0,1,z], ggt[0,1,z], ggw[0,1,z], ggwb[0,1,z]};
NO[x_] := x;

gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === ww[0,1,z]) := ww[gr-3/2, 1, z];
gen[op_, gr_] /; (gr >= 1 && IntegerQ[gr] && op === jj0[0,1,z]) := jj0[gr-1, 1, z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === wwb[0,1,z]) := wwb[gr-3/2, 1, z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === tt[0,1,z]) := tt[gr-2, 1, z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === gg[0,1,z]) := gg[gr-3/2, 1, z];
gen[op_, gr_] /; (gr >= 3/2 && OddQ[2 gr] && op === ggt[0,1,z]) := ggt[gr-3/2, 1, z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === ggw[0,1,z]) := ggw[gr-2, 1, z];
gen[op_, gr_] /; (gr >= 2 && IntegerQ[gr] && op === ggwb[0,1,z]) := ggwb[gr-2, 1, z];
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

part = 1/2 IntegerPartitions[2 lev];
generated = Select[GenOp[part], !NumberQ[#] &] /. -1*x_ :> x;
generated = DeleteDuplicates[generated];
actual = Table[
  If[Head[generated[[i]]] === NO,
    generated[[i]] /. NO -> Nl,
    Nl[generated[[i]]]
  ],
  {i, Length[generated]}
];

(* Normalize Nl expressions by sorting arguments *)
normalizeNl[expr_] := Module[{args},
  If[Head[expr] === Nl,
    args = Sort[List @@ expr];
    Nl @@ args,
    expr
  ]
];

expectedNorm = normalizeNl /@ expected;
actualNorm = normalizeNl /@ actual;

Print["Expected count: ", Length[expected]];
Print["Actual count: ", Length[actual]];
Print["\n=== COMPARISON (with normalization) ==="];

If[Length[expected] == Length[actual],
  Print["✓ Counts match!"],
  Print["✗ Counts differ: expected ", Length[expected], ", got ", Length[actual]]
];

(* Check if all expected elements are in actual *)
missingInActual = Complement[expectedNorm, actualNorm];
extraInActual = Complement[actualNorm, expectedNorm];

If[Length[missingInActual] == 0,
  Print["✓ All expected elements found in actual result"],
  Print["✗ Missing ", Length[missingInActual], " elements in actual result:"];
  Print[missingInActual]
];

If[Length[extraInActual] == 0,
  Print["✓ No extra elements in actual result"],
  Print["✗ Found ", Length[extraInActual], " extra elements in actual result:"];
  Print[extraInActual]
];

If[Length[missingInActual] == 0 && Length[extraInActual] == 0,
  Print["\n✓✓✓ RESULTS MATCH PERFECTLY! ✓✓✓"],
  Print["\n✗✗✗ RESULTS DO NOT MATCH ✗✗✗"]
];
