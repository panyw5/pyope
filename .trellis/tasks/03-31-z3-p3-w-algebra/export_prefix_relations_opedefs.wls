#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

outPath = "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/prefix_relations_opedefs.json";

ClearAll[betaff, gammaff, bff, cff, SG0, SG1, SG2, SG3, SG4, SG5, SG6, SG7, terms, operatorQ, splitTerm, serializeExpr, results];

Bosonic[betaff, gammaff];
Fermionic[bff, cff];
OPE[bff, cff] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[betaff, gammaff] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

SG0 = -NO[Derivative[1][bff], cff]-2*NO[bff, Derivative[1][cff]]-3/2*NO[betaff, Derivative[1][gammaff]]-1/2*NO[Derivative[1][betaff], gammaff];
SG1 = 2*NO[bff, cff]+3*NO[betaff, gammaff];
SG2 = NO[gammaff, bff];
SG3 = 2*NO[Derivative[1][betaff], cff]+3*NO[betaff, Derivative[1][cff]];
SG4 = betaff;
SG5 = -4*NO[betaff, NO[Derivative[1][gammaff], gammaff]]+2*NO[betaff, NO[gammaff, NO[gammaff, NO[bff, cff]]]]-8/3*NO[Derivative[1][gammaff], NO[bff, cff]]-4/3*NO[gammaff, NO[bff, Derivative[1][cff]]]+2/3*NO[gammaff, NO[Derivative[1][bff], cff]]+2/3*NO[Derivative[1][betaff], NO[gammaff, gammaff]]+10/9*Derivative[2][gammaff]+NO[betaff, NO[betaff, NO[gammaff, NO[gammaff, gammaff]]]];
SG6 = bff;
SG7 = -3*NO[betaff, NO[betaff, NO[gammaff, NO[gammaff, Derivative[1][cff]]]]]-2*NO[Derivative[1][betaff], NO[betaff, NO[gammaff, NO[gammaff, cff]]]]+4*NO[betaff, NO[gammaff, NO[bff, NO[Derivative[1][cff], cff]]]]+4*NO[betaff, NO[gammaff, Derivative[2][cff]]]+4*NO[betaff, NO[Derivative[1][gammaff], Derivative[1][cff]]]-10/9*Derivative[3][cff]-8/3*NO[bff, NO[Derivative[2][cff], cff]]-2/3*NO[Derivative[2][betaff], NO[gammaff, cff]]+2/3*NO[Derivative[1][bff], NO[Derivative[1][cff], cff]]+8/3*NO[Derivative[1][betaff], NO[Derivative[1][gammaff], cff]];

terms = {
  NO[SG0, NO[SG0, NO[SG0, SG0]]],
  - 2/3 NO[SG0, NO[SG0, NO[SG0, NO[SG1, SG1]]]],
  1/3 NO[SG0, NO[SG0, NO[SG0, Derivative[1][SG1]]]],
  - 2/3 NO[SG0, NO[SG0, NO[SG1, NO[SG2, SG3]]]],
  1/2 NO[SG0, NO[SG0, NO[SG2, Derivative[1][SG3]]]],
  - 9/4 NO[SG0, NO[SG0, NO[SG4, Derivative[1][SG5]]]],
  - 1/12 NO[SG0, NO[SG0, NO[Derivative[1][SG1], NO[SG1, SG1]]]],
  - 1/3 NO[SG0, NO[SG0, NO[Derivative[1][SG1], Derivative[1][SG1]]]],
  - 1/6 NO[SG0, NO[SG0, NO[Derivative[1][SG2], SG3]]],
  15/4 NO[SG0, NO[SG0, NO[Derivative[1][SG4], SG5]]],
  3/4 NO[SG0, NO[SG0, NO[Derivative[2][SG1], SG1]]]
};

operatorQ[expr_] := ! FreeQ[expr, NO | Derivative | betaff | gammaff | bff | cff];

splitTerm[term_] := Module[{factors, scalar = 1, ops = {}},
  factors = If[Head[term] === Times, List @@ term, {term}];
  Do[If[operatorQ[f], AppendTo[ops, f], scalar = scalar * f], {f, factors}];
  <|"monomial" -> ToString[If[ops === {}, 1, Times @@ ops], InputForm], "coeff" -> ToString[scalar, InputForm]|>
];

serializeExpr[expr_] := Module[{expanded, parts},
  expanded = Expand[expr];
  If[expanded === 0, {}, parts = If[Head[expanded] === Plus, List @@ expanded, {expanded}]; splitTerm /@ parts]
];

results = Association[];
Do[
  results[ToString[n]] = serializeExpr[Total[Take[terms, n]]],
  {n, 1, Length[terms]}
];

Export[outPath, results, "JSON"];
Print["Exported prefix_relations_opedefs.json"];
