#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

ClearAll[beta, gamma, b, c, T, J, G, Gbar, W, GW, Wbar, GWbar];

Bosonic[beta, gamma];
Fermionic[b, c];

OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[beta, gamma] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

J = 2 NO[b, c] + 3 NO[beta, gamma];
G = NO[gamma, b];
Gbar = 2 NO[Derivative[1][beta], c] + 3 NO[beta, Derivative[1][c]];
T = -2 NO[b, Derivative[1][c]] - 3/2 NO[beta, Derivative[1][gamma]] - NO[Derivative[1][b], c] - 1/2 NO[Derivative[1][beta], gamma];
W = beta;
GW = b;
Wbar = NO[beta, NO[beta, NO[gamma, NO[gamma, gamma]]]]
  + 2 NO[beta, NO[gamma, NO[gamma, NO[b, c]]]]
  - 4 NO[beta, NO[Derivative[1][gamma], gamma]]
  - 4/3 NO[gamma, NO[b, Derivative[1][c]]]
  + 2/3 NO[gamma, NO[Derivative[1][b], c]]
  + 2/3 NO[Derivative[1][beta], NO[gamma, gamma]]
  - 8/3 NO[Derivative[1][gamma], NO[b, c]]
  + 10/9 Derivative[2][gamma];
GWbar = -10/9 Derivative[3][c]
  - 8/3 NO[b, NO[Derivative[2][c], c]]
  - 3 NO[beta, NO[beta, NO[gamma, NO[gamma, Derivative[1][c]]]]]
  + 4 NO[beta, NO[gamma, NO[b, NO[Derivative[1][c], c]]]]
  + 4 NO[beta, NO[gamma, Derivative[2][c]]]
  + 4 NO[beta, NO[Derivative[1][gamma], Derivative[1][c]]]
  - 2/3 NO[Derivative[2][beta], NO[gamma, c]]
  + 2/3 NO[Derivative[1][b], NO[Derivative[1][c], c]]
  - 2 NO[Derivative[1][beta], NO[beta, NO[gamma, NO[gamma, c]]]]
  + 8/3 NO[Derivative[1][beta], NO[Derivative[1][gamma], c]];

terms = <|
  "term0_T4" -> NO[T, NO[T, NO[T, T]]],
  "term1_T3JJ" -> -2/3 NO[T, NO[T, NO[T, NO[J, J]]]],
  "term2_T3dJ" -> 1/3 NO[T, NO[T, NO[T, Derivative[1][J]]]],
  "term3_T2JGGB" -> -2/3 NO[T, NO[T, NO[J, NO[G, Gbar]]]],
  "term4_T2GdGB" -> 1/2 NO[T, NO[T, NO[G, Derivative[1][Gbar]]]],
  "term5_T2WdWB" -> -9/4 NO[T, NO[T, NO[W, Derivative[1][Wbar]]]],
  "term6_T2dJJJ" -> -1/12 NO[T, NO[T, NO[Derivative[1][J], NO[J, J]]]],
  "term7_T2dJdJ" -> -1/3 NO[T, NO[T, NO[Derivative[1][J], Derivative[1][J]]]],
  "term8_T2dGGB" -> -1/6 NO[T, NO[T, NO[Derivative[1][G], Gbar]]],
  "term9_T2dWWB" -> 15/4 NO[T, NO[T, NO[Derivative[1][W], Wbar]]],
  "term10_T2d2JJ" -> 3/4 NO[T, NO[T, NO[Derivative[2][J], J]]]
|>;

operatorQ[expr_] := ! FreeQ[expr, NO | Derivative | beta | gamma | b | c];

splitTerm[term_] := Module[{factors, scalar = 1, ops = {}},
  factors = If[Head[term] === Times, List @@ term, {term}];
  Do[
    If[operatorQ[f], AppendTo[ops, f], scalar = scalar * f],
    {f, factors}
  ];
  <|
    "monomial" -> ToString[If[ops === {}, 1, Times @@ ops], InputForm],
    "coeff" -> ToString[scalar, InputForm]
  |>
];

serializeExpr[expr_] := Module[{expanded, parts},
  expanded = Expand[expr];
  If[expanded === 0,
    {},
    parts = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
    splitTerm /@ parts
  ]
];

results = Association[];
Do[
  results[key] = serializeExpr[terms[key]],
  {key, Keys[terms]}
];

Export[
  "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/compare_nt_terms_opedefs.json",
  results,
  "JSON"
];

Print["Exported compare_nt_terms_opedefs.json"];
