#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

ClearAll[beta, gamma, b, c, T, W, Wbar];

Bosonic[beta, gamma];
Fermionic[b, c];
OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[beta, gamma] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

T = -2 NO[b, Derivative[1][c]] - 3/2 NO[beta, Derivative[1][gamma]] - NO[Derivative[1][b], c] - 1/2 NO[Derivative[1][beta], gamma];
W = beta;
Wbar = (
  NO[beta, NO[beta, NO[gamma, NO[gamma, gamma]]]]
  + 2 NO[beta, NO[gamma, NO[gamma, NO[b, c]]]]
  - 4 NO[beta, NO[Derivative[1][gamma], gamma]]
  - 4/3 NO[gamma, NO[b, Derivative[1][c]]]
  + 2/3 NO[gamma, NO[Derivative[1][b], c]]
  + 2/3 NO[Derivative[1][beta], NO[gamma, gamma]]
  - 8/3 NO[Derivative[1][gamma], NO[b, c]]
  + 10/9 Derivative[2][gamma]
);

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

results = Association[
  "Wbar" -> serializeExpr[Wbar],
  "dWbar" -> serializeExpr[Derivative[1][Wbar]],
  "NO_W_dWbar" -> serializeExpr[NO[W, Derivative[1][Wbar]]],
  "NO_dW_Wbar" -> serializeExpr[NO[Derivative[1][W], Wbar]],
  "term5" -> serializeExpr[-9/4 NO[T, NO[T, NO[W, Derivative[1][Wbar]]]]],
  "term9" -> serializeExpr[15/4 NO[T, NO[T, NO[Derivative[1][W], Wbar]]]],
  "term5_plus_term9" -> serializeExpr[-9/4 NO[T, NO[T, NO[W, Derivative[1][Wbar]]]] + 15/4 NO[T, NO[T, NO[Derivative[1][W], Wbar]]]]
];

Export[
  "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/compare_term5_term9_opedefs.json",
  results,
  "JSON"
];

Print["Exported compare_term5_term9_opedefs.json"];
