#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

data = Import[
  "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/omega_weight4_monomials_opedefs.json",
  "RawJSON"
];

outPath = "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/omega_weight4_t_poles_opedefs.json";

ClearAll[betaff, gammaff, bff, cff, FG0, FG1, FG2, FG3, FG4, FG5, FG6, FG7, T];

Bosonic[betaff, gammaff];
Fermionic[bff, cff];
OPE[bff, cff] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[betaff, gammaff] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

Scan[
  Function[name,
    ToExpression[name <> " = " <> data["generator_defs"][name] <> ";"]
  ],
  Keys[data["generator_defs"]]
];

T = FG0;

operatorQ[expr_] := ! FreeQ[expr, NO | Derivative | betaff | gammaff | bff | cff];

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
  exprStr = item["expr"];
  coeff = ToExpression[item["coeff"]];
  mmaStr = item["mma"];
  op = coeff * ToExpression[mmaStr];
  results[exprStr] = <|
    "pole4" -> serializeExpr[OPE[T, op][4]],
    "pole3" -> serializeExpr[OPE[T, op][3]],
    "pole2" -> serializeExpr[OPE[T, op][2]],
    "pole1" -> serializeExpr[OPE[T, op][1]]
  |>,
  {item, data["terms"]}
];

Export[outPath, results, "JSON"];
Print["Exported omega_weight4_t_poles_opedefs.json"];
