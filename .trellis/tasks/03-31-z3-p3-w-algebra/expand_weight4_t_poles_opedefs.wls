#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

data = Import[
  "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/weight4_quasiprimary_candidates.json",
  "RawJSON"
];

outPath = "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/weight4_t_poles_opedefs.json";

ClearAll[betaff, gammaff, bff, cff, SG0, SG1, SG2, SG3, SG4, SG5, SG6, SG7];

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

T = SG0;

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
  mmaStr = item["mma"];
  op = ToExpression[mmaStr];
  results[exprStr] = <|
    "pole4" -> serializeExpr[OPE[T, op][4]],
    "pole3" -> serializeExpr[OPE[T, op][3]],
    "pole2" -> serializeExpr[OPE[T, op][2]],
    "pole1" -> serializeExpr[OPE[T, op][1]]
  |>,
  {item, data["candidates"]}
];

Export[outPath, results, "JSON"];
Print["Exported weight4_t_poles_opedefs.json"];
