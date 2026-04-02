#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

data = Import[
  "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/weight8_candidates_for_opedefs.json",
  "RawJSON"
];

cachePath = "/Users/lelouch/pyope/.trellis/tasks/03-31-z3-p3-w-algebra/weight8_candidates_opedefs_cache.json";

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

cache = If[FileExistsQ[cachePath], Import[cachePath, "RawJSON"], <||>];
entries = If[KeyExistsQ[cache, "entries"], Association[cache["entries"]], <||>];

candidates = data["candidates"];
remaining = Select[candidates, ! KeyExistsQ[entries, #["expr"]] &];

Print["cached candidates: " <> ToString[Length[entries]] <> "/" <> ToString[Length[candidates]]];

Do[
  exprStr = item["expr"];
  mmaStr = item["mma"];
  Export[cachePath, <|"entries" -> entries, "in_progress" -> exprStr|>, "JSON"];
  Print["starting: " <> exprStr];
  timing = AbsoluteTiming[result = serializeExpr[ToExpression[mmaStr]];][[1]];
  entries[exprStr] = <|"seconds" -> timing, "terms" -> result|>;
  Export[cachePath, <|"entries" -> entries, "last_completed" -> exprStr|>, "JSON"];
  Print["done in " <> ToString[timing] <> " sec"];
  ,
  {item, remaining}
];

Export[cachePath, <|"entries" -> entries, "completed" -> True|>, "JSON"];
Print["completed all candidates"];
