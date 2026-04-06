AppendTo[$Path, Directory[]];

Quiet[
  Get["OPEdefs`"],
  {FrontEndObject::notavail, Get::noopen, Needs::nocont}
];

ClearAll[
  PySymbolName,
  PyEncode,
  PyEncodeProduct,
  PyProtocolError,
  PyScalarNormalize,
  PyNonScalarQ,
  PyOperatorSymbolQ,
  PyRegisteredOperatorQ,
  PyEvalExpr,
  PyMain
];

PyProtocolError[expr_] := Throw[
  "Unsupported Wolfram bridge expression: " <> ToString[Unevaluated[expr], InputForm],
  "PYOPE_PROTOCOL_ERROR"
];

PyRegisteredOperatorQ[s_Symbol] := Module[{raw, names},
  raw = Environment["PYOPE_WL_OPERATORS"];
  names = If[StringLength[raw] > 0, StringSplit[raw, ","], {}];
  MemberQ[names, SymbolName[Unevaluated[s]]]
];

PyOperatorSymbolQ[s_Symbol] := PyRegisteredOperatorQ[s] || Head[Unevaluated[s]] === Operator;

PyNonScalarQ[expr_] := MatchQ[expr, _NO | Derivative[_Integer][_]];
PyNonScalarQ[expr_Symbol] := PyOperatorSymbolQ[expr];
PyNonScalarQ[expr_Times] := Count[List @@ expr, factor_ /; PyNonScalarQ[factor]] >= 1;
PyNonScalarQ[_] := False;

PyScalarNormalize[expr_] := FixedPoint[
  Simplify[# /. {
    If[test_, a_, b_] :> If[TrueQ[test], a, b],
    Power[If[test_, a_, b_], n_Integer] :> If[TrueQ[test], a, b]^n
  }] &,
  expr
];

PySymbolName[s_Symbol] := SymbolName[Unevaluated[s]];

PyEncode[Derivative[n_Integer][x_]] := "dn(" <> ToString[n, InputForm] <> ", " <> PyEncode[x] <> ")";
PyEncode[NO[a_, b_]] := "NO(" <> PyEncode[a] <> ", " <> PyEncode[b] <> ")";
PyEncode[OPEData[a_List]] := "MakeOPE([" <> StringRiffle[PyEncode /@ a, ", "] <> "])";
PyEncode[expr_Plus] := "(" <> StringRiffle[PyEncode /@ (List @@ expr), " + "] <> ")";
PyEncode[expr_Times] := PyEncodeProduct[List @@ expr, expr];
PyEncode[expr_Rational] := "(" <> ToString[Numerator[expr], InputForm] <> "/" <> ToString[Denominator[expr], InputForm] <> ")";
PyEncode[expr_Integer] := ToString[expr, InputForm];
PyEncode[expr_Real] := Module[{r = Rationalize[expr, 0]},
  If[Head[r] === Rational,
    PyEncode[r],
    ToString[expr, InputForm]
  ]
];
PyEncode[0] := "Zero";
PyEncode[One] := "One";
PyEncode[expr_Symbol] := PySymbolName[expr];
PyEncode[expr_] := PyProtocolError[expr];

PyEncodeProduct[factors_List, original_] := Module[{nonscalar, scalarFactors, nonscalarFactor},
  nonscalar = Select[factors, PyNonScalarQ];
  If[Length[nonscalar] > 1, PyProtocolError[original]];
  If[Length[nonscalar] == 0,
    Return["(" <> StringRiffle[PyEncode /@ (PyScalarNormalize /@ factors), " * "] <> ")"]
  ];
  nonscalarFactor = First[nonscalar];
  scalarFactors = PyScalarNormalize /@ Select[factors, Not @* PyNonScalarQ];
  If[Length[scalarFactors] == 0, Return[PyEncode[nonscalarFactor]]];
  "(" <> StringRiffle[Join[PyEncode /@ scalarFactors, {PyEncode[nonscalarFactor]}], " * "] <> ")"
];

PyApplyRegistry[] := Module[{code = Environment["PYOPE_WL_REGISTRY"]},
  If[StringLength[code] > 0, ToExpression[code]];
];

PyEvalExpr[] := Module[{exprText, expr},
  exprText = Environment["PYOPE_WL_EXPR"];
  If[StringLength[exprText] == 0, PyProtocolError["missing PYOPE_WL_EXPR"]];
  expr = ToExpression[exprText];
  Expand[expr]
];

PyMain[] := Module[{op, left, right, result, encoded},
  Catch[
    op = Environment["PYOPE_WL_OPERATION"];
    PyApplyRegistry[];
    result = Switch[op,
      "EVAL", PyEvalExpr[],
      "EXPAND", PyEvalExpr[],
      _,
      left = ToExpression[Environment["PYOPE_WL_LEFT"]];
      right = ToExpression[Environment["PYOPE_WL_RIGHT"]];
      Switch[op,
      "OPE", OPE[left, right],
      "NO", NO[left, right],
      _, $Failed
      ]
    ];
    encoded = PyEncode[result];
    Print["PYOPE_RESULT: " <> encoded],
    "PYOPE_PROTOCOL_ERROR",
    (Print["PYOPE_PROTOCOL_ERROR: " <> #1]; Exit[3]) &
  ]
];

PyMain[];
