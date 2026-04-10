AppendTo[$Path, Directory[]];

Quiet[
  Get["OPEdefs`"],
  {FrontEndObject::notavail, Get::noopen, Needs::nocont}
];

ClearAll[
  PySymbolName,
  PyEncode,
  PyEncodeChunkedList,
  PyEncodeAssociation,
  PyEncodeNOChain,
  PyEncodeProduct,
  PyProtocolError,
  PyReadPayload,
  PyReadPositiveIntegerPayload,
  PyWriteResult,
  PyScalarNormalize,
  PyNonScalarQ,
  PyOperatorSymbolQ,
  PyRegisteredOperatorQ,
  PyEvalExpr,
  PyEvalExprList,
  PyEvaluateValue,
  PyCanonicalizeExpr,
  PyCanonicalizeExprList,
  PyCanonicalizeValue,
  PyMain
];

PyProtocolError[expr_] := Throw[
  "Unsupported Wolfram bridge expression: " <> ToString[Unevaluated[expr], InputForm],
  "PYOPE_PROTOCOL_ERROR"
];

PyRegisteredOperatorQ[expr_] := Module[{raw, names, name},
  raw = PyReadPayload["PYOPE_WL_OPERATORS"];
  names = If[StringLength[raw] > 0, StringSplit[raw, ","], {}];
  name = ToString[Unevaluated[expr], InputForm];
  MemberQ[names, name]
];

PyOperatorSymbolQ[expr_] := PyRegisteredOperatorQ[expr] || Head[Unevaluated[expr]] === Operator;
PyDict[list_List] := PyDict[Association[Rule @@@ list]];

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

PySymbolName[expr_] := ToString[Unevaluated[expr], InputForm];

PyEncode[Derivative[n_Integer][x_]] := "dn(" <> ToString[n, InputForm] <> ", " <> PyEncode[x] <> ")";
PyEncodeNOChain[expr_NO] := Module[{items = {}, current = expr},
  While[Head[current] === NO,
    AppendTo[items, current[[1]]];
    current = current[[2]];
  ];
  AppendTo[items, current];
  If[Length[items] >= 3,
    "PyNOChain([" <> StringRiffle[PyEncode /@ items, ", "] <> "])",
    "NO(" <> PyEncode[expr[[1]]] <> ", " <> PyEncode[expr[[2]]] <> ")"
  ]
];
PyEncode[NO[a_, b_]] := PyEncodeNOChain[NO[a, b]];
PyEncode[OPEData[a_List]] := "MakeOPE([" <> StringRiffle[PyEncode /@ a, ", "] <> "])";
PyEncode[PyTuple[list_List]] := "PyTuple([" <> StringRiffle[PyEncode /@ list, ", "] <> "])";
PyEncode[PyDict[assoc_Association]] := PyEncodeAssociation[assoc];
PyEncode[list_List] := "[" <> StringRiffle[PyEncode /@ list, ", "] <> "]";
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
PyEncode[str_String] := ToString[InputForm[str]];
PyEncode[expr_ /; PyOperatorSymbolQ[Unevaluated[expr]]] := PySymbolName[Unevaluated[expr]];
PyEncode[expr_Symbol] := PySymbolName[expr];
PyEncode[expr_] := PyProtocolError[expr];

PyEncodeChunkedList[list_List, chunkSize_Integer:64] := Module[
  {safeChunkSize, chunks, encodedChunks},
  safeChunkSize = Max[1, chunkSize];
  chunks = Partition[list, UpTo[safeChunkSize]];
  encodedChunks = PyEncode /@ chunks;
  "PYOPE_CHUNKED_LIST\n" <> StringRiffle[encodedChunks, "\n"]
];

PyEncodeAssociation[assoc_Association] := Module[{items},
  items = KeyValueMap[
    "(" <> PyEncode[#1] <> ", " <> PyEncode[#2] <> ")" &,
    assoc
  ];
  "PyDict([" <> StringRiffle[items, ", "] <> "])"
];

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

PyReadPayload[name_String] := Module[{path, stream, text},
  path = Environment[name <> "_PATH"];
  If[!StringQ[path], path = ""];
  If[StringLength[path] == 0,
    text = Environment[name];
    If[!StringQ[text], text = ""];
    Return[text]
  ];
  stream = OpenRead[path, BinaryFormat -> False, CharacterEncoding -> "UTF-8"];
  If[stream === $Failed, PyProtocolError["unable to open payload file for " <> name]];
  text = ReadString[stream];
  Close[stream];
  If[!StringQ[text], text = ""];
  text
];

PyReadPositiveIntegerPayload[name_String, default_Integer] := Module[{raw, parsed},
  raw = StringTrim[PyReadPayload[name]];
  If[StringLength[raw] == 0, Return[default]];
  parsed = Quiet[Check[ToExpression[raw], $Failed]];
  If[IntegerQ[parsed] && parsed > 0, parsed, default]
];

PyWriteResult[text_String] := Module[{path, stream},
  path = Environment["PYOPE_WL_RESULT_PATH"];
  If[!StringQ[path], path = ""];
  If[StringLength[path] == 0,
    Print["PYOPE_RESULT: " <> text];
    Return[text]
  ];
  stream = OpenWrite[path, BinaryFormat -> False, CharacterEncoding -> "UTF-8"];
  If[stream === $Failed, PyProtocolError["unable to open result payload file"]];
  WriteString[stream, text];
  Close[stream];
  text
];

PyApplyRegistry[] := Module[{code = PyReadPayload["PYOPE_WL_REGISTRY"]},
  If[StringLength[code] > 0, ToExpression[code]];
];

PyEvalExpr[] := Module[{exprText, expr},
  exprText = PyReadPayload["PYOPE_WL_EXPR"];
  If[StringLength[exprText] == 0, PyProtocolError["missing PYOPE_WL_EXPR"]];
  expr = ToExpression[exprText];
  PyEvaluateValue[expr]
];

PyEvaluateValue[value_List] := PyEvaluateValue /@ value;
PyEvaluateValue[PyTuple[list_List]] := PyTuple[PyEvaluateValue /@ list];
PyEvaluateValue[PyDict[assoc_Association]] := PyDict[Map[PyEvaluateValue, assoc]];
PyEvaluateValue[expr_] := Expand[expr];

PyEvalExprList[] := Module[{exprText, exprs},
  exprText = PyReadPayload["PYOPE_WL_EXPR_LIST"];
  If[StringLength[exprText] == 0, PyProtocolError["missing PYOPE_WL_EXPR_LIST"]];
  exprs = ToExpression[exprText];
  If[Head[exprs] =!= List, PyProtocolError["PYOPE_WL_EXPR_LIST must decode to List"]];
  PyEvaluateValue /@ exprs
];

PyCanonicalizeExpr[expr_] := OPESimplify[Expand[expr]];

PyCanonicalizeValue[value_List] := PyCanonicalizeValue /@ value;
PyCanonicalizeValue[PyTuple[list_List]] := PyTuple[PyCanonicalizeValue /@ list];
PyCanonicalizeValue[PyDict[assoc_Association]] := PyDict[Map[PyCanonicalizeValue, assoc]];
PyCanonicalizeValue[expr_] := PyCanonicalizeExpr[expr];

PyCanonicalizeExprList[] := Module[{exprText, exprs},
  exprText = PyReadPayload["PYOPE_WL_EXPR_LIST"];
  If[StringLength[exprText] == 0, PyProtocolError["missing PYOPE_WL_EXPR_LIST"]];
  exprs = ToExpression[exprText];
  If[Head[exprs] =!= List, PyProtocolError["PYOPE_WL_EXPR_LIST must decode to List"]];
  PyCanonicalizeValue /@ exprs
];

PyMain[] := Module[{op, left, right, result, encoded},
  Catch[
    op = PyReadPayload["PYOPE_WL_OPERATION"];
    PyApplyRegistry[];
    result = Switch[op,
      "EVAL", PyEvalExpr[],
      "EVAL_LIST", PyEvalExprList[],
      "EXPAND", PyEvalExpr[],
      "CANONICALIZE_LIST", PyCanonicalizeExprList[],
      _,
      left = ToExpression[PyReadPayload["PYOPE_WL_LEFT"]];
      right = ToExpression[PyReadPayload["PYOPE_WL_RIGHT"]];
      Switch[op,
      "OPE", OPE[left, right],
      "NO", NO[left, right],
      _, $Failed
      ]
    ];
    encoded = If[
      StringMatchQ[op, "*_LIST"] && Head[result] === List,
      PyEncodeChunkedList[
        result,
        PyReadPositiveIntegerPayload["PYOPE_WL_RESULT_CHUNK_SIZE", 64]
      ],
      PyEncode[result]
    ];
    PyWriteResult[encoded],
    "PYOPE_PROTOCOL_ERROR",
    (Print["PYOPE_PROTOCOL_ERROR: " <> #1]; Exit[3]) &
  ]
];

PyMain[];
