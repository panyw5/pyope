#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

Unprotect[Print];
Print[_] := Null;
Protect[Print];

ClearAll[pp, showOPE];
pp[label_, expr_] := WriteString[$Output, label <> ": " <> ToString[InputForm[expr]] <> "\n"];

Bosonic[beta, gamma];
Fermionic[b, c];

OPE[beta, gamma] = MakeOPE[{-One}];
OPE[b, c] = MakeOPE[{One}];

J = NO[beta, gamma];
T = 1/2 (NO[beta', gamma] - NO[beta, gamma']) - NO[b, c'];
G = NO[b, gamma];
Gb = NO[beta, c'];
W = NO[beta, beta, beta];
Wb = NO[gamma, gamma, gamma];

N2a = NO[beta, beta, gamma, gamma];
N2b = NO[beta', gamma];
N2c = NO[beta, gamma'];
N2d = NO[b, c'];

N3a = NO[beta, beta, beta, gamma, gamma, gamma];
N3b = NO[beta', beta, gamma, gamma];
N3c = NO[beta, beta, gamma', gamma];
N3d = NO[beta, gamma, b, c'];
N3e = NO[beta'', gamma];
N3f = NO[beta', gamma'];
N3g = NO[beta, gamma''];
N3h = NO[b', c'];
N3i = NO[b, c''];

showOPE[name_, a_, b_] := Module[{ope},
  ope = OPESimplify[OPE[a, b], Expand];
  pp["=== " <> name <> " ===", ope];
  pp[name <> " pole5", OPEPole[5][ope]];
  pp[name <> " pole4", OPEPole[4][ope]];
  pp[name <> " pole3", OPEPole[3][ope]];
  pp[name <> " pole2", OPEPole[2][ope]];
  pp[name <> " pole1", OPEPole[1][ope]];
];

pp["Neutral weight-2 basis", {N2a,N2b,N2c,N2d}];
pp["Neutral weight-3 basis", {N3a,N3b,N3c,N3d,N3e,N3f,N3g,N3h,N3i}];

showOPE["W Wb", W, Wb];
showOPE["W N2a", W, N2a];
showOPE["W N2b", W, N2b];
showOPE["W N2c", W, N2c];
showOPE["W N2d", W, N2d];
showOPE["Wb N2a", Wb, N2a];
showOPE["Wb N2b", Wb, N2b];
showOPE["Wb N2c", Wb, N2c];
showOPE["Wb N2d", Wb, N2d];
showOPE["N2a N2a", N2a, N2a];
showOPE["N2a N2d", N2a, N2d];
showOPE["N2d N2d", N2d, N2d];
