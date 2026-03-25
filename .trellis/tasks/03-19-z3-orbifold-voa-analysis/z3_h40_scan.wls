#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

Unprotect[Print];
Print[_] := Null;
Protect[Print];

ClearAll[pp, showPoleData, scanField];
pp[label_, expr_] := WriteString[$Output, label <> ": " <> ToString[InputForm[expr]] <> "\n"];

Bosonic[beta, gamma];
Fermionic[b, c];

OPE[beta, gamma] = MakeOPE[{-One}];
OPE[b, c] = MakeOPE[{One}];

J = NO[beta, gamma];
T = 1/2 (NO[beta', gamma] - NO[beta, gamma']) - NO[b, c'];
G = NO[b, gamma];
Gb = NO[beta, c'];

showPoleData[name_, ope_] := Module[{},
  pp[name <> " maxpole", MaxPole[ope]];
  pp[name <> " pole6", OPEPole[6][ope]];
  pp[name <> " pole5", OPEPole[5][ope]];
  pp[name <> " pole4", OPEPole[4][ope]];
  pp[name <> " pole3", OPEPole[3][ope]];
  pp[name <> " pole2", OPEPole[2][ope]];
  pp[name <> " pole1", OPEPole[1][ope]];
];

scanField[name_, phi_, h_, q_] := Module[{oj, ot, og, ogb},
  pp["", ""];
  pp["################", name];
  pp[name <> " field", phi];
  oj = OPESimplify[OPE[J, phi], Expand];
  ot = OPESimplify[OPE[T, phi], Expand];
  og = OPESimplify[OPE[G, phi], Expand];
  ogb = OPESimplify[OPE[Gb, phi], Expand];
  showPoleData[name <> " J", oj];
  pp[name <> " J pole1 - q phi", OPESimplify[OPEPole[1][oj] - q phi, Expand]];
  showPoleData[name <> " T", ot];
  pp[name <> " T pole2 - h phi", OPESimplify[OPEPole[2][ot] - h phi, Expand]];
  pp[name <> " T pole1 - d phi", OPESimplify[OPEPole[1][ot] - phi', Expand]];
  showPoleData[name <> " G", og];
  showPoleData[name <> " Gb", ogb];
];

v1 = NO[beta''', gamma];
v2 = NO[beta'', gamma'];
v3 = NO[beta', gamma''];
v4 = NO[beta, gamma'''];
v5 = NO[beta'', beta, gamma, gamma];
v6 = NO[beta', beta', gamma, gamma];
v7 = NO[beta', beta, gamma', gamma];
v8 = NO[beta, beta, gamma'', gamma];
v9 = NO[beta, beta, gamma', gamma'];
v10 = NO[beta, beta, beta, gamma, gamma, gamma, gamma];
v11 = NO[b'', c'];
v12 = NO[b', c''];
v13 = NO[b, c'''];
v14 = NO[beta', gamma, b, c'];
v15 = NO[beta, gamma', b, c'];
v16 = NO[beta, gamma, b', c'];
v17 = NO[beta, gamma, b, c''];
v18 = NO[beta, beta, gamma, gamma, b, c'];

Phi0 = s1*v1 + s2*v2 + s3*v3 + s4*v4 + s5*v5 + s6*v6 + s7*v7 + s8*v8 + s9*v9 + s10*v10 + s11*v11 + s12*v12 + s13*v13 + s14*v14 + s15*v15 + s16*v16 + s17*v17 + s18*v18;

pp["Basis (h=4, neutral)", {v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18}];

scanField["Phi0", Phi0, 4, 0];

pp["Descendant d3J", J'''];
pp["Descendant d2T", T''];
pp["Composite JJ", NO[J, J]];
pp["Composite JT", NO[J, T]];
pp["Composite TT", NO[T, T]];
