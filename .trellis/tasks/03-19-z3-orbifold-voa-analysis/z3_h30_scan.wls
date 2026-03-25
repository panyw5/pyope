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

W = NO[beta, beta, beta];
Wb = NO[gamma, gamma, gamma];
Xp = NO[beta, beta, b];
Xm = NO[gamma, gamma, c'];
Yp = NO[beta', c'] - NO[beta, c''];
Ym = NO[gamma, b'] - NO[gamma', b];

showPoleData[name_, ope_] := Module[{},
  pp[name <> " maxpole", MaxPole[ope]];
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

n1 = NO[beta'', gamma];
n2 = NO[beta', gamma'];
n3 = NO[beta, gamma''];
n4 = NO[beta', beta, gamma, gamma];
n5 = NO[beta, beta, gamma', gamma];
n6 = NO[b', c'];
n7 = NO[b, c''];
n8 = NO[b, c', beta, gamma];
n9 = NO[beta, beta, beta, gamma, gamma, gamma];

Phi0 = t1 n1 + t2 n2 + t3 n3 + t4 n4 + t5 n5 + t6 n6 + t7 n7 + t8 n8 + t9 n9;

pp["Basis (h=3,q=0)", {n1, n2, n3, n4, n5, n6, n7, n8, n9}];

scanField["Phi0", Phi0, 3, 0];

pp["Descendant d2J", J''];
pp["Descendant dT", T'];
pp["Descendant JJ", NO[J, J]];
pp["Descendant JT", NO[J, T]];
pp["Descendant GGb", NO[G, Gb]];
pp["Descendant WWb", NO[W, Wb]];
