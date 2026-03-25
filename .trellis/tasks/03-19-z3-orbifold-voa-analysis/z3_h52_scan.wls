#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

Unprotect[Print];
Print[_] := Null;
Protect[Print];

ClearAll[pp, showPoleData, showOPE, scanField];
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
  pp[name <> " pole4", OPEPole[4][ope]];
  pp[name <> " pole3", OPEPole[3][ope]];
  pp[name <> " pole2", OPEPole[2][ope]];
  pp[name <> " pole1", OPEPole[1][ope]];
];

showOPE[name_, a_, b_] := Module[{ope},
  ope = OPESimplify[OPE[a, b], Expand];
  pp["=== " <> name <> " ===", ope];
  showPoleData[name, ope];
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

phiP1 = beta'';
phiP2 = NO[beta', beta, gamma];
phiP3 = NO[beta, beta, gamma'];
phiP4 = NO[b, c', beta];

phiM1 = gamma'';
phiM2 = NO[gamma', beta, gamma];
phiM3 = NO[gamma, gamma, beta'];
phiM4 = NO[b, c', gamma];

PhiP = a1 phiP1 + a2 phiP2 + a3 phiP3 + a4 phiP4;
PhiM = m1 phiM1 + m2 phiM2 + m3 phiM3 + m4 phiM4;

pp["Basis (h=5/2,q=+1)", {phiP1, phiP2, phiP3, phiP4}];
pp["Basis (h=5/2,q=-1)", {phiM1, phiM2, phiM3, phiM4}];

scanField["PhiP", PhiP, 5/2, 1];
scanField["PhiM", PhiM, 5/2, -1];

showOPE["Gb Xp", Gb, Xp];
showOPE["G Xm", G, Xm];
showOPE["T Xp", T, Xp];
showOPE["T Xm", T, Xm];
showOPE["J Xp", J, Xp];
showOPE["J Xm", J, Xm];
showOPE["G Yp", G, Yp];
showOPE["Gb Ym", Gb, Ym];
