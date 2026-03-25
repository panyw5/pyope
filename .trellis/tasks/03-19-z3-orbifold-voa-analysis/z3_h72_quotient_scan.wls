#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

Unprotect[Print];
Print[_] := Null;
Protect[Print];

ClearAll[pp, showPoleData, scanField, expandBasis];
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

u1 = NO[beta''', beta, gamma];
u2 = NO[beta'', beta, gamma'];
u3 = NO[beta'', beta', gamma];
u4 = NO[beta', beta, gamma''];
u5 = NO[beta', beta', gamma'];
u6 = NO[beta', beta, beta, gamma, gamma];
u7 = NO[beta, beta, beta, gamma', gamma];
u8 = NO[beta, beta, b', c'];
u9 = NO[beta, beta, b, c''];
u10 = NO[beta', beta, b, c'];
u11 = NO[beta, gamma, gamma, b, c'];

v1 = Xp';
v2 = NO[J, Xp];
v3 = NO[W, Ym];
v4 = NO[Yp, J];

pp["Unique-ish basis (h=7/2,q=+1)", {u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11}];
pp["Known descendant/composite span", {v1,v2,v3,v4}];

PhiRaw = r1*u1 + r2*u2 + r3*u3 + r4*u4 + r5*u5 + r6*u6 + r7*u7 + r8*u8 + r9*u9 + r10*u10 + r11*u11;
PhiQuot = PhiRaw - q1*v1 - q2*v2 - q3*v3 - q4*v4;

scanField["PhiQuotPlus", PhiQuot, 7/2, 1];

pp["Descendant Xp'", OPESimplify[v1, Expand]];
pp["Composite JXp", OPESimplify[v2, Expand]];
pp["Composite WYm", OPESimplify[v3, Expand]];
pp["Composite YpJ", OPESimplify[v4, Expand]];
