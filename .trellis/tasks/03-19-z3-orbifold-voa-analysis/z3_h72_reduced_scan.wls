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

u2 = NO[beta'', NO[beta, gamma']];
u3 = NO[beta'', NO[beta', gamma]];
u4 = NO[beta', NO[beta, gamma'']];
u5 = NO[beta', NO[beta', gamma']];
u6 = NO[beta', NO[beta, NO[beta, NO[gamma, gamma]]]];
u7 = NO[beta, NO[beta, NO[beta, NO[gamma', gamma]]]];
u8 = NO[beta, NO[beta, NO[b', c']]];
u9 = NO[beta, NO[beta, NO[b, c'']]];
u10 = NO[beta', NO[beta, NO[b, c']]];
u11 = NO[beta, NO[gamma, NO[gamma, NO[b, c']]]];

(* reduced by previous necessary conditions: r1=0, r2=-r4, r3=-r4-r5, q3=q4=0 *)
PhiRed = rr4*(u4 - u2 - u3) + rr5*(u5 - u3) + rr6*u6 + rr7*u7 + rr8*u8 + rr9*u9 + rr10*u10 + rr11*u11 - qq1*(NO[beta, NO[beta, b']] + 2*NO[beta', NO[beta, b]]) - qq2*(NO[beta, NO[beta, NO[beta, NO[gamma, b]]]] + 2*NO[beta', NO[beta, b]]);

scanField["PhiRedPlus", PhiRed, 7/2, 1];
