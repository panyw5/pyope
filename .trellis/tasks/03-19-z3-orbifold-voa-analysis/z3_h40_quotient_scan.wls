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
  pp[name <> " J pole1", OPEPole[1][oj]];
  showPoleData[name <> " T", ot];
  pp[name <> " T pole2 - 4 phi", OPESimplify[OPEPole[2][ot] - 4 phi, Expand]];
  pp[name <> " T pole1 - d phi", OPESimplify[OPEPole[1][ot] - phi', Expand]];
  showPoleData[name <> " G", og];
  showPoleData[name <> " Gb", ogb];
];

u1 = NO[beta''', gamma];
u2 = NO[beta'', gamma'];
u3 = NO[beta', gamma''];
u4 = NO[beta, gamma'''];
u5 = NO[beta'', beta, gamma, gamma];
u6 = NO[beta', beta', gamma, gamma];
u7 = NO[beta', beta, gamma', gamma];
u8 = NO[beta, beta, gamma'', gamma];
u9 = NO[beta, beta, gamma', gamma'];
u10 = NO[b'', c'];
u11 = NO[b', c''];
u12 = NO[b, c'''];
u13 = NO[beta', gamma, b, c'];
u14 = NO[beta, gamma', b, c'];
u15 = NO[beta, gamma, b', c'];
u16 = NO[beta, gamma, b, c''];
u17 = NO[beta, beta, gamma, gamma, b, c'];

v1 = J''';
v2 = T'';
v3 = NO[J, J]';
v4 = NO[J, T];
v5 = NO[T, T];
v6 = OPEPole[1][OPESimplify[OPE[W, Wb], Expand]]';

pp["Raw basis h4 neutral", {u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17}];
pp["Quotient generators", {v1,v2,v3,v4,v5,v6}];

PhiRaw = a1*u1 + a2*u2 + a3*u3 + a4*u4 + a5*u5 + a6*u6 + a7*u7 + a8*u8 + a9*u9 + a10*u10 + a11*u11 + a12*u12 + a13*u13 + a14*u14 + a15*u15 + a16*u16 + a17*u17;
PhiQ = PhiRaw - q1*v1 - q2*v2 - q3*v3 - q4*v4 - q5*v5 - q6*v6;

scanField["PhiQ_h4_neutral", PhiQ, 4, 0];
