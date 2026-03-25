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

p1 = NO[beta''', beta, gamma];
p2 = NO[beta'', beta, gamma'];
p3 = NO[beta'', beta', gamma];
p4 = NO[beta', beta, gamma''];
p5 = NO[beta', beta', gamma'];
p6 = NO[beta'', beta', gamma];
p7 = NO[beta, beta, gamma'''];
p8 = NO[beta', beta, gamma''];
p9 = NO[beta'', beta, gamma'];
p10 = NO[beta''', beta, gamma];
p11 = NO[beta', beta, beta, gamma, gamma];
p12 = NO[beta, beta, beta, gamma', gamma];
p13 = NO[beta, beta, b', c'];
p14 = NO[beta, beta, b, c''];
p15 = NO[beta', beta, b, c'];
p16 = NO[beta, gamma, gamma, b, c'];

n1 = NO[gamma''', beta, gamma];
n2 = NO[gamma'', beta', gamma];
n3 = NO[gamma'', beta, gamma'];
n4 = NO[gamma', beta'', gamma];
n5 = NO[gamma', beta', gamma'];
n6 = NO[gamma', beta, gamma''];
n7 = NO[gamma, beta''', gamma];
n8 = NO[gamma, beta'', gamma'];
n9 = NO[gamma, beta', gamma''];
n10 = NO[gamma, beta, gamma'''];
n11 = NO[gamma', gamma, gamma, beta, beta];
n12 = NO[gamma, gamma, gamma, beta', beta];
n13 = NO[gamma, gamma, b', c'];
n14 = NO[gamma, gamma, b, c''];
n15 = NO[gamma', gamma, b, c'];
n16 = NO[gamma, beta, beta, b, c'];

PhiP = ap1*p1 + ap2*p2 + ap3*p3 + ap4*p4 + ap5*p5 + ap6*p6 + ap7*p7 + ap8*p8 + ap9*p9 + ap10*p10 + ap11*p11 + ap12*p12 + ap13*p13 + ap14*p14 + ap15*p15 + ap16*p16;
PhiM = am1*n1 + am2*n2 + am3*n3 + am4*n4 + am5*n5 + am6*n6 + am7*n7 + am8*n8 + am9*n9 + am10*n10 + am11*n11 + am12*n12 + am13*n13 + am14*n14 + am15*n15 + am16*n16;

pp["Basis (h=7/2,q=+1)", {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16}];
pp["Basis (h=7/2,q=-1)", {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16}];

scanField["PhiP", PhiP, 7/2, 1];
scanField["PhiM", PhiM, 7/2, -1];

pp["Descendant dXp", Xp'];
pp["Descendant dXm", Xm'];
pp["Composite JXp", NO[J, Xp]];
pp["Composite JXm", NO[J, Xm]];
