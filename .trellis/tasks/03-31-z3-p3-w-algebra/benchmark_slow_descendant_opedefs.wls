#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

ClearAll[beta, gamma, b, c, T, J, G, Gbar, W, GW, Wbar, GWbar, expr, expanded, timing];

Bosonic[beta, gamma];
Fermionic[b, c];
OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[beta, gamma] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

T = -2 NO[b, Derivative[1][c]] - 3/2 NO[beta, Derivative[1][gamma]] - NO[Derivative[1][b], c] - 1/2 NO[Derivative[1][beta], gamma];
J = 2 NO[b, c] + 3 NO[beta, gamma];
G = NO[gamma, b];
Gbar = 2 NO[Derivative[1][beta], c] + 3 NO[beta, Derivative[1][c]];
W = beta;
GW = b;
Wbar = (
  NO[beta, NO[beta, NO[gamma, NO[gamma, gamma]]]]
  + 2 NO[beta, NO[gamma, NO[gamma, NO[b, c]]]]
  - 4 NO[beta, NO[Derivative[1][gamma], gamma]]
  - 4/3 NO[gamma, NO[b, Derivative[1][c]]]
  + 2/3 NO[gamma, NO[Derivative[1][b], c]]
  + 2/3 NO[Derivative[1][beta], NO[gamma, gamma]]
  - 8/3 NO[Derivative[1][gamma], NO[b, c]]
  + 10/9 Derivative[2][gamma]
);
GWbar = (
  -10/9 Derivative[3][c]
  - 8/3 NO[b, NO[Derivative[2][c], c]]
  - 3 NO[beta, NO[beta, NO[gamma, NO[gamma, Derivative[1][c]]]]]
  + 4 NO[beta, NO[gamma, NO[b, NO[Derivative[1][c], c]]]]
  + 4 NO[beta, NO[gamma, Derivative[2][c]]]
  + 4 NO[beta, NO[Derivative[1][gamma], Derivative[1][c]]]
  - 2/3 NO[Derivative[2][beta], NO[gamma, c]]
  + 2/3 NO[Derivative[1][b], NO[Derivative[1][c], c]]
  - 2 NO[Derivative[1][beta], NO[beta, NO[gamma, NO[gamma, c]]]]
  + 8/3 NO[Derivative[1][beta], NO[Derivative[1][gamma], c]]
);

expr = 4 NO[T, NO[J, NO[J, NO[W, Derivative[1][Wbar]]]]]/9
  + 4 NO[T, NO[J, NO[J, NO[Derivative[1][W], Wbar]]]]/9
  + 8 NO[T, NO[Derivative[1][J], NO[J, NO[W, Wbar]]]]/9
  - NO[J, NO[J, NO[J, NO[J, NO[W, Derivative[1][Wbar]]]]]]/27
  - NO[J, NO[J, NO[J, NO[J, NO[Derivative[1][W], Wbar]]]]]/27
  + 2 NO[J, NO[G, NO[Gbar, NO[W, Derivative[1][Wbar]]]]]/9
  + 2 NO[J, NO[G, NO[Gbar, NO[Derivative[1][W], Wbar]]]]/9
  + 2 NO[J, NO[G, NO[Derivative[1][Gbar], NO[W, Wbar]]]]/9
  + 2 NO[J, NO[W, NO[W, NO[Derivative[1][Wbar], Wbar]]]]
  + 2 NO[J, NO[Derivative[1][G], NO[Gbar, NO[W, Wbar]]]]/9
  + 2 NO[J, NO[Derivative[1][W], NO[W, NO[Wbar, Wbar]]]]
  + 4 NO[Derivative[1][T], NO[J, NO[J, NO[W, Wbar]]]]/9
  - 4 NO[Derivative[1][J], NO[J, NO[J, NO[J, NO[W, Wbar]]]]]/27
  + 2 NO[Derivative[1][J], NO[G, NO[Gbar, NO[W, Wbar]]]]/9
  + NO[Derivative[1][J], NO[W, NO[W, NO[Wbar, Wbar]]]];

timing = AbsoluteTiming[expanded = Expand[expr];][[1]];
Print[ToString[<|"engine" -> "OPEdefs", "seconds" -> timing, "zeroQ" -> (expanded === 0), "length" -> StringLength[ToString[expanded, InputForm]]|>, InputForm]];
