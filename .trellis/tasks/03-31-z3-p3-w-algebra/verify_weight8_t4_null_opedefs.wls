#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

ClearAll[β, γ, b, c, T, J, G, Gbar, W, GW, Wbar, GWbar, NT, expanded];

Bosonic[β, γ];
Fermionic[b, c];

OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[β, γ] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

J = 2 NO[b, c] + 3 NO[β, γ];
G = NO[γ, b];
Gbar = 2 NO[Derivative[1][β], c] + 3 NO[β, Derivative[1][c]];
T = -2 NO[b, Derivative[1][c]] - 3/2 NO[β, Derivative[1][γ]] - NO[Derivative[1][b], c] - 1/2 NO[Derivative[1][β], γ];
W = β;
GW = b;
Wbar = (
  NO[β, NO[β, NO[γ, NO[γ, γ]]]]
  + 2 NO[β, NO[γ, NO[γ, NO[b, c]]]]
  - 4 NO[β, NO[Derivative[1][γ], γ]]
  - 4/3 NO[γ, NO[b, Derivative[1][c]]]
  + 2/3 NO[γ, NO[Derivative[1][b], c]]
  + 2/3 NO[Derivative[1][β], NO[γ, γ]]
  - 8/3 NO[Derivative[1][γ], NO[b, c]]
  + 10/9 Derivative[2][γ]
);
GWbar = (
  -10/9 Derivative[3][c]
  - 8/3 NO[b, NO[Derivative[2][c], c]]
  - 3 NO[β, NO[β, NO[γ, NO[γ, Derivative[1][c]]]]]
  + 4 NO[β, NO[γ, NO[b, NO[Derivative[1][c], c]]]]
  + 4 NO[β, NO[γ, Derivative[2][c]]]
  + 4 NO[β, NO[Derivative[1][γ], Derivative[1][c]]]
  - 2/3 NO[Derivative[2][β], NO[γ, c]]
  + 2/3 NO[Derivative[1][b], NO[Derivative[1][c], c]]
  - 2 NO[Derivative[1][β], NO[β, NO[γ, NO[γ, c]]]]
  + 8/3 NO[Derivative[1][β], NO[Derivative[1][γ], c]]
);

NT = NO[T, NO[T, NO[T, T]]]
  - 2/3 NO[T, NO[T, NO[T, NO[J, J]]]]
  + 1/3 NO[T, NO[T, NO[T, Derivative[1][J]]]]
  - 2/3 NO[T, NO[T, NO[J, NO[G, Gbar]]]]
  + 1/2 NO[T, NO[T, NO[G, Derivative[1][Gbar]]]]
  - 9/4 NO[T, NO[T, NO[W, Derivative[1][Wbar]]]]
  - 1/12 NO[T, NO[T, NO[Derivative[1][J], NO[J, J]]]]
  - 1/3 NO[T, NO[T, NO[Derivative[1][J], Derivative[1][J]]]]
  - 1/6 NO[T, NO[T, NO[Derivative[1][G], Gbar]]]
  + 15/4 NO[T, NO[T, NO[Derivative[1][W], Wbar]]]
  + 3/4 NO[T, NO[T, NO[Derivative[2][J], J]]];

expanded = Expand[NT];

Print["Expanded N_T:"];
Print[expanded];
Print[""];
Print["NullQ: ", expanded === 0];
Print["InputForm:"];
Print[ToString[expanded, InputForm]];
