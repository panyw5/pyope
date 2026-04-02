#!/usr/bin/env wolframscript

Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

ClearAll[betaff, gammaff, bff, cff, T, J, G, Gbar, W, GW, Wbar, GWbar, NT, expanded, parts];

Bosonic[betaff, gammaff];
Fermionic[bff, cff];

OPE[bff, cff] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[betaff, gammaff] = MakeOPE[-One[w]/(z - w) + Ord[z, w, 0]];

J = 2 NO[bff, cff] + 3 NO[betaff, gammaff];
G = NO[gammaff, bff];
Gbar = 2 NO[Derivative[1][betaff], cff] + 3 NO[betaff, Derivative[1][cff]];
T = -2 NO[bff, Derivative[1][cff]] - 3/2 NO[betaff, Derivative[1][gammaff]] - NO[Derivative[1][bff], cff] - 1/2 NO[Derivative[1][betaff], gammaff];
W = betaff;
GW = bff;
Wbar = (
  NO[betaff, NO[betaff, NO[gammaff, NO[gammaff, gammaff]]]]
  + 2 NO[betaff, NO[gammaff, NO[gammaff, NO[bff, cff]]]]
  - 4 NO[betaff, NO[Derivative[1][gammaff], gammaff]]
  - 4/3 NO[gammaff, NO[bff, Derivative[1][cff]]]
  + 2/3 NO[gammaff, NO[Derivative[1][bff], cff]]
  + 2/3 NO[Derivative[1][betaff], NO[gammaff, gammaff]]
  - 8/3 NO[Derivative[1][gammaff], NO[bff, cff]]
  + 10/9 Derivative[2][gammaff]
);
GWbar = (
  -10/9 Derivative[3][cff]
  - 8/3 NO[bff, NO[Derivative[2][cff], cff]]
  - 3 NO[betaff, NO[betaff, NO[gammaff, NO[gammaff, Derivative[1][cff]]]]]
  + 4 NO[betaff, NO[gammaff, NO[bff, NO[Derivative[1][cff], cff]]]]
  + 4 NO[betaff, NO[gammaff, Derivative[2][cff]]]
  + 4 NO[betaff, NO[Derivative[1][gammaff], Derivative[1][cff]]]
  - 2/3 NO[Derivative[2][betaff], NO[gammaff, cff]]
  + 2/3 NO[Derivative[1][bff], NO[Derivative[1][cff], cff]]
  - 2 NO[Derivative[1][betaff], NO[betaff, NO[gammaff, NO[gammaff, cff]]]]
  + 8/3 NO[Derivative[1][betaff], NO[Derivative[1][gammaff], cff]]
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
Print["NullQ: ", expanded === 0];
Print["Length: ", StringLength[ToString[expanded, InputForm]]];
parts = If[expanded === 0, {}, If[Head[expanded] === Plus, List @@ expanded, {expanded}]];
Print["TermCount: ", Length[parts]];
If[expanded =!= 0,
  Print["FirstTerms:"];
  Do[Print[ToString[parts[[i]], InputForm]], {i, 1, Min[10, Length[parts]]}]
];
