#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
Z3 W-代数 OPE 参考计算
使用 OPEdefs.m 计算 Z3 W-代数生成元之间的 OPE
*)

(* 加载 OPEdefs 包 *)
<< OPEdefs`

(* 定义自由场 *)
Bosonic[\[Beta]1, \[Gamma]1];
Fermionic[b1, c1];

(* 定义自由场 OPE *)
OPE[\[Beta]1, \[Gamma]1] = MakeOPE[-One (z-w)^-1 + Ord[z, w, 0]];
OPE[b1, c1] = MakeOPE[One (z-w)^-1 + Ord[z, w, 0]];

(* 定义 Z3 W-代数生成元 *)
w = \[Beta]1;
j0 = 2 NO[b1, c1] + 3 NO[\[Beta]1, \[Gamma]1];
t = -2 NO[b1, c1'] - 3/2 NO[\[Beta]1, \[Gamma]1'] - NO[b1', c1] - 1/2 NO[\[Beta]1', \[Gamma]1];
g = NO[\[Gamma]1, b1];
gt = 2 NO[\[Beta]1', c1] + 3 NO[\[Beta]1, c1'];
gw = b1;

(* 计算关键 OPE *)
Print["计算 Z3 W-代数生成元的 OPE..."];
Print[""];

(* 1. j0(z)j0(w) *)
Print["1. j0(z)j0(w):"];
opeJ0J0 = OPE[j0, j0];
Print["  OPE = ", opeJ0J0];
Print["  pole(2) = ", GetCoefficient[opeJ0J0, 2]];
Print[""];

(* 2. t(z)t(w) *)
Print["2. t(z)t(w):"];
opeTT = OPE[t, t];
Print["  OPE = ", opeTT];
Print["  pole(4) = ", GetCoefficient[opeTT, 4]];
Print["  pole(2) = ", GetCoefficient[opeTT, 2]];
Print[""];

(* 3. t(z)j0(w) *)
Print["3. t(z)j0(w):"];
opeTJ0 = OPE[t, j0];
Print["  OPE = ", opeTJ0];
Print["  pole(2) = ", GetCoefficient[opeTJ0, 2]];
Print["  pole(1) = ", GetCoefficient[opeTJ0, 1]];
Print[""];

(* 4. t(z)w(w) *)
Print["4. t(z)w(w):"];
opeTW = OPE[t, w];
Print["  OPE = ", opeTW];
Print["  pole(2) = ", GetCoefficient[opeTW, 2]];
Print["  pole(1) = ", GetCoefficient[opeTW, 1]];
Print[""];

(* 5. t(z)g(w) *)
Print["5. t(z)g(w):"];
opeTG = OPE[t, g];
Print["  OPE = ", opeTG];
Print["  pole(2) = ", GetCoefficient[opeTG, 2]];
Print["  pole(1) = ", GetCoefficient[opeTG, 1]];
Print[""];

(* 导出结果为 JSON *)
results = Association[
  "j0_j0" -> Association[
    "pole2" -> ToString[GetCoefficient[opeJ0J0, 2], InputForm]
  ],
  "t_t" -> Association[
    "pole4" -> ToString[GetCoefficient[opeTT, 4], InputForm],
    "pole2" -> ToString[GetCoefficient[opeTT, 2], InputForm]
  ],
  "t_j0" -> Association[
    "pole2" -> ToString[GetCoefficient[opeTJ0, 2], InputForm],
    "pole1" -> ToString[GetCoefficient[opeTJ0, 1], InputForm]
  ],
  "t_w" -> Association[
    "pole2" -> ToString[GetCoefficient[opeTW, 2], InputForm],
    "pole1" -> ToString[GetCoefficient[opeTW, 1], InputForm]
  ],
  "t_g" -> Association[
    "pole2" -> ToString[GetCoefficient[opeTG, 2], InputForm],
    "pole1" -> ToString[GetCoefficient[opeTG, 1], InputForm]
  ]
];

Export["tests/ref_z3_w_algebra_opes_output.json", results, "JSON"];
Print["结果已导出到 tests/ref_z3_w_algebra_opes_output.json"];
