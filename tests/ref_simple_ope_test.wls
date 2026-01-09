#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
简单测试：验证基本的正规序 OPE
*)

(* 加载 OPEdefs 包 *)
<< OPEdefs`

(* 定义自由场 *)
Fermionic[b1, c1];
Bosonic[\[Beta]1, \[Gamma]1];

(* 定义自由场 OPE *)
OPE[b1, c1] = MakeOPE[One (z-w)^-1 + Ord[z, w, 0]];
OPE[\[Beta]1, \[Gamma]1] = MakeOPE[-One (z-w)^-1 + Ord[z, w, 0]];

Print["=== 基本 OPE 测试 ==="];
Print[""];

(* 测试 1: :bc: 与 :bc: 的 OPE *)
Print["1. OPE[:bc:, :bc:]:"];
bc = NO[b1, c1];
opeBCBC = OPE[bc, bc];
Print["  结果: ", opeBCBC];
Print[""];

(* 测试 2: :βγ: 与 :βγ: 的 OPE *)
Print["2. OPE[:βγ:, :βγ:]:"];
bg = NO[\[Beta]1, \[Gamma]1];
opeBGBG = OPE[bg, bg];
Print["  结果: ", opeBGBG];
Print[""];

(* 测试 3: (2:bc: + 3:βγ:) 与 (2:bc: + 3:βγ:) 的 OPE *)
Print["3. OPE[2:bc: + 3:βγ:, 2:bc: + 3:βγ:]:"];
j0 = 2*bc + 3*bg;
opeJ0J0 = OPE[j0, j0];
Print["  结果: ", opeJ0J0];
Print[""];

(* 测试 4: :b∂c: 与 :b∂c: 的 OPE *)
Print["4. OPE[:b∂c:, :b∂c:]:"];
bdc = NO[b1, c1'];
opeBDCBDC = OPE[bdc, bdc];
Print["  结果: ", opeBDCBDC];
Print[""];
