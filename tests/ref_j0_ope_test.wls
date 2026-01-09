#!/usr/bin/env wolframscript
(* 验证 j0 的 OPE 计算 *)

<< OPEdefs`

(* 定义自由场 *)
Bosonic[b, c];
Fermionic[β, γ];

(* 定义基本 OPE *)
OPE[b, c] = MakeOPE[One (z - w)^-1 + Ord[z, w, 0]];
OPE[β, γ] = MakeOPE[-One (z - w)^-1 + Ord[z, w, 0]];

(* 定义 j0 *)
j0 = 2 NO[b, c] + 3 NO[β, γ];

Print["j0 = ", j0];
Print[];

(* 计算 OPE(j0, j0) *)
Print["计算 OPE(j0, j0):"];
ope_j0_j0 = OPE[j0, j0];
Print["OPE(j0, j0) = ", ope_j0_j0];
Print[];

(* 提取极点 *)
Print["极点分析:"];
pole2 = GetPole[ope_j0_j0, 2];
Print["pole(2) = ", pole2];
Print[];

pole1 = GetPole[ope_j0_j0, 1];
Print["pole(1) = ", pole1];
