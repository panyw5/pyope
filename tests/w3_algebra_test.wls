#!/usr/bin/env wolframscript
(* W3 代数测试 - Mathematica 参考实现 *)

(* 设置工作目录并导入 OPEdefs *)
SetDirectory[DirectoryName[$InputFileName]];
SetDirectory[ParentDirectory[]];
SetDirectory["OPEdefs"];
Import["OPEdefs.m"];

(* ============================================================================ *)
(* W3 代数定义 *)
(* ============================================================================ *)

Print["[W3 Algebra Test - Mathematica Reference]"];
Print["=========================================="];
Print[];

(* 声明算符类型 *)
Bosonic[T, W];

(* 定义辅助算符 Λ *)
Λ[w_] := NO[T, T][w] - (3/10) * T''[w];

(* 定义 OPE 规则 *)
Print["Defining W3 algebra OPEs..."];

(* T-T OPE *)
OPE[T, T] := MakeOPE[
    c * One[w] / (2 * (z - w)^4) +
    (2 * T[w]) / (z - w)^2 +
    T'[w] / (z - w) +
    Ord[z, w, 0]
];

(* T-W OPE *)
OPE[T, W] := MakeOPE[
    (3 * W[w]) / (z - w)^2 +
    W'[w] / (z - w) +
    Ord[z, w, 0]
];

(* W-W OPE *)
OPE[W, W] := MakeOPE[
    (c * One[w]) / (z - w)^6 +
    (2 * T[w]) / (z - w)^4 +
    T'[w] / (z - w)^3 +
    (1 / (z - w)^2) * (2 * β * Λ[w] + (3/10) * T''[w]) +
    (1 / (z - w)) * (β * Λ'[w] + (1/15) * T'''[w]) +
    Ord[z, w, 0]
];

Print["OPE definitions completed"];
Print[];

(* ============================================================================ *)
(* Test 1: 基本 OPE 计算 *)
(* ============================================================================ *)

Print["Test 1: T-T OPE"];
Print["================"];
Print["OPE[T, T] = ", OPE[T, T]];
Print[];

Print["Test 2: T-W OPE"];
Print["================"];
Print["OPE[T, W] = ", OPE[T, W]];
Print[];

Print["Test 3: W-W OPE"];
Print["================"];
Print["OPE[W, W] = ", OPE[W, W]];
Print[];

(* ============================================================================ *)
(* Test 4: 导数算符 OPE *)
(* ============================================================================ *)

Print["Test 4: OPE[T, NO[T,T] - (3/10)T'']"];
Print["===================================="];
Print["Result = ", OPE[T, NO[T, T] - (3/10) * T'']];
Print[];

Print["Test 5: OPE[T, W'']"];
Print["==================="];
Print["Result = ", OPE[T, W'']];

(* ============================================================================ *)
(* Test 6: 复合正规序 OPE *)
(* ============================================================================ *)

Print["Test 6: NO[T', NO[W'', T]]"];
Print["==========================="];
Print["Result = ", Expand[NO[T', NO[W'', T]]]];
Print[];

Print["Test 7: NO[W', NO[W'', T]]"];
Print["==========================="];
Print["Result = ", Expand[NO[W', NO[W'', T]]]];
Print[];

(* ============================================================================ *)
(* Test 8: 数值检验（设定 c=100, β=1/10） *)
(* ============================================================================ *)

Print["Test 8: Numerical verification"];
Print["================"];
Print["Setting c=100, β=1/10"];

(* 替换符号参数 *)
Print["OPE[W,W] numerical result = ", OPE[W, W] /. {c -> 100, β -> 1/10}];
Print[];

(* ============================================================================ *)
(* 输出结果摘要 *)
(* ============================================================================ *)

Print["Tests completed!"];
Print["All test results have been output and can be used to compare with Python implementation"];
