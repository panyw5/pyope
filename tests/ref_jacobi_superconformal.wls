#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
 * Jacobi 恒等式参考测试 - N=1 超共形代数
 *
 * 本脚本使用 OPEdefs.m 计算 N=1 超共形代数的 Jacobi 恒等式，
 * 作为 pyope 实现的参考结果。
 *
 * N=1 超共形代数包含：
 * - T(z): 能动张量，h=2，玻色子（parity=0）
 * - G(z): 超流，h=3/2，费米子（parity=1）
 *
 * OPE 公式：
 * 1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)
 * 2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)
 * 3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)
 *
 * 测试内容：
 * 1. OPEJacobi[G, G, G] - 最关键，测试费米子符号
 * 2. OPEJacobi[T, G, G] - 混合玻色-费米
 * 3. OPEJacobi[T, T, G] - 验证 G 的变换性质
 *
 * 参考：
 * - OPEdefs.m 第 1601-1637 行：OPEJacobi 实现
 * - OPEdefs.m 第 740-743 行：Fermionic 声明
 * - OPEdefs.m 第 1604 行：费米子符号因子 sign = (-1)^(OPEParity[A] OPEParity[B])
 *)

Print["========================================"];
Print["N=1 Superconformal Algebra Jacobi Test"];
Print["========================================"];
Print[""];

(* 加载 OPEdefs.m *)
Print["Loading OPEdefs.m..."];
SetDirectory[FileNameJoin[{Directory[], "..", "OPEdefs"}]];
<< OPEdefs.m
Print["OPEdefs.m loaded successfully!"];
Print[""];

(* 定义符号 *)
Print["Defining symbols..."];
c = Symbol["c"];  (* 中心荷 *)
Print["  c: central charge"];
Print[""];

(* 定义算符 *)
Print["Defining operators..."];
Print["  T: Energy-momentum tensor (bosonic, h=2)"];
Bosonic[T];
Print["    OPEParity[T] = ", OPEParity[T]];

Print["  G: Supercurrent (fermionic, h=3/2)"];
Fermionic[G];
Print["    OPEParity[G] = ", OPEParity[G]];
Print[""];

(* 定义 OPE *)
Print["Defining OPEs..."];

Print["  1. T(z)T(w) = c/2/(z-w)^4 + 2T(w)/(z-w)^2 + ∂T(w)/(z-w)"];
OPE[T, T] = OPEData[{c/2 One, 0, 2 T, Derivative[1][T]}];

Print["  2. T(z)G(w) = (3/2)G(w)/(z-w)^2 + ∂G(w)/(z-w)"];
OPE[T, G] = OPEData[{0, 3/2 G, Derivative[1][G]}];

Print["  3. G(z)G(w) = (2c/3)/(z-w)^3 + 2T(w)/(z-w)"];
OPE[G, G] = OPEData[{2*c/3 One, 0, 2 T}];

Print[""];

(* 测试 Jacobi 恒等式 *)
Print["========================================"];
Print["Testing Jacobi Identities"];
Print["========================================"];
Print[""];

(* 测试 1: OPEJacobi[G, G, G] *)
Print["Test 1: OPEJacobi[G, G, G]"];
Print["  This is the most critical test for fermionic sign factors"];
Print["  Sign factor: (-1)^(OPEParity[G] * OPEParity[G]) = (-1)^(1*1) = -1"];
Print[""];
Print["  Computing..."];
jacobiGGG = OPEJacobi[G, G, G];
Print["  Result dimensions: ", Dimensions[jacobiGGG]];
Print["  Result:"];
Print[MatrixForm[jacobiGGG]];

(* 检查是否全为 0 *)
allZeroGGG = AllTrue[Flatten[jacobiGGG], # === 0 &];
Print["  All entries zero? ", allZeroGGG];
If[allZeroGGG,
    Print["  SUCCESS: Jacobi identity holds for G-G-G!"],
    Print["  FAILURE: Some entries are non-zero!"]
];
Print[""];

(* 测试 2: OPEJacobi[T, G, G] *)
Print["Test 2: OPEJacobi[T, G, G]"];
Print["  Mixed bosonic-fermionic test"];
Print["  Sign factor: (-1)^(OPEParity[T] * OPEParity[G]) = (-1)^(0*1) = 1"];
Print[""];
Print["  Computing..."];
jacobiTGG = OPEJacobi[T, G, G];
Print["  Result dimensions: ", Dimensions[jacobiTGG]];
Print["  Result:"];
Print[MatrixForm[jacobiTGG]];

(* 检查是否全为 0 *)
allZeroTGG = AllTrue[Flatten[jacobiTGG], # === 0 &];
Print["  All entries zero? ", allZeroTGG];
If[allZeroTGG,
    Print["  SUCCESS: Jacobi identity holds for T-G-G!"],
    Print["  FAILURE: Some entries are non-zero!"]
];
Print[""];

(* 测试 3: OPEJacobi[T, T, G] *)
Print["Test 3: OPEJacobi[T, T, G]"];
Print["  Verifies transformation properties of G"];
Print["  Sign factor: (-1)^(OPEParity[T] * OPEParity[T]) = (-1)^(0*0) = 1"];
Print[""];
Print["  Computing..."];
jacobiTTG = OPEJacobi[T, T, G];
Print["  Result dimensions: ", Dimensions[jacobiTTG]];
Print["  Result:"];
Print[MatrixForm[jacobiTTG]];

(* 检查是否全为 0 *)
allZeroTTG = AllTrue[Flatten[jacobiTTG], # === 0 &];
Print["  All entries zero? ", allZeroTTG];
If[allZeroTTG,
    Print["  SUCCESS: Jacobi identity holds for T-T-G!"],
    Print["  FAILURE: Some entries are non-zero!"]
];
Print[""];

(* 总结 *)
Print["========================================"];
Print["Summary"];
Print["========================================"];
Print["  OPEJacobi[G, G, G]: ", If[allZeroGGG, "PASS", "FAIL"]];
Print["  OPEJacobi[T, G, G]: ", If[allZeroTGG, "PASS", "FAIL"]];
Print["  OPEJacobi[T, T, G]: ", If[allZeroTTG, "PASS", "FAIL"]];
Print[""];

allPass = allZeroGGG && allZeroTGG && allZeroTTG;
If[allPass,
    Print["ALL TESTS PASSED!"],
    Print["SOME TESTS FAILED!"]
];
Print[""];

(* 导出结果供 Python 测试使用 *)
Print["Exporting results for Python comparison..."];
results = Association[
    "jacobi_GGG" -> jacobiGGG,
    "jacobi_TGG" -> jacobiTGG,
    "jacobi_TTG" -> jacobiTTG,
    "all_zero_GGG" -> allZeroGGG,
    "all_zero_TGG" -> allZeroTGG,
    "all_zero_TTG" -> allZeroTTG,
    "dimensions_GGG" -> Dimensions[jacobiGGG],
    "dimensions_TGG" -> Dimensions[jacobiTGG],
    "dimensions_TTG" -> Dimensions[jacobiTTG]
];

(* 将结果转换为可导出的格式 *)
exportData = Association[
    "jacobi_GGG" -> Map[ToString, jacobiGGG, {2}],
    "jacobi_TGG" -> Map[ToString, jacobiTGG, {2}],
    "jacobi_TTG" -> Map[ToString, jacobiTTG, {2}],
    "all_zero_GGG" -> allZeroGGG,
    "all_zero_TGG" -> allZeroTGG,
    "all_zero_TTG" -> allZeroTTG,
    "dimensions_GGG" -> Dimensions[jacobiGGG],
    "dimensions_TGG" -> Dimensions[jacobiTGG],
    "dimensions_TTG" -> Dimensions[jacobiTTG]
];

outputFile = FileNameJoin[{Directory[], "..", "tests", "ref_jacobi_superconformal_output.json"}];
Export[outputFile, exportData, "JSON"];
Print["Results exported to: ", outputFile];
Print[""];

Print["========================================"];
Print["Reference test completed!"];
Print["========================================"];
