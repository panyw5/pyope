#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
Jacobi 恒等式测试 - Virasoro 代数参考实现

本脚本使用 OPEdefs.m 计算 Virasoro 算符 T 的 Jacobi 恒等式，
作为 pyope 实现的参考结果。

测试内容：
1. 定义 Virasoro 算符 T，共形权重为 2
2. 使用 OPEJacobi[T, T, T] 检查 Jacobi 恒等式
3. 验证结果全为 0（Jacobi 恒等式成立）

参考：
- OPEdefs.m 第 1598-1637 行：OPEJacobi 实现
- OPEconf.m 第 417-432 行：Virasoro 定义
*)

(* 加载 OPEdefs 包 *)
Print["Loading OPEdefs.m..."];
SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{ParentDirectory[], "OPEdefs", "OPEdefs.m"}]];

(* 定义 Virasoro 算符 T *)
Print["Defining Virasoro operator T with central charge c..."];
c = Symbol["c"];  (* 中心荷作为符号参数 *)
Bosonic[T];
OPE[T, T] = MakeOPE[{c/2 One, 0, 2 T, T'}];

Print["OPE[T,T] = ", OPE[T,T]];
Print[""];

(* 验证 T 的 OPE 结构 *)
Print["Verifying T(z)T(w) OPE structure:"];
Print["  Pole 4: ", OPEPole[4][T,T], " (should be c/2*One)"];
Print["  Pole 3: ", OPEPole[3][T,T], " (should be 0)"];
Print["  Pole 2: ", OPEPole[2][T,T], " (should be 2*T)"];
Print["  Pole 1: ", OPEPole[1][T,T], " (should be T')"];
Print[""];

(* 计算 Jacobi 恒等式 OPEJacobi[T, T, T] *)
Print["Computing Jacobi identity OPEJacobi[T, T, T]..."];
Print["This checks: [T, [T, T]_q]_m - [T, [T, T]_m]_q - Sum[...] = 0"];
Print[""];

jacobiResult = OPEJacobi[T, T, T];

Print["Jacobi identity result:"];
Print["  Type: ", Head[jacobiResult]];
Print["  Dimensions: ", Dimensions[jacobiResult]];
Print["  Full result: ", jacobiResult];
Print[""];

(* 检查结果是否全为 0 *)
Print["Checking if all entries are zero..."];
allZero = True;

(* OPEJacobi 返回一个列表的列表 *)
If[Head[jacobiResult] === List,
  Do[
    If[Head[jacobiResult[[i]]] === List,
      (* 二维数组 *)
      Do[
        entry = jacobiResult[[i, j]];
        If[entry =!= 0,
          Print["  Non-zero entry at [", i, ",", j, "]: ", entry];
          allZero = False;
        ],
        {j, Length[jacobiResult[[i]]]}
      ],
      (* 一维数组 *)
      entry = jacobiResult[[i]];
      If[entry =!= 0,
        Print["  Non-zero entry at [", i, "]: ", entry];
        allZero = False;
      ]
    ],
    {i, Length[jacobiResult]}
  ],
  (* 不是列表 *)
  If[jacobiResult =!= 0,
    Print["  Result is non-zero: ", jacobiResult];
    allZero = False;
  ]
];

If[allZero,
  Print["SUCCESS: All entries are zero. Jacobi identity holds!"],
  Print["FAILURE: Some entries are non-zero. Jacobi identity violated!"]
];
Print[""];

(* 简化输出 *)
Print["Summary:"];
Print["  Test: Jacobi identity for Virasoro algebra"];
Print["  Operator: T (conformal weight 2, bosonic)"];
Print["  Central charge: c (symbolic)"];
Print["  Result: ", If[allZero, "PASS", "FAIL"]];
Print[""];
Print["Test completed!"];
