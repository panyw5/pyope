#!/usr/bin/env wolframscript
(* 
 * Wolfram Script to test nested NO expansion in OPEdefs.m
 * This script generates reference results for pyope's expand_nested_no function
 *)

(* Load OPEdefs.m - use absolute path *)
Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

(* Helper function to format output *)
FormatResult[label_, expr_, result_] := Module[{},
    Print["========================================"];
    Print["Test: ", label];
    Print["Input: ", expr];
    Print["Output: ", result];
    Print[""];
];

(* ============================================ *)
(* Test 1: sl(2) Kac-Moody algebra *)
(* ============================================ *)
Print[""];
Print["######## SL(2) KAC-MOODY ALGEBRA ########"];
Print[""];

(* Define sl(2) generators using Declare (matching your Mathematica code) *)
Declare[{}][{Jplus, Jzero, Jminus}][{}];

(* Define OPEs for sl(2) at level k (matching your Mathematica code) *)
k = Symbol["k"];  (* level parameter *)

(* Use MakeOPE format to match your code *)
OPE[Jplus, Jzero] = MakeOPE[-((2 Jplus[w])/(z - w)) + Ord[z, w, 0]];
OPE[Jzero, Jzero] = MakeOPE[(2 k One[w])/(z - w)^2 + Ord[z, w, 0]];
OPE[Jzero, Jminus] = MakeOPE[(-2 Jminus[w])/(z - w) + Ord[z, w, 0]];
OPE[Jplus, Jminus] = MakeOPE[(k One[w])/(z - w)^2 + Jzero[w]/(z - w) + Ord[z, w, 0]];

(* Test Case 1.1: NO(NO(Jplus, Jzero), Jminus) - left nested *)
expr1 = NO[NO[Jplus, Jzero], Jminus];
result1 = Expand[expr1];
FormatResult["SL2-1: NO(NO(J+, J0), J-)", expr1, result1];

(* Test Case 1.2: NO(Jplus, NO(Jzero, Jminus)) - right nested *)
expr2 = NO[Jplus, NO[Jzero, Jminus]];
result2 = Expand[expr2];
FormatResult["SL2-2: NO(J+, NO(J0, J-))", expr2, result2];

(* Test Case 1.3: NO(NO(Jminus, Jplus), Jzero) - left nested with non-trivial OPE *)
expr3 = NO[NO[Jminus, Jplus], Jzero];
result3 = Expand[expr3];
FormatResult["SL2-3: NO(NO(J-, J+), J0)", expr3, result3];

(* Test Case 1.4: NO(Jminus, NO(Jplus, Jzero)) - right nested (matching your test) *)
expr4 = NO[Jminus, NO[Jplus, Jzero]];
result4 = Expand[expr4];
FormatResult["SL2-4: NO(J-, NO(J+, J0))", expr4, result4];

(* Test Case 1.5: Triple nested - NO(NO(Jminus, NO(Jplus, Jzero))) *)
expr5 = NO[NO[Jminus, NO[Jplus, Jzero]]];
result5 = Expand[expr5];
FormatResult["SL2-5: NO(NO(J-, NO(J+, J0)))", expr5, result5];

(* Test Case 1.6: NO(NO(Jzero, NO(Jplus, Jminus))) - matching your second test *)
expr6 = NO[NO[Jzero, NO[Jplus, Jminus]]];
result6 = Expand[expr6];
FormatResult["SL2-6: NO(NO(J0, NO(J+, J-)))", expr6, result6];

(* ============================================ *)
(* Test 2: Virasoro algebra *)
(* ============================================ *)
Print[""];
Print["######## VIRASORO ALGEBRA ########"];
Print[""];

(* Define Virasoro generator *)
Declare[{}][{T}][{}];

(* Define central charge *)
c = Symbol["c"];

(* Define Virasoro OPE: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w) *)
OPE[T, T] = MakeOPE[(c One[w])/(2(z - w)^4) + (2 T[w])/(z - w)^2 + Derivative[1][T][w]/(z - w) + Ord[z, w, 0]];

(* Test Case 2.1: NO(NO(T, T), T) - left nested *)
expr7 = NO[NO[T, T], T];
result7 = Expand[expr7];
FormatResult["VIR-1: NO(NO(T, T), T)", expr7, result7];

(* Test Case 2.2: NO(T, NO(T, T)) - right nested *)
expr8 = NO[T, NO[T, T]];
result8 = Expand[expr8];
FormatResult["VIR-2: NO(T, NO(T, T))", expr8, result8];

(* Test Case 2.3: NO(NO(T, Derivative[T]), T) - with derivative *)
expr9 = NO[NO[T, Derivative[1][T]], T];
result9 = Expand[expr9];
FormatResult["VIR-3: NO(NO(T, dT), T)", expr9, result9];

(* Test Case 2.4: NO(T, NO(Derivative[T], T)) - derivative in nested *)
expr10 = NO[T, NO[Derivative[1][T], T]];
result10 = Expand[expr10];
FormatResult["VIR-4: NO(T, NO(dT, T))", expr10, result10];

(* Test Case 2.5: NO(NO(Derivative[T], T), Derivative[T]) - both derivatives *)
expr11 = NO[NO[Derivative[1][T], T], Derivative[1][T]];
result11 = Expand[expr11];
FormatResult["VIR-5: NO(NO(dT, T), dT)", expr11, result11];

(* ============================================ *)
(* Test 3: Mixed sl(2) and Virasoro *)
(* ============================================ *)
Print[""];
Print["######## MIXED SL(2) AND VIRASORO ########"];
Print[""];

(* Assume T and J+ commute (no OPE defined) *)
(* Test Case 3.1: NO(NO(T, Jplus), Jminus) *)
expr12 = NO[NO[T, Jplus], Jminus];
result12 = Expand[expr12];
FormatResult["MIX-1: NO(NO(T, J+), J-)", expr12, result12];

(* Test Case 3.2: NO(T, NO(Jplus, Jminus)) *)
expr13 = NO[T, NO[Jplus, Jminus]];
result13 = Expand[expr13];
FormatResult["MIX-2: NO(T, NO(J+, J-))", expr13, result13];

Print[""];
Print["######## ALL TESTS COMPLETED ########"];
Print[""];
