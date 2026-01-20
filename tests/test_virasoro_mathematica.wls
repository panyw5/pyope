#!/usr/bin/env wolframscript
(* 
 * Virasoro Algebra Nested NO Expansion Tests
 * Reference results for pyope testing
 *)

(* Load OPEdefs.m *)
Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];

(* Helper function to format output *)
FormatResult[label_, result_] := Module[{},
    Print["========================================"];
    Print["Test: ", label];
    Print["Output: ", result];
    Print[""];
];

Print[""];
Print["######## VIRASORO ALGEBRA ########"];
Print[""];

(* Define Virasoro generator *)
Bosonic[T]

(* Define central charge *)
c = Symbol["c"];

(* Define Virasoro OPE: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w) *)
OPE[T, T] = MakeOPE[(c One[w])/(2(z - w)^4) + (2 T[w])/(z - w)^2 + Derivative[1][T][w]/(z - w) + Ord[z, w, 0]];

(* ============================================ *)
(* Test Case 1: Simple nested NO *)
(* ============================================ *)

(* Test 1.1: NO(NO(T, T), T) - left nested *)
expr1 = NO[NO[T, T], T];
result1 = Expand[expr1];
FormatResult["VIR-1.1: NO(NO(T, T), T)", result1];

(* Test 1.2: NO(T, NO(T, T)) - right nested *)
expr2 = NO[T, NO[T, T]];
result2 = Expand[expr2];
FormatResult["VIR-1.2: NO(T, NO(T, T))", result2];

(* ============================================ *)
(* Test Case 2: Nested NO with derivatives *)
(* ============================================ *)

(* Test 2.1: NO(NO(T, ∂T), T) *)
expr3 = NO[NO[T, Derivative[1][T]], T];
result3 = Expand[expr3];
FormatResult["VIR-2.1: NO(NO(T, dT), T)", result3];

(* Test 2.2: NO(T, NO(∂T, T)) *)
expr4 = NO[T, NO[Derivative[1][T], T]];
result4 = Expand[expr4];
FormatResult["VIR-2.2: NO(T, NO(dT, T))", result4];

(* Test 2.3: NO(NO(∂T, T), ∂T) *)
expr5 = NO[NO[Derivative[1][T], T], Derivative[1][T]];
result5 = Expand[expr5];
FormatResult["VIR-2.3: NO(NO(dT, T), dT)", result5];

(* Test 2.4: NO(∂T, NO(T, ∂T)) *)
expr6 = NO[Derivative[1][T], NO[T, Derivative[1][T]]];
result6 = Expand[expr6];
FormatResult["VIR-2.4: NO(dT, NO(T, dT))", result6];

(* ============================================ *)
(* Test Case 3: Triple nested NO *)
(* ============================================ *)

(* Test 3.1: NO(NO(NO(T, T), T), T) *)
expr7 = NO[NO[NO[T, T], T], T];
result7 = Expand[expr7];
FormatResult["VIR-3.1: NO(NO(NO(T, T), T), T)", result7];

(* Test 3.2: NO(T, NO(T, NO(T, T))) *)
expr8 = NO[T, NO[T, NO[T, T]]];
result8 = Expand[expr8];
FormatResult["VIR-3.2: NO(T, NO(T, NO(T, T)))", result8];

(* Test 3.3: NO(NO(T, NO(T, T)), T) *)
expr9 = NO[NO[T, NO[T, T]], T];
result9 = Expand[expr9];
FormatResult["VIR-3.3: NO(NO(T, NO(T, T)), T)", result9];

(* ============================================ *)
(* Test Case 4: Higher derivatives *)
(* ============================================ *)

(* Test 4.1: NO(NO(∂²T, T), T) *)
expr10 = NO[NO[Derivative[2][T], T], T];
result10 = Expand[expr10];
FormatResult["VIR-4.1: NO(NO(d²T, T), T)", result10];

(* Test 4.2: NO(T, NO(T, ∂²T)) *)
expr11 = NO[T, NO[T, Derivative[2][T]]];
result11 = Expand[expr11];
FormatResult["VIR-4.2: NO(T, NO(T, d²T))", result11];

(* Test 4.3: NO(NO(∂T, ∂T), T) *)
expr12 = NO[NO[Derivative[1][T], Derivative[1][T]], T];
result12 = Expand[expr12];
FormatResult["VIR-4.3: NO(NO(dT, dT), T)", result12];

(* Test 4.4: NO(∂²T, NO(∂T, T)) *)
expr13 = NO[Derivative[2][T], NO[Derivative[1][T], T]];
result13 = Expand[expr13];
FormatResult["VIR-4.4: NO(d²T, NO(dT, T))", result13];

(* ============================================ *)
(* Test Case 5: Complex nested structures *)
(* ============================================ *)

(* Test 5.1: NO(NO(T, ∂T), NO(T, ∂T)) *)
expr14 = NO[NO[T, Derivative[1][T]], NO[T, Derivative[1][T]]];
result14 = Expand[expr14];
FormatResult["VIR-5.1: NO(NO(T, dT), NO(T, dT))", result14];

(* Test 5.2: NO(NO(NO(∂T, T), ∂T), T) *)
expr15 = NO[NO[NO[Derivative[1][T], T], Derivative[1][T]], T];
result15 = Expand[expr15];
FormatResult["VIR-5.2: NO(NO(NO(dT, T), dT), T)", result15];

Print[""];
Print["######## ALL VIRASORO TESTS COMPLETED ########"];
Print[""];
