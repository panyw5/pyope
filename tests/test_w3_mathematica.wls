#!/usr/bin/env wolframscript
(* 
 * W3 Algebra Nested NO Expansion Tests
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
Print["######## W3 ALGEBRA ########"];
Print[""];

(* Define operators *)
Bosonic[T, W]

(* Define central charge *)
c = Symbol["c"];

(* Define parameter β *)
β = Symbol["β"];

(* Define Λ operator *)
Λ = NO[T, T] - 3/10 Derivative[1][Derivative[1][T]];

(* Define W3 algebra OPEs *)
OPE[T, T] = MakeOPE[(c One[w])/(2(z - w)^4) + (2 T[w])/(z - w)^2 + Derivative[1][T][w]/(z - w) + Ord[z, w, 0]];

OPE[T, W] = MakeOPE[(3 W[w])/(z - w)^2 + Derivative[1][W][w]/(z - w) + Ord[z, w, 0]];

OPE[W, W] = MakeOPE[(c One[w])/(z - w)^6 + (2 T[w])/(z - w)^4 + Derivative[1][T][w]/(z - w)^3 + 
    1/(z - w)^2 (2 β Λ[w] + 3/10 Derivative[1][Derivative[1][T]][w]) +
    1/(z - w) (β Derivative[1][Λ][w] + 1/15 Derivative[1][Derivative[1][Derivative[1][T]]][w]) + Ord[z, w, 0]];

(* ============================================ *)
(* Test Case 1: Simple nested NO with T and W *)
(* ============================================ *)

(* Test 1.1: NO(NO(T, W), T) *)
expr1 = NO[NO[T, W], T];
result1 = Expand[expr1];
FormatResult["W3-1.1: NO(NO(T, W), T)", result1];

(* Test 1.2: NO(T, NO(W, T)) *)
expr2 = NO[T, NO[W, T]];
result2 = Expand[expr2];
FormatResult["W3-1.2: NO(T, NO(W, T))", result2];

(* Test 1.3: NO(NO(W, T), W) *)
expr3 = NO[NO[W, T], W];
result3 = Expand[expr3];
FormatResult["W3-1.3: NO(NO(W, T), W)", result3];

(* Test 1.4: NO(W, NO(T, W)) *)
expr4 = NO[W, NO[T, W]];
result4 = Expand[expr4];
FormatResult["W3-1.4: NO(W, NO(T, W))", result4];

(* Test 1.5: NO(NO(W, W), T) *)
expr5 = NO[NO[W, W], T];
result5 = Expand[expr5];
FormatResult["W3-1.5: NO(NO(W, W), T)", result5];

(* Test 1.6: NO(W, NO(W, T)) *)
expr6 = NO[W, NO[W, T]];
result6 = Expand[expr6];
FormatResult["W3-1.6: NO(W, NO(W, T))", result6];

(* ============================================ *)
(* Test Case 2: Nested NO with T only *)
(* ============================================ *)

(* Test 2.1: NO(NO(T, T), W) *)
expr7 = NO[NO[T, T], W];
result7 = Expand[expr7];
FormatResult["W3-2.1: NO(NO(T, T), W)", result7];

(* Test 2.2: NO(T, NO(T, W)) *)
expr8 = NO[T, NO[T, W]];
result8 = Expand[expr8];
FormatResult["W3-2.2: NO(T, NO(T, W))", result8];

(* Test 2.3: NO(NO(NO(T, T), T), W) *)
expr9 = NO[NO[NO[T, T], T], W];
result9 = Expand[expr9];
FormatResult["W3-2.3: NO(NO(NO(T, T), T), W)", result9];

(* Test 2.4: NO(T, NO(T, NO(T, W))) *)
expr10 = NO[T, NO[T, NO[T, W]]];
result10 = Expand[expr10];
FormatResult["W3-2.4: NO(T, NO(T, NO(T, W)))", result10];

(* ============================================ *)
(* Test Case 3: Nested NO with derivatives *)
(* ============================================ *)

(* Test 3.1: NO(NO(∂T, W), T) *)
expr11 = NO[NO[Derivative[1][T], W], T];
result11 = Expand[expr11];
FormatResult["W3-3.1: NO(NO(dT, W), T)", result11];

(* Test 3.2: NO(T, NO(∂W, T)) *)
expr12 = NO[T, NO[Derivative[1][W], T]];
result12 = Expand[expr12];
FormatResult["W3-3.2: NO(T, NO(dW, T))", result12];

(* Test 3.3: NO(NO(T, ∂W), T) *)
expr13 = NO[NO[T, Derivative[1][W]], T];
result13 = Expand[expr13];
FormatResult["W3-3.3: NO(NO(T, dW), T)", result13];

(* Test 3.4: NO(NO(∂W, T), W) *)
expr14 = NO[NO[Derivative[1][W], T], W];
result14 = Expand[expr14];
FormatResult["W3-3.4: NO(NO(dW, T), W)", result14];

(* Test 3.5: NO(∂W, NO(T, W)) *)
expr15 = NO[Derivative[1][W], NO[T, W]];
result15 = Expand[expr15];
FormatResult["W3-3.5: NO(dW, NO(T, W))", result15];

(* Test 3.6: NO(NO(∂T, ∂W), T) *)
expr16 = NO[NO[Derivative[1][T], Derivative[1][W]], T];
result16 = Expand[expr16];
FormatResult["W3-3.6: NO(NO(dT, dW), T)", result16];

(* ============================================ *)
(* Test Case 4: Triple nested NO *)
(* ============================================ *)

(* Test 4.1: NO(NO(NO(T, W), T), W) *)
expr17 = NO[NO[NO[T, W], T], W];
result17 = Expand[expr17];
FormatResult["W3-4.1: NO(NO(NO(T, W), T), W)", result17];

(* Test 4.2: NO(T, NO(W, NO(T, W))) *)
expr18 = NO[T, NO[W, NO[T, W]]];
result18 = Expand[expr18];
FormatResult["W3-4.2: NO(T, NO(W, NO(T, W)))", result18];

(* Test 4.3: NO(NO(T, NO(W, T)), W) *)
expr19 = NO[NO[T, NO[W, T]], W];
result19 = Expand[expr19];
FormatResult["W3-4.3: NO(NO(T, NO(W, T)), W)", result19];

(* Test 4.4: NO(NO(NO(W, W), T), T) *)
expr20 = NO[NO[NO[W, W], T], T];
result20 = Expand[expr20];
FormatResult["W3-4.4: NO(NO(NO(W, W), T), T)", result20];

(* ============================================ *)
(* Test Case 5: Higher derivatives *)
(* ============================================ *)

(* Test 5.1: NO(NO(∂²T, W), T) *)
expr21 = NO[NO[Derivative[2][T], W], T];
result21 = Expand[expr21];
FormatResult["W3-5.1: NO(NO(d²T, W), T)", result21];

(* Test 5.2: NO(T, NO(∂²W, T)) *)
expr22 = NO[T, NO[Derivative[2][W], T]];
result22 = Expand[expr22];
FormatResult["W3-5.2: NO(T, NO(d²W, T))", result22];

(* Test 5.3: NO(NO(∂W, ∂W), T) *)
expr23 = NO[NO[Derivative[1][W], Derivative[1][W]], T];
result23 = Expand[expr23];
FormatResult["W3-5.3: NO(NO(dW, dW), T)", result23];

(* Test 5.4: NO(∂²W, NO(T, W)) *)
expr24 = NO[Derivative[2][W], NO[T, W]];
result24 = Expand[expr24];
FormatResult["W3-5.4: NO(d²W, NO(T, W))", result24];

(* ============================================ *)
(* Test Case 6: Complex nested structures *)
(* ============================================ *)

(* Test 6.1: NO(NO(T, T), NO(W, W)) *)
expr25 = NO[NO[T, T], NO[W, W]];
result25 = Expand[expr25];
FormatResult["W3-6.1: NO(NO(T, T), NO(W, W))", result25];

(* Test 6.2: NO(NO(NO(T, W), NO(T, W)), T) *)
expr26 = NO[NO[NO[T, W], NO[T, W]], T];
result26 = Expand[expr26];
FormatResult["W3-6.2: NO(NO(NO(T, W), NO(T, W)), T)", result26];

(* Test 6.3: NO(NO(∂T, NO(W, T)), ∂W) *)
expr27 = NO[NO[Derivative[1][T], NO[W, T]], Derivative[1][W]];
result27 = Expand[expr27];
FormatResult["W3-6.3: NO(NO(dT, NO(W, T)), dW)", result27];

Print[""];
Print["######## ALL W3 TESTS COMPLETED ########"];
Print[""];
