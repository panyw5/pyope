#!/usr/bin/env wolframscript
(* 
 * U(1) Current Algebra Nested NO Expansion Tests
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
Print["######## U(1) CURRENT ALGEBRA ########"];
Print[""];

(* Define generators *)
Bosonic[T, J]

(* Define central charge *)
c = Symbol["c"];

(* Define OPEs *)
OPE[J, J] = MakeOPE[One[w]/(z - w)^2 + Ord[z, w, 0]];
OPE[T, J] = MakeOPE[J[w]/(z - w)^2 + Derivative[1][J][w]/(z - w) + Ord[z, w, 0]];
OPE[T, T] = MakeOPE[(c/2 One[w])/(z - w)^4 + (2 T[w])/(z - w)^2 + Derivative[1][T][w]/(z - w) + Ord[z, w, 0]];

(* ============================================ *)
(* Test Case 1: Simple nested NO with J *)
(* ============================================ *)

(* Test 1.1: NO(NO(J, J), J) *)
expr1 = NO[NO[J, J], J];
result1 = Expand[expr1];
FormatResult["U1-1.1: NO(NO(J, J), J)", result1];

(* Test 1.2: NO(J, NO(J, J)) *)
expr2 = NO[J, NO[J, J]];
result2 = Expand[expr2];
FormatResult["U1-1.2: NO(J, NO(J, J))", result2];

(* Test 1.3: NO(NO(NO(J, J), J), J) *)
expr3 = NO[NO[NO[J, J], J], J];
result3 = Expand[expr3];
FormatResult["U1-1.3: NO(NO(NO(J, J), J), J)", result3];

(* Test 1.4: NO(J, NO(J, NO(J, J))) *)
expr4 = NO[J, NO[J, NO[J, J]]];
result4 = Expand[expr4];
FormatResult["U1-1.4: NO(J, NO(J, NO(J, J)))", result4];

(* ============================================ *)
(* Test Case 2: Mixed T and J *)
(* ============================================ *)

(* Test 2.1: NO(NO(T, J), J) *)
expr5 = NO[NO[T, J], J];
result5 = Expand[expr5];
FormatResult["U1-2.1: NO(NO(T, J), J)", result5];

(* Test 2.2: NO(T, NO(J, J)) *)
expr6 = NO[T, NO[J, J]];
result6 = Expand[expr6];
FormatResult["U1-2.2: NO(T, NO(J, J))", result6];

(* Test 2.3: NO(NO(J, T), J) *)
expr7 = NO[NO[J, T], J];
result7 = Expand[expr7];
FormatResult["U1-2.3: NO(NO(J, T), J)", result7];

(* Test 2.4: NO(J, NO(T, J)) *)
expr8 = NO[J, NO[T, J]];
result8 = Expand[expr8];
FormatResult["U1-2.4: NO(J, NO(T, J))", result8];

(* Test 2.5: NO(NO(T, T), J) *)
expr9 = NO[NO[T, T], J];
result9 = Expand[expr9];
FormatResult["U1-2.5: NO(NO(T, T), J)", result9];

(* Test 2.6: NO(T, NO(T, J)) *)
expr10 = NO[T, NO[T, J]];
result10 = Expand[expr10];
FormatResult["U1-2.6: NO(T, NO(T, J))", result10];

(* ============================================ *)
(* Test Case 3: Nested NO with derivatives of J *)
(* ============================================ *)

(* Test 3.1: NO(NO(∂J, J), J) *)
expr11 = NO[NO[Derivative[1][J], J], J];
result11 = Expand[expr11];
FormatResult["U1-3.1: NO(NO(dJ, J), J)", result11];

(* Test 3.2: NO(J, NO(∂J, J)) *)
expr12 = NO[J, NO[Derivative[1][J], J]];
result12 = Expand[expr12];
FormatResult["U1-3.2: NO(J, NO(dJ, J))", result12];

(* Test 3.3: NO(NO(J, ∂J), J) *)
expr13 = NO[NO[J, Derivative[1][J]], J];
result13 = Expand[expr13];
FormatResult["U1-3.3: NO(NO(J, dJ), J)", result13];

(* Test 3.4: NO(NO(∂J, ∂J), J) *)
expr14 = NO[NO[Derivative[1][J], Derivative[1][J]], J];
result14 = Expand[expr14];
FormatResult["U1-3.4: NO(NO(dJ, dJ), J)", result14];

(* Test 3.5: NO(∂J, NO(J, ∂J)) *)
expr15 = NO[Derivative[1][J], NO[J, Derivative[1][J]]];
result15 = Expand[expr15];
FormatResult["U1-3.5: NO(dJ, NO(J, dJ))", result15];

(* Test 3.6: NO(NO(∂²J, J), J) *)
expr16 = NO[NO[Derivative[2][J], J], J];
result16 = Expand[expr16];
FormatResult["U1-3.6: NO(NO(d²J, J), J)", result16];

(* ============================================ *)
(* Test Case 4: Mixed T, J with derivatives *)
(* ============================================ *)

(* Test 4.1: NO(NO(T, ∂J), J) *)
expr17 = NO[NO[T, Derivative[1][J]], J];
result17 = Expand[expr17];
FormatResult["U1-4.1: NO(NO(T, dJ), J)", result17];

(* Test 4.2: NO(T, NO(∂J, J)) *)
expr18 = NO[T, NO[Derivative[1][J], J]];
result18 = Expand[expr18];
FormatResult["U1-4.2: NO(T, NO(dJ, J))", result18];

(* Test 4.3: NO(NO(∂T, J), J) *)
expr19 = NO[NO[Derivative[1][T], J], J];
result19 = Expand[expr19];
FormatResult["U1-4.3: NO(NO(dT, J), J)", result19];

(* Test 4.4: NO(∂T, NO(J, J)) *)
expr20 = NO[Derivative[1][T], NO[J, J]];
result20 = Expand[expr20];
FormatResult["U1-4.4: NO(dT, NO(J, J))", result20];

(* Test 4.5: NO(NO(T, J), ∂J) *)
expr21 = NO[NO[T, J], Derivative[1][J]];
result21 = Expand[expr21];
FormatResult["U1-4.5: NO(NO(T, J), dJ)", result21];

(* ============================================ *)
(* Test Case 5: Complex nested structures *)
(* ============================================ *)

(* Test 5.1: NO(NO(NO(T, J), J), J) *)
expr22 = NO[NO[NO[T, J], J], J];
result22 = Expand[expr22];
FormatResult["U1-5.1: NO(NO(NO(T, J), J), J)", result22];

(* Test 5.2: NO(T, NO(J, NO(J, J))) *)
expr23 = NO[T, NO[J, NO[J, J]]];
result23 = Expand[expr23];
FormatResult["U1-5.2: NO(T, NO(J, NO(J, J)))", result23];

(* Test 5.3: NO(NO(T, NO(J, J)), J) *)
expr24 = NO[NO[T, NO[J, J]], J];
result24 = Expand[expr24];
FormatResult["U1-5.3: NO(NO(T, NO(J, J)), J)", result24];

(* Test 5.4: NO(NO(NO(∂J, J), ∂J), J) *)
expr25 = NO[NO[NO[Derivative[1][J], J], Derivative[1][J]], J];
result25 = Expand[expr25];
FormatResult["U1-5.4: NO(NO(NO(dJ, J), dJ), J)", result25];

(* Test 5.5: NO(NO(T, ∂T), NO(J, ∂J)) *)
expr26 = NO[NO[T, Derivative[1][T]], NO[J, Derivative[1][J]]];
result26 = Expand[expr26];
FormatResult["U1-5.5: NO(NO(T, dT), NO(J, dJ))", result26];

Print[""];
Print["######## ALL U(1) TESTS COMPLETED ########"];
Print[""];
