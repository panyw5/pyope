#!/usr/bin/env wolframscript
(* 
 * sl(2) Kac-Moody Algebra Nested NO Expansion Tests
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
Print["######## SL(2) KAC-MOODY ALGEBRA ########"];
Print[""];

(* Define sl(2) generators *)
Bosonic[Jplus, Jzero, Jminus]

(* Define level parameter *)
k = Symbol["k"];

(* Define OPEs for sl(2) at level k *)
OPE[Jplus, Jzero] = MakeOPE[Jplus[w]/(z - w) + Ord[z, w, 0]];
OPE[Jzero, Jzero] = MakeOPE[(k/2 One[w])/(z - w)^2 + Ord[z, w, 0]];
OPE[Jzero, Jminus] = MakeOPE[-Jminus[w]/(z - w) + Ord[z, w, 0]];
OPE[Jplus, Jminus] = MakeOPE[(k One[w])/(z - w)^2 + (2 Jzero[w])/(z - w) + Ord[z, w, 0]];

(* ============================================ *)
(* Test Case 1: Simple nested NO *)
(* ============================================ *)

(* Test 1.1: NO(NO(Jplus, Jzero), Jminus) *)
expr1 = NO[NO[Jplus, Jzero], Jminus];
result1 = Expand[expr1];
FormatResult["SL2-1.1: NO(NO(J+, J0), J-)", result1];

(* Test 1.2: NO(Jplus, NO(Jzero, Jminus)) *)
expr2 = NO[Jplus, NO[Jzero, Jminus]];
result2 = Expand[expr2];
FormatResult["SL2-1.2: NO(J+, NO(J0, J-))", result2];

(* Test 1.3: NO(NO(Jminus, Jplus), Jzero) *)
expr3 = NO[NO[Jminus, Jplus], Jzero];
result3 = Expand[expr3];
FormatResult["SL2-1.3: NO(NO(J-, J+), J0)", result3];

(* Test 1.4: NO(Jminus, NO(Jplus, Jzero)) *)
expr4 = NO[Jminus, NO[Jplus, Jzero]];
result4 = Expand[expr4];
FormatResult["SL2-1.4: NO(J-, NO(J+, J0))", result4];

(* Test 1.5: NO(NO(Jzero, Jplus), Jminus) *)
expr5 = NO[NO[Jzero, Jplus], Jminus];
result5 = Expand[expr5];
FormatResult["SL2-1.5: NO(NO(J0, J+), J-)", result5];

(* Test 1.6: NO(Jzero, NO(Jminus, Jplus)) *)
expr6 = NO[Jzero, NO[Jminus, Jplus]];
result6 = Expand[expr6];
FormatResult["SL2-1.6: NO(J0, NO(J-, J+))", result6];

(* ============================================ *)
(* Test Case 2: Triple nested NO *)
(* ============================================ *)

(* Test 2.1: NO(NO(Jminus, NO(Jplus, Jzero))) *)
expr7 = NO[NO[Jminus, NO[Jplus, Jzero]]];
result7 = Expand[expr7];
FormatResult["SL2-2.1: NO(NO(J-, NO(J+, J0)))", result7];

(* Test 2.2: NO(NO(Jzero, NO(Jplus, Jminus))) *)
expr8 = NO[NO[Jzero, NO[Jplus, Jminus]]];
result8 = Expand[expr8];
FormatResult["SL2-2.2: NO(NO(J0, NO(J+, J-)))", result8];

(* Test 2.3: NO(NO(NO(Jplus, Jzero), Jminus), Jzero) *)
expr9 = NO[NO[NO[Jplus, Jzero], Jminus], Jzero];
result9 = Expand[expr9];
FormatResult["SL2-2.3: NO(NO(NO(J+, J0), J-), J0)", result9];

(* Test 2.4: NO(Jplus, NO(Jzero, NO(Jminus, Jplus))) *)
expr10 = NO[Jplus, NO[Jzero, NO[Jminus, Jplus]]];
result10 = Expand[expr10];
FormatResult["SL2-2.4: NO(J+, NO(J0, NO(J-, J+)))", result10];

(* ============================================ *)
(* Test Case 3: Nested NO with derivatives *)
(* ============================================ *)

(* Test 3.1: NO(NO(∂Jplus, Jzero), Jminus) *)
expr11 = NO[NO[Derivative[1][Jplus], Jzero], Jminus];
result11 = Expand[expr11];
FormatResult["SL2-3.1: NO(NO(dJ+, J0), J-)", result11];

(* Test 3.2: NO(Jplus, NO(∂Jzero, Jminus)) *)
expr12 = NO[Jplus, NO[Derivative[1][Jzero], Jminus]];
result12 = Expand[expr12];
FormatResult["SL2-3.2: NO(J+, NO(dJ0, J-))", result12];

(* Test 3.3: NO(NO(Jplus, ∂Jzero), Jminus) *)
expr13 = NO[NO[Jplus, Derivative[1][Jzero]], Jminus];
result13 = Expand[expr13];
FormatResult["SL2-3.3: NO(NO(J+, dJ0), J-)", result13];

(* Test 3.4: NO(NO(Jminus, Jplus), ∂Jzero) *)
expr14 = NO[NO[Jminus, Jplus], Derivative[1][Jzero]];
result14 = Expand[expr14];
FormatResult["SL2-3.4: NO(NO(J-, J+), dJ0)", result14];

(* Test 3.5: NO(∂Jminus, NO(Jplus, Jzero)) *)
expr15 = NO[Derivative[1][Jminus], NO[Jplus, Jzero]];
result15 = Expand[expr15];
FormatResult["SL2-3.5: NO(dJ-, NO(J+, J0))", result15];

(* ============================================ *)
(* Test Case 4: Sugawara construction *)
(* ============================================ *)

(* Test 4.1: NO(NO(Jplus, Jminus), NO(Jplus, Jminus)) *)
expr16 = NO[NO[Jplus, Jminus], NO[Jplus, Jminus]];
result16 = Expand[expr16];
FormatResult["SL2-4.1: NO(NO(J+, J-), NO(J+, J-))", result16];

(* Test 4.2: NO(NO(Jzero, Jzero), NO(Jplus, Jminus)) *)
expr17 = NO[NO[Jzero, Jzero], NO[Jplus, Jminus]];
result17 = Expand[expr17];
FormatResult["SL2-4.2: NO(NO(J0, J0), NO(J+, J-))", result17];

(* Test 4.3: NO(NO(Jplus, Jminus), NO(Jzero, Jzero)) *)
expr18 = NO[NO[Jplus, Jminus], NO[Jzero, Jzero]];
result18 = Expand[expr18];
FormatResult["SL2-4.3: NO(NO(J+, J-), NO(J0, J0))", result18];

(* ============================================ *)
(* Test Case 5: Complex nested structures *)
(* ============================================ *)

(* Test 5.1: NO(NO(NO(Jplus, Jminus), Jzero), Jplus) *)
expr19 = NO[NO[NO[Jplus, Jminus], Jzero], Jplus];
result19 = Expand[expr19];
FormatResult["SL2-5.1: NO(NO(NO(J+, J-), J0), J+)", result19];

(* Test 5.2: NO(NO(Jplus, NO(Jminus, Jzero)), Jplus) *)
expr20 = NO[NO[Jplus, NO[Jminus, Jzero]], Jplus];
result20 = Expand[expr20];
FormatResult["SL2-5.2: NO(NO(J+, NO(J-, J0)), J+)", result20];

(* Test 5.3: NO(NO(∂Jplus, ∂Jminus), Jzero) *)
expr21 = NO[NO[Derivative[1][Jplus], Derivative[1][Jminus]], Jzero];
result21 = Expand[expr21];
FormatResult["SL2-5.3: NO(NO(dJ+, dJ-), J0)", result21];

(* Test 5.4: NO(NO(Jplus, ∂Jminus), ∂Jzero) *)
expr22 = NO[NO[Jplus, Derivative[1][Jminus]], Derivative[1][Jzero]];
result22 = Expand[expr22];
FormatResult["SL2-5.4: NO(NO(J+, dJ-), dJ0)", result22];

Print[""];
Print["######## ALL SL(2) TESTS COMPLETED ########"];
Print[""];
