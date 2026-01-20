#!/usr/bin/env wolframscript
(* 
 * N=4 Small Superconformal Algebra Nested NO Expansion Tests
 * Reference results for pyope testing
 *)

(* Load OPEdefs.m *)
Import["OPEdefs.m"];
Print["OPEdefs.m loaded successfully!"];
Print[""];

(* Helper function to format output *)
FormatResult[label_, result_] := Module[{},
    Print["========================================"];
    Print["Test: ", label];
    Print["Output: ", result];
    Print[""];
];

Print[""];
Print["######## N=4 SMALL SUPERCONFORMAL ALGEBRA ########"];
Print[""];

(* Define operators *)
Bosonic[β, γ];
Fermionic[b, c];

(* Define OPEs for bc and βγ systems *)
OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[β, γ] = MakeOPE[-(One[w]/(z - w)) + Ord[z, w, 0]];

(* Define N=4 small SCFA generators *)
Jplus = β;
J0 = NO[b, c] + 2 NO[β, γ];
Jminus = -NO[β, NO[γ, γ]] - NO[γ, NO[b, c]] + 3/2 Derivative[1][γ];

Gplus = b;
Gminus = NO[b, γ];
Gtplus = NO[c, Derivative[1][β]] + 2 NO[Derivative[1][c], β];
Gtminus = -NO[b, Derivative[1][c], c] + 2 NO[β, γ, Derivative[1][c]] + NO[Derivative[1][β], γ, c] - 3/2 Derivative[1][Derivative[1][c]];

T = -3/2 NO[b, Derivative[1][c]] - NO[β, Derivative[1][γ]] - 1/2 NO[Derivative[1][b], c];

(* ============================================ *)
(* Test Case 1: Nested NO with J generators *)
(* ============================================ *)

(* Test 1.1: NO(NO(Jplus, J0), Jminus) *)
expr1 = NO[NO[Jplus, J0], Jminus];
result1 = Expand[expr1];
FormatResult["N4-1.1: NO(NO(J+, J0), J-)", result1];

(* Test 1.2: NO(Jplus, NO(J0, Jminus)) *)
expr2 = NO[Jplus, NO[J0, Jminus]];
result2 = Expand[expr2];
FormatResult["N4-1.2: NO(J+, NO(J0, J-))", result2];

(* Test 1.3: NO(NO(J0, J0), Jplus) *)
expr3 = NO[NO[J0, J0], Jplus];
result3 = Expand[expr3];
FormatResult["N4-1.3: NO(NO(J0, J0), J+)", result3];

(* Test 1.4: NO(J0, NO(Jplus, Jminus)) *)
expr4 = NO[J0, NO[Jplus, Jminus]];
result4 = Expand[expr4];
FormatResult["N4-1.4: NO(J0, NO(J+, J-))", result4];

(* ============================================ *)
(* Test Case 2: Nested NO with G generators *)
(* ============================================ *)

(* Test 2.1: NO(NO(Gplus, Gminus), Gplus) *)
expr5 = NO[NO[Gplus, Gminus], Gplus];
result5 = Expand[expr5];
FormatResult["N4-2.1: NO(NO(G+, G-), G+)", result5];

(* Test 2.2: NO(Gplus, NO(Gminus, Gplus)) *)
expr6 = NO[Gplus, NO[Gminus, Gplus]];
result6 = Expand[expr6];
FormatResult["N4-2.2: NO(G+, NO(G-, G+))", result6];

(* Test 2.3: NO(NO(Gtplus, Gtminus), Gtplus) *)
expr7 = NO[NO[Gtplus, Gtminus], Gtplus];
result7 = Expand[expr7];
FormatResult["N4-2.3: NO(NO(Gt+, Gt-), Gt+)", result7];

(* Test 2.4: NO(Gtplus, NO(Gtminus, Gtplus)) *)
expr8 = NO[Gtplus, NO[Gtminus, Gtplus]];
result8 = Expand[expr8];
FormatResult("N4-2.4: NO(Gt+, NO(Gt-, Gt+))", result8];

(* ============================================ *)
(* Test Case 3: Mixed J and G generators *)
(* ============================================ *)

(* Test 3.1: NO(NO(Jplus, Gplus), Gminus) *)
expr9 = NO[NO[Jplus, Gplus], Gminus];
result9 = Expand[expr9];
FormatResult["N4-3.1: NO(NO(J+, G+), G-)", result9];

(* Test 3.2: NO(Jplus, NO(Gplus, Gminus)) *)
expr10 = NO[Jplus, NO[Gplus, Gminus]];
result10 = Expand[expr10];
FormatResult["N4-3.2: NO(J+, NO(G+, G-))", result10];

(* Test 3.3: NO(NO(J0, Gplus), Gminus) *)
expr11 = NO[NO[J0, Gplus], Gminus];
result11 = Expand[expr11];
FormatResult["N4-3.3: NO(NO(J0, G+), G-)", result11];

(* Test 3.4: NO(J0, NO(Gtplus, Gtminus)) *)
expr12 = NO[J0, NO[Gtplus, Gtminus]];
result12 = Expand[expr12];
FormatResult["N4-3.4: NO(J0, NO(Gt+, Gt-))", result12];

(* ============================================ *)
(* Test Case 4: Nested NO with T *)
(* ============================================ *)

(* Test 4.1: NO(NO(T, Jplus), Jminus) *)
expr13 = NO[NO[T, Jplus], Jminus];
result13 = Expand[expr13];
FormatResult["N4-4.1: NO(NO(T, J+), J-)", result13];

(* Test 4.2: NO(T, NO(Jplus, Jminus)) *)
expr14 = NO[T, NO[Jplus, Jminus]];
result14 = Expand[expr14];
FormatResult["N4-4.2: NO(T, NO(J+, J-))", result14];

(* Test 4.3: NO(NO(T, Gplus), Gminus) *)
expr15 = NO[NO[T, Gplus], Gminus];
result15 = Expand[expr15];
FormatResult["N4-4.3: NO(NO(T, G+), G-)", result15];

(* Test 4.4: NO(T, NO(Gplus, Gminus)) *)
expr16 = NO[T, NO[Gplus, Gminus]];
result16 = Expand[expr16];
FormatResult["N4-4.4: NO(T, NO(G+, G-))", result16];

(* Test 4.5: NO(NO(T, T), Jplus) *)
expr17 = NO[NO[T, T], Jplus];
result17 = Expand[expr17];
FormatResult["N4-4.5: NO(NO(T, T), J+)", result17];

(* ============================================ *)
(* Test Case 5: Triple nested structures *)
(* ============================================ *)

(* Test 5.1: NO(NO(NO(Jplus, J0), Jminus), J0) *)
expr18 = NO[NO[NO[Jplus, J0], Jminus], J0];
result18 = Expand[expr18];
FormatResult["N4-5.1: NO(NO(NO(J+, J0), J-), J0)", result18];

(* Test 5.2: NO(Jplus, NO(J0, NO(Jminus, J0))) *)
expr19 = NO[Jplus, NO[J0, NO[Jminus, J0]]];
result19 = Expand[expr19];
FormatResult["N4-5.2: NO(J+, NO(J0, NO(J-, J0)))", result19];

(* Test 5.3: NO(NO(T, NO(Jplus, J0)), Jminus) *)
expr20 = NO[NO[T, NO[Jplus, J0]], Jminus];
result20 = Expand[expr20];
FormatResult["N4-5.3: NO(NO(T, NO(J+, J0)), J-)", result20];

(* Test 5.4: NO(NO(NO(Gplus, Gminus), Gtplus), Gtminus) *)
expr21 = NO[NO[NO[Gplus, Gminus], Gtplus], Gtminus];
result21 = Expand[expr21];
FormatResult["N4-5.4: NO(NO(NO(G+, G-), Gt+), Gt-)", result21];

(* ============================================ *)
(* Test Case 6: Complex nested with derivatives *)
(* ============================================ *)

(* Test 6.1: NO(NO(∂Jplus, J0), Jminus) *)
expr22 = NO[NO[Derivative[1][Jplus], J0], Jminus];
result22 = Expand[expr22];
FormatResult["N4-6.1: NO(NO(dJ+, J0), J-)", result22];

(* Test 6.2: NO(Jplus, NO(∂J0, Jminus)) *)
expr23 = NO[Jplus, NO[Derivative[1][J0], Jminus]];
result23 = Expand[expr23];
FormatResult["N4-6.2: NO(J+, NO(dJ0, J-))", result23];

(* Test 6.3: NO(NO(∂T, Gplus), Gminus) *)
expr24 = NO[NO[Derivative[1][T], Gplus], Gminus];
result24 = Expand[expr24];
FormatResult["N4-6.3: NO(NO(dT, G+), G-)", result24];

(* Test 6.4: NO(T, NO(∂Gplus, Gminus)) *)
expr25 = NO[T, NO[Derivative[1][Gplus], Gminus]];
result25 = Expand[expr25];
FormatResult["N4-6.4: NO(T, NO(dG+, G-))", result25];

(* ============================================ *)
(* Test Case 7: Sugawara-like constructions *)
(* ============================================ *)

(* Test 7.1: NO(NO(Jplus, Jminus), NO(J0, J0)) *)
expr26 = NO[NO[Jplus, Jminus], NO[J0, J0]];
result26 = Expand[expr26];
FormatResult["N4-7.1: NO(NO(J+, J-), NO(J0, J0))", result26];

(* Test 7.2: NO(NO(Gplus, Gminus), NO(Gtplus, Gtminus)) *)
expr27 = NO[NO[Gplus, Gminus], NO[Gtplus, Gtminus]];
result27 = Expand[expr27];
FormatResult["N4-7.2: NO(NO(G+, G-), NO(Gt+, Gt-))", result27];

(* Test 7.3: NO(NO(T, T), NO(J0, J0)) *)
expr28 = NO[NO[T, T], NO[J0, J0]];
result28 = Expand[expr28];
FormatResult["N4-7.3: NO(NO(T, T), NO(J0, J0))", result28];

Print[""];
Print["######## ALL N=4 TESTS COMPLETED ########"];
Print[""];
