#!/usr/bin/env wolframscript
(* 
 * bc-βγ Free Field System Nested NO Expansion Tests
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
Print["######## bc-βγ FREE FIELD SYSTEM ########"];
Print[""];

(* Define operators *)
Bosonic[β, γ];
Fermionic[b, c];

(* Define OPEs *)
OPE[b, c] = MakeOPE[One[w]/(z - w) + Ord[z, w, 0]];
OPE[β, γ] = MakeOPE[-(One[w]/(z - w)) + Ord[z, w, 0]];

(* ============================================ *)
(* Test Case 1: Simple nested NO with bc *)
(* ============================================ *)

(* Test 1.1: NO(NO(b, c), b) *)
expr1 = NO[NO[b, c], b];
result1 = Expand[expr1];
FormatResult["BC-1.1: NO(NO(b, c), b)", result1];

(* Test 1.2: NO(b, NO(c, b)) *)
expr2 = NO[b, NO[c, b]];
result2 = Expand[expr2];
FormatResult["BC-1.2: NO(b, NO(c, b))", result2];

(* Test 1.3: NO(NO(c, b), c) *)
expr3 = NO[NO[c, b], c];
result3 = Expand[expr3];
FormatResult["BC-1.3: NO(NO(c, b), c)", result3];

(* Test 1.4: NO(c, NO(b, c)) *)
expr4 = NO[c, NO[b, c]];
result4 = Expand[expr4];
FormatResult["BC-1.4: NO(c, NO(b, c))", result4];

(* Test 1.5: NO(NO(NO(b, c), b), c) *)
expr5 = NO[NO[NO[b, c], b], c];
result5 = Expand[expr5];
FormatResult["BC-1.5: NO(NO(NO(b, c), b), c)", result5];

(* ============================================ *)
(* Test Case 2: Simple nested NO with βγ *)
(* ============================================ *)

(* Test 2.1: NO(NO(β, γ), β) *)
expr6 = NO[NO[β, γ], β];
result6 = Expand[expr6];
FormatResult["BG-2.1: NO(NO(β, γ), β)", result6];

(* Test 2.2: NO(β, NO(γ, β)) *)
expr7 = NO[β, NO[γ, β]];
result7 = Expand[expr7];
FormatResult["BG-2.2: NO(β, NO(γ, β))", result7];

(* Test 2.3: NO(NO(γ, β), γ) *)
expr8 = NO[NO[γ, β], γ];
result8 = Expand[expr8];
FormatResult["BG-2.3: NO(NO(γ, β), γ)", result8];

(* Test 2.4: NO(γ, NO(β, γ)) *)
expr9 = NO[γ, NO[β, γ]];
result9 = Expand[expr9];
FormatResult["BG-2.4: NO(γ, NO(β, γ))", result9];

(* Test 2.5: NO(NO(NO(β, γ), β), γ) *)
expr10 = NO[NO[NO[β, γ], β], γ];
result10 = Expand[expr10];
FormatResult["BG-2.5: NO(NO(NO(β, γ), β), γ)", result10];

(* ============================================ *)
(* Test Case 3: Nested NO with derivatives *)
(* ============================================ *)

(* Test 3.1: NO(NO(∂b, c), b) *)
expr11 = NO[NO[Derivative[1][b], c], b];
result11 = Expand[expr11];
FormatResult["BC-3.1: NO(NO(db, c), b)", result11];

(* Test 3.2: NO(b, NO(∂c, b)) *)
expr12 = NO[b, NO[Derivative[1][c], b]];
result12 = Expand[expr12];
FormatResult["BC-3.2: NO(b, NO(dc, b))", result12];

(* Test 3.3: NO(NO(b, ∂c), b) *)
expr13 = NO[NO[b, Derivative[1][c]], b];
result13 = Expand[expr13];
FormatResult["BC-3.3: NO(NO(b, dc), b)", result13];

(* Test 3.4: NO(NO(∂β, γ), β) *)
expr14 = NO[NO[Derivative[1][β], γ], β];
result14 = Expand[expr14];
FormatResult["BG-3.4: NO(NO(dβ, γ), β)", result14];

(* Test 3.5: NO(β, NO(∂γ, β)) *)
expr15 = NO[β, NO[Derivative[1][γ], β]];
result15 = Expand[expr15];
FormatResult["BG-3.5: NO(β, NO(dγ, β))", result15];

(* Test 3.6: NO(NO(β, ∂γ), β) *)
expr16 = NO[NO[β, Derivative[1][γ]], β];
result16 = Expand[expr16];
FormatResult["BG-3.6: NO(NO(β, dγ), β)", result16];

(* ============================================ *)
(* Test Case 4: Mixed bc and βγ *)
(* ============================================ *)

(* Test 4.1: NO(NO(b, c), β) *)
expr17 = NO[NO[b, c], β];
result17 = Expand[expr17];
FormatResult["MIX-4.1: NO(NO(b, c), β)", result17];

(* Test 4.2: NO(b, NO(c, β)) *)
expr18 = NO[b, NO[c, β]];
result18 = Expand[expr18];
FormatResult["MIX-4.2: NO(b, NO(c, β))", result18];

(* Test 4.3: NO(NO(β, γ), b) *)
expr19 = NO[NO[β, γ], b];
result19 = Expand[expr19];
FormatResult["MIX-4.3: NO(NO(β, γ), b)", result19];

(* Test 4.4: NO(β, NO(γ, b)) *)
expr20 = NO[β, NO[γ, b]];
result20 = Expand[expr20];
FormatResult["MIX-4.4: NO(β, NO(γ, b))", result20];

(* Test 4.5: NO(NO(b, β), NO(c, γ)) *)
expr21 = NO[NO[b, β], NO[c, γ]];
result21 = Expand[expr21];
FormatResult["MIX-4.5: NO(NO(b, β), NO(c, γ))", result21];

(* ============================================ *)
(* Test Case 5: Higher derivatives *)
(* ============================================ *)

(* Test 5.1: NO(NO(∂²b, c), b) *)
expr22 = NO[NO[Derivative[2][b], c], b];
result22 = Expand[expr22];
FormatResult["BC-5.1: NO(NO(d²b, c), b)", result22];

(* Test 5.2: NO(NO(∂b, ∂c), b) *)
expr23 = NO[NO[Derivative[1][b], Derivative[1][c]], b];
result23 = Expand[expr23];
FormatResult["BC-5.2: NO(NO(db, dc), b)", result23];

(* Test 5.3: NO(NO(∂²β, γ), β) *)
expr24 = NO[NO[Derivative[2][β], γ], β];
result24 = Expand[expr24];
FormatResult["BG-5.3: NO(NO(d²β, γ), β)", result24];

(* Test 5.4: NO(NO(∂β, ∂γ), β) *)
expr25 = NO[NO[Derivative[1][β], Derivative[1][γ]], β];
result25 = Expand[expr25];
FormatResult["BG-5.4: NO(NO(dβ, dγ), β)", result25];

(* Test 5.5: NO(∂b, NO(∂c, NO(b, c))) *)
expr26 = NO[Derivative[1][b], NO[Derivative[1][c], NO[b, c]]];
result26 = Expand[expr26];
FormatResult["BC-5.5: NO(db, NO(dc, NO(b, c)))", result26];

(* ============================================ *)
(* Test Case 6: Complex nested structures *)
(* ============================================ *)

(* Test 6.1: NO(NO(NO(b, c), NO(b, c)), b) *)
expr27 = NO[NO[NO[b, c], NO[b, c]], b];
result27 = Expand[expr27];
FormatResult["BC-6.1: NO(NO(NO(b, c), NO(b, c)), b)", result27];

(* Test 6.2: NO(NO(NO(β, γ), NO(β, γ)), β) *)
expr28 = NO[NO[NO[β, γ], NO[β, γ]], β];
result28 = Expand[expr28];
FormatResult["BG-6.2: NO(NO(NO(β, γ), NO(β, γ)), β)", result28];

(* Test 6.3: NO(NO(∂b, NO(c, b)), ∂c) *)
expr29 = NO[NO[Derivative[1][b], NO[c, b]], Derivative[1][c]];
result29 = Expand[expr29];
FormatResult["BC-6.3: NO(NO(db, NO(c, b)), dc)", result29];

(* Test 6.4: NO(NO(∂β, NO(γ, β)), ∂γ) *)
expr30 = NO[NO[Derivative[1][β], NO[γ, β]], Derivative[1][γ]];
result30 = Expand[expr30];
FormatResult["BG-6.4: NO(NO(dβ, NO(γ, β)), dγ)", result30];

Print[""];
Print["######## ALL bc-βγ TESTS COMPLETED ########"];
Print[""];
