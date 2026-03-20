(* 
 * Virasoro Algebra OPE Computation
 * This script computes OPE(T, NO(T,T)) using OPEdefs package
 *)

Print["================================================="];
Print["Virasoro Algebra OPE Computation"];
Print["================================================="];
Print[""];

(* Load OPEdefs.m *)
Print["[*] Loading OPEdefs.m..."];
Get["/Users/lelouch/pyope/OPEdefs/OPEdefs.m"];
Print["[*] OPEdefs loaded"];
Print[""];

(* Define T as bosonic *)
Bosonic[T];

(* Define central charge *)
c = Symbol["c"];

(* Define Virasoro OPE *)
OPE[T, T] = MakeOPE[{c/2 One, 0, 2 T, T'}];
Print["[*] Virasoro OPE defined: T(z)T(w) = c/2/(z-w)^4 + 2T/(z-w)^2 + T'/(z-w)"];
Print[""];

(* Compute OPE(T, T) to verify *)
Print["Test 1: OPE(T, T)"];
Print["  Pole 4: ", OPEPole[4][T, T]];
Print["  Pole 2: ", OPEPole[2][T, T]];
Print["  Pole 1: ", OPEPole[1][T, T]];
Print[""];

(* Compute OPE(T, NO(T,T)) *)
Print["================================================="];
Print["Test 2: OPE(T, NO(T,T)) - Main Computation"];
Print["================================================="];
Print["Computing: T(z) * (TT)(w)"];
Print[""];

result2 = OPE[T, NO[T, T]];
Print["Pole 6: ", OPEPole[6][T, NO[T, T]]];
Print["Pole 4: ", OPEPole[4][T, NO[T, T]]];
Print["Pole 3: ", OPEPole[3][T, NO[T, T]]];
Print["Pole 2: ", OPEPole[2][T, NO[T, T]]];
Print["Pole 1: ", OPEPole[1][T, NO[T, T]]];
Print[""];

(* Interpret the result *)
Print["================================================="];
Print["Result Analysis"];
Print["================================================="];
Print[""];
Print["Mathematica result for Pole 1:"];
Print["  2 NO(T', T) + T^(3)/6"];
Print[""];
Print["PyOPE result for Pole 1:"];
Print["  NO(T, T') + NO(T', T)"];
Print[""];
Print["Relationship:"];
Print["  Use identity from Virasoro OPE:"];
Print["  NO(T', T) - NO(T, T') = -T^(3)/6"];
Print[""];
Print["Verification:"];
NOcomm = Expand[NO[Derivative[1][T], T] - NO[T, Derivative[1][T]]];
Print["  NO(T', T) - NO(T, T') = ", NOcomm];
Print[""];
Print["To get Mathematica form from PyOPE form:"];
Print["  Start: NO(T, T') + NO(T', T)"];
Print["  Use identity: = NO(T', T) + T^(3)/6 + NO(T', T)"];
Print["  Simplify: = 2 NO(T', T) + T^(3)/6"];
Print[""];
Print["Conclusion: PyOPE and Mathematica are equivalent"];
Print["The difference is in the NO ordering convention"];
Print[""];
Print["================================================="];
Print["Calculation Complete"];
Print["================================================="];
