#!/usr/bin/env wolframscript
(* Explanation: Why Nl[gg, gg] is excluded *)

Print["========================================"];
Print["Why Nl[gg[0,1,z], gg[0,1,z]] does NOT appear in final physical states"];
Print["========================================"];
Print[""];

Print["=== 1. Field Generation Stage ==="];
Print["At finalgen generation, the code DOES generate Nl[gg[0,1,z], gg[0,1,z]]"];
Print["This is because the code mechanically generates all possible operator combinations."];
Print["Position: finalgen[[19]] = Nl[gg[0,1,z], gg[0,1,z]]"];
Print[""];

Print["=== 2. Why this term is mathematically ZERO ==="];
Print["gg is a fermionic operator"];
Print["Fermions satisfy anticommutation: {psi(z), psi(w)} = 0 when z != w"];
Print["At the same point, normal ordered product of fermions:"];
Print["  :psi(z) psi(z): = (1/2){psi(z), psi(z)} = 0"];
Print[""];
Print["Therefore: Nl[gg[0,1,z], gg[0,1,z]] = :gg(z)*gg(z): = 0"];
Print[""];

Print["=== 3. How the code excludes this term ==="];
Print["The code does NOT exclude it at generation stage,"];
Print["but in subsequent steps:"];
Print[""];

Print["Method A: Automatic zero in OPE calculations"];
Print["  When computing OPEs involving Nl[gg, gg],"];
Print["  due to OPE[gg, gg] definition (fermionic anticommutation),"];
Print["  all terms containing Nl[gg, gg] automatically simplify to zero."];
Print[""];

Print["Method B: Identified in null state analysis"];
Print["  The 'Calculating Physical States' section will:"];
Print["  1. Generate all possible states (including Nl[gg, gg])"];
Print["  2. Compute relations under OPE action"];
Print["  3. Identify null states via matrix rank analysis"];
Print["  4. Nl[gg, gg] is identified as a null state (because it equals zero)"];
Print["  5. Physical states = naive states - null states"];
Print[""];

Print["=== 4. Code verification ==="];
Print["Check the OPE definition at line 231:"];
Print["  Fermionic[b1,c1,gg[i_,j_,z],ggt[i_,j_,z],ggw[i_,j_,z],ggwb[i_,j_,z]];"];
Print[""];
Print["This tells OPEdefs package: gg, ggt, ggw, ggwb are all fermions."];
Print["When computing NO[gg, gg], OPEdefs automatically applies anticommutation."];
Print[""];

Print["=== 5. Summary ==="];
Print["Nl[gg[0,1,z], gg[0,1,z]] EXISTS in finalgen (position 19)"];
Print["But it is excluded in subsequent steps because:"];
Print["  * It equals zero mathematically (fermionic anticommutation)"];
Print["  * OPEdefs package automatically handles fermionic statistics"];
Print["  * In physical state calculation, it is identified as a null state"];
Print[""];
Print["This is RUNTIME exclusion, not GENERATION-TIME exclusion."];
Print["========================================"];
