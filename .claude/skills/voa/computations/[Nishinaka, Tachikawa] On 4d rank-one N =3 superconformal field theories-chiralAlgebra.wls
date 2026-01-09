#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*SetOptions[EvaluationNotebook[],Background->RGBColor[36/256,34/256,40/256]]*)


(* ::Section:: *)
(*Setting up*)


(* ::Text:: *)
(*In this notebook we construct the 2d chiral algebras associated to 4d rank-1 N=3 SCFTs.*)
(*Note that within this notebook what is referred by script l is denoted by k.*)
(*We use the Mathematica code "SOPEN2defs" developed by Krivonos and Thielemans [hep-th/9512029].*)
(*It can be obtained from one of its authors K. Thielemans (k.thielemans@ucl.ac.uk).*)
(*We assume it is saved in the same directory/folder as this notebook.*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*<<SOPEN2defs.m*)


(* ::Text:: *)
(*As  the helpful output of SOPEN2defs are sometimes annoying, we disable it:*)


(* ::Input:: *)
(*Unprotect[Print];Print[_]:={};Protect[Print];*)


(* ::Text:: *)
(*First we set up the OPE of the super energy-momentum tensor multiplet. This is done by first declaring the symbol to be used, and then listing the residues of the poles under Subscript[Z, 12] using MakeOPE:*)


(* ::Input:: *)
(*Bosonic[J]*)


(* ::Input:: *)
(*OPE[J,J]=MakeOPE[{c/3 One + t12**tb12**J,tb12**DBT[J]-t12**DT[J]+t12**tb12**J'}]*)


(* ::Text:: *)
(*Here, t12 is Subscript[\[Theta], 12],  tb12 is Subscript[\!\(\*OverscriptBox[\(\[Theta]\), \(_\)]\), 12], DT is the superderivative D, and DBT is \!\(\*OverscriptBox[\(D\), \(_\)]\). DDB which is [D,\!\(\*OverscriptBox[\(D\), \(_\)]\)], will also be used below.  Later it will be useful to have a function to extract the coefficients of Subscript[\[Theta], 12] etc:*)


(* ::Input:: *)
(*coeffOft12tb12[x_]:=x/. NonCommutativeMultiply[t12,tb12,y_]:> y*)


(* ::Input:: *)
(*coeffOft12[x_]:=coeffOft12tb12[tb12**x]*)


(* ::Input:: *)
(*coeffOftb12[x_]:=coeffOft12tb12[t12**x]*)


(* ::Input:: *)
(*coeffOf1[x_]:=coeffOft12tb12[t12**tb12**x]*)


(* ::Section:: *)
(*k=3*)


(* ::Text:: *)
(*Let us first consider the case k=3.*)


(* ::Text:: *)
(*A superprimary has the following OPE in general:*)


(* ::Input:: *)
(*primaryope[op_]:=MakeOPE[{\[CapitalDelta][op] t12**tb12**op, -F[op]op - t12** DT[op]+tb12**DBT[op]+ t12**tb12**op'}]*)


(* ::Text:: *)
(*Here we declare that W and \!\(\*OverscriptBox[\(W\), \(_\)]\) (denoted as WB below) are a chiral primary and an antichiral primary with the required dimensions. Note that somehow Krivonos and Thielemans defined (anti)chiralness oppositely from the usual convention.*)


(* ::Input:: *)
(*AChiralBosonic[W];\[CapitalDelta][W]=3/2;F[W]=3;*)


(* ::Input:: *)
(*ChiralBosonic[WB];\[CapitalDelta][WB]=3/2;F[WB]=-3;*)


(* ::Input:: *)
(*OPE[J,W]=primaryope[W]//OPESimplify*)


(* ::Input:: *)
(*OPE[J,WB]=primaryope[WB]//OPESimplify*)


(* ::Text:: *)
(*Now we set up the OPE of W and WB. As discussed in the paper, it should have the following form:*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(3/2) a  t12**tb12**One,*)
(*a One+(2/2) t12**tb12**J,*)
(* J +  t12**DT[J]*)
(*+(1/2) t12**tb12**(u1 NO[J,J]+u2 DDB[J]+u3 J'),*)
(*u1 NO[J,J]+u2  DDB[J]+u3  J' *)
(*+  t12**  DT[u1 NO[J,J]+u2 DDB[J]+u3 J']}]//OPESimplify;*)


(* ::Text:: *)
(*We fix the coefficients by demanding that the Jacobi identity for J, W and WB be satisfied. The failure is given by:*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Union*)


(* ::Text:: *)
(*Here and below, it is useful to ignore numerical coefficients before the terms that should vanish. So let us define a function to remove them:*)


(* ::Input:: *)
(*Prune[x_]:=PruneItem/@x*)


(* ::Input:: *)
(*PruneItem[x_]:=x*)


(* ::Input:: *)
(*PruneItem[Times[n_?NumberQ,xs__]]:=Times[xs]*)


(* ::Text:: *)
(*Then the failure of the Jacobi is given by*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Prune//Union*)


(* ::Text:: *)
(*which is much shorter. All terms can be made zero by the following choice:*)


(* ::Input:: *)
(*solution=Solve[GetCoefficients[%]==0]*)


(* ::Text:: *)
(*As explained in Krivonos-Thielemans, we then need to check the Jacobi for every permutation:*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{J,W,WB}]/.solution[[1]]//Simplify//Flatten//Union*)


(* ::Text:: *)
(*Good. Then we demand the Jacobi of W, W and WB to vanish:*)


(* ::Input:: *)
(*OPEJacobi[W,W,WB,Function->Factor]/.solution[[1]]//Simplify//Prune//Union*)


(* ::Text:: *)
(*The terms proportional to W can be made to vanish only when  c = \[Dash]15. We also see that we need to set the following operator as null:*)


(* ::Input:: *)
(*null=DT[Derivative[1][W]]-NO[J,DT[W]]+3 NO[DT[J],W];\[CapitalDelta][null]=3;F[null]=2;*)


(* ::Text:: *)
(*Now we study the Jacobi for all permutations of W, W and WB:*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{W,W,WB}]/.solution[[1]]/.{c-> -15,DT[W']-> Expand[DT[W']-null]}//Simplify//Flatten//Union*)


(* ::Text:: *)
(*The reminders are both descendants of null, as you can easily see:*)


(* ::Input:: *)
(*DT[null]/2*)


(* ::Input:: *)
(*2DBT[null]//Simplify*)


(* ::Text:: *)
(*We also need to check that the OPE of the null field with J, W, WB only produces operators with larger or equal dimensions. First we can check that this null field is a super Virasoro primary:*)


(* ::Input:: *)
(*OPE[J,null]-primaryope[null]/.solution[[1]]/.{c-> -15,DT[W']-> Expand[DT[W']-null]}//Simplify*)


(* ::Text:: *)
(*Second its OPE with W vanishes: *)


(* ::Input:: *)
(*OPE[W,null]/.solution[[1]]/.{c-> -15,DT[W']-> Expand[DT[W']-null]}//Simplify*)


(* ::Text:: *)
(*Finally its OPE with WB is:*)


(* ::Input:: *)
(*OPE[WB,null]/.solution[[1]]/.{c-> -15,DT[W']-> Expand[DT[W']-null]}//Simplify*)


(* ::Text:: *)
(*It only consists of operators with dimensions larger than that of the original null field. In particular, the coefficient of Subscript[\!\(\*OverscriptBox[\(\[Theta]\), \(_\)]\), 12]of the leading pole is*)


(* ::Input:: *)
(*relation=coeffOftb12[OPEPole[2][%]]//Expand*)


(* ::Text:: *)
(*which automatically realizes the expected Higgs branch relation (J^3)~ W \!\(\*OverscriptBox[\(W\), \(_\)]\) . Including the coefficients, we have J^3/4 = 9 W \!\(\*OverscriptBox[\(\(W\)\(\ \)\), \(_\)]\) up to descendants. *)


(* ::Section:: *)
(*An aside*)


(* ::Text:: *)
(*Let us reproduce the results of Odake and Inami-Matsuo-Yamanaka, where almost the same problem was analyzed; they took the bottom components of W and \!\(\*OverscriptBox[\(W\), \(_\)]\) to be fermionic, probably because fields with half-integral spin was assumed to be so. As OPEdefs do not like operators with the same name with different Grassmann assignment, we call them Z and ZB:*)


(* ::Input:: *)
(*AChiralFermionic[Z];\[CapitalDelta][Z]=3/2;F[Z]=3;*)


(* ::Input:: *)
(*ChiralFermionic[ZB];\[CapitalDelta][ZB]=3/2;F[ZB]=-3;*)


(* ::Input:: *)
(*OPE[J,Z]=primaryope[Z]//OPESimplify;*)


(* ::Input:: *)
(*OPE[J,ZB]=primaryope[ZB]//OPESimplify;*)


(* ::Text:: *)
(*The OPE of Z and ZB should be exactly as before:*)


(* ::Input:: *)
(*OPE[Z,ZB]=MakeOPE[{(3/2) a  t12**tb12**One,*)
(*a One+(2/2) t12**tb12**J,*)
(* J +  t12**DT[J]*)
(*+(1/2) t12**tb12**(u1 NO[J,J]+u2 DDB[J]+u3 J'),*)
(*u1 NO[J,J]+u2  DDB[J]+u3  J' *)
(*+  t12**  DT[u1 NO[J,J]+u2 DDB[J]+u3 J']}]/.solution[[1]]//OPESimplify;*)


(* ::Text:: *)
(*Just to be sure, let's check the Jacobi:*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{J,Z,ZB}]//Simplify//Flatten//Union*)


(* ::Text:: *)
(*Good. The the Jacobi of Z, Z and ZB gives*)


(* ::Input:: *)
(*OPEJacobi[Z,Z,ZB,Function->Factor]//Factor//Prune//Union*)


(* ::Text:: *)
(*We see we need to set c=9, as was shown in old papers:*)


(* ::Input:: *)
(*OPEJacobi[Z,Z,ZB,Function->Factor]/.c-> 9//Factor//Prune//Union*)


(* ::Text:: *)
(*The null field in this case should be*)


(* ::Input:: *)
(*nullF=NO[J,Z]+Z';\[CapitalDelta][nullF]=5/2;F[nullF]=3;*)


(* ::Text:: *)
(*The Jacobi identities in other permutations give:*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{Z,Z,ZB}]//.{c-> 9,Z'-> Expand[Z'-nullF],Z''-> Expand[Z'-nullF]' }//Factor//Flatten//Prune//Union*)


(* ::Text:: *)
(*We see that the remainders are indeed descendants of the null field:*)


(* ::Input:: *)
(*DBT[nullF]*)


(* ::Input:: *)
(*DDB[nullF]//.{Z'-> Expand[Z'-nullF],Z''-> Expand[Z'-nullF]' }//Simplify*)


(* ::Text:: *)
(*The null field is super Virasoro primary as before:*)


(* ::Input:: *)
(*OPE[J,nullF]-primaryope[nullF]/.c-> 9//Simplify*)


(* ::Text:: *)
(*And the OPE with Z and ZB only gives higher-dimension fields:*)


(* ::Input:: *)
(*OPE[Z,nullF]//Simplify*)


(* ::Input:: *)
(*OPE[ZB,nullF]/.c-> 9//Simplify*)


(* ::Section:: *)
(*k=2*)


(* ::Text:: *)
(*Before analyzing k = 4 and 6, it is instructive to consider the case k=2. First let's clear the internal state of OPEdefs:*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Text:: *)
(*We redefine the charges as appropriate for k = 2 :*)


(* ::Input:: *)
(*\[CapitalDelta][W]=1;F[W]=2;*)


(* ::Input:: *)
(*OPE[J,W]=primaryope[W]//OPESimplify*)


(* ::Input:: *)
(*\[CapitalDelta][WB]=1;F[WB]=-2;*)


(* ::Input:: *)
(*OPE[J,WB]=primaryope[WB]//OPESimplify*)


(* ::Text:: *)
(*The W - WB OPE  is very simple:*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(2/2)a t12**tb12**One,*)
(*a One+(1/2)t12**tb12**J,*)
(* J +  t12**DT[J]*)
(*}];*)


(* ::Text:: *)
(*The Jacobi for J, W, WB immediately fixes the coefficients.*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Prune//Union*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(2/2)a t12**tb12**One,*)
(*a One+(1/2)t12**tb12**J,*)
(* J +  t12**DT[J]*)
(*}]/.a-> -c/6//OPESimplify*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{J,W,WB}]//Simplify//Flatten//Union*)


(* ::Text:: *)
(*And then the Jacobi for W, W and WB is automatically satisfied, as it should be.*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{W,W,WB}]//Simplify//Flatten//Union*)


(* ::Text:: *)
(*Now let us demand that there is a relation of the form W \!\(\*OverscriptBox[\(W\), \(_\)]\)~ (J^2) .  The most general form is*)


(* ::Input:: *)
(*op=NO[W,WB]+a NO[J,J]+b J'+ d DDB[J];\[CapitalDelta][op]=2;F[op]=0;*)


(* ::Text:: *)
(*We demand that it is a super Virasoro primary :*)


(* ::Input:: *)
(*OPE[J,op]-primaryope[op]//GetCoefficients//Simplify*)


(* ::Input:: *)
(*solution=Solve[%==0,{a,b,d}]*)


(* ::Input:: *)
(*op=op/.solution[[1]]*)


(* ::Text:: *)
(*Its OPE with W is*)


(* ::Input:: *)
(*OPE[W,op]//FullSimplify*)


(* ::Text:: *)
(*from which we immediately see that we need c=\[Minus]9 or +6 for the "op" not to force W itself to be zero. When c=+6, we also see JW + (descendants) =0, which we do not want. So we conclude that in order to just impose the relation W \!\(\*OverscriptBox[\(W\), \(_\)]\)~ (J^2), we are required to take c=\[Minus]9, and this is the value that follows from the 4d analysis. The final relation then is*)


(* ::Input:: *)
(*op=op/.c-> -9*)


(* ::Text:: *)
(*So the relation is J^2 =4W \!\(\*OverscriptBox[\(W\), \(_\)]\)  up to descendants.*)


(* ::Section:: *)
(*k=4*)


(* ::Text:: *)
(*Now let us proceed to k=4.*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Text:: *)
(*It would be convenient to have a function that generates all possible operators of a given dimension constructed out of J only. A quick hack is as follows :*)


(* ::Input:: *)
(*X=(1-t J)(1-t^2 J') (1-t^3 J'') (1-t^4 J''') (1-t^5 J'''')(1-t^6 J''''') (1-t^(3/2) q DT[J])(1-t^(5/2) q DT[J'])(1-t^(7/2) q DT[J'']) **)
(*(1-t^(3/2) q^(-1) DBT[J])(1-t^(5/2) q^(-1) DBT[J'])(1-t^(7/2) q^(-1) DBT[J''])**)
(*(1-t^2 DDB[J])(1-t^3 DDB[J'])(1-t^4 DDB[J''])(1-t^5 DDB[J'''])(1-t^6 DDB[J'''']);*)


(* ::Input:: *)
(*XX=1/X+O[t]^7//Normal;*)


(* ::Input:: *)
(*monom[n_]:=MonomialList[Residue[Coefficient[XX,t^n] /q, {q, 0}]]*)


(* ::Input:: *)
(*NOy[x_]:=x;NOy[x__]:=NO[x];*)


(* ::Input:: *)
(*monomToNO[{x__}]:=NOy@@ (Join@@(monomToNOX/@{x}));*)


(* ::Input:: *)
(*monomToNO[x_]:=monomToNO[{x}];*)


(* ::Input:: *)
(*monomToNOX[Power[x_,n_]]:= Table[x,n];*)


(* ::Input:: *)
(*monomToNOX[x_]:= {x};*)


(* ::Input:: *)
(*SetOPEOptions[NOOrdering,1];*)


(* ::Input:: *)
(*ops[n_]:=ops[n]=monomToNO /@ (monom[n]/.Times-> List)*)


(* ::Text:: *)
(*We can test it:*)


(* ::Input:: *)
(*ops[1]*)


(* ::Input:: *)
(*ops[2]*)


(* ::Input:: *)
(*ops[3]*)


(* ::Text:: *)
(*Then the following gives a generic linear combination of dimension-n composite operators constructed of J:*)


(* ::Input:: *)
(*gL[n_,a_]:=Sum[ops[n][[i]] a[i],{i,1,Length[ops[n]]}]*)


(* ::Text:: *)
(*Now we set the dimensions of W and WB as appropriate for k=4:*)


(* ::Input:: *)
(*\[CapitalDelta][W]=2;F[W]=4;*)


(* ::Input:: *)
(*OPE[J,W]=primaryope[W]//OPESimplify*)


(* ::Input:: *)
(*\[CapitalDelta][WB]=2;F[WB]=-4;*)


(* ::Input:: *)
(*OPE[J,WB]=primaryope[WB]//OPESimplify*)


(* ::Text:: *)
(*The OPE of W and WB has the general form*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(4/2)a t12**tb12**One,*)
(*a One+(3/2)t12**tb12**J,*)
(*J +  t12**DT[J]*)
(*+(2/2)t12**tb12**gL[2,A],*)
(*gL[2,A] +  t12**  DT[gL[2,A]]*)
(*+(1/2)t12**tb12**gL[3,AA],*)
(*gL[3,AA]+t12**DT[gL[3,AA]]}]//OPESimplify;*)


(* ::Text:: *)
(*from which we can fix the coefficients by demanding the J-W-WB Jacobi:*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Union;*)


(* ::Input:: *)
(*solution=Solve[GetCoefficients[%]==0]*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{J,W,WB}]/.solution[[1]]//Simplify//Flatten//Union*)


(* ::Text:: *)
(*We then move on to the Jacobi of W, W and WB:*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]/.solution[[1]]//Factor//Prune//Union*)


(* ::Text:: *)
(*We do not study all the possibilities; for now let us content ourselves with the fact that the value c=\[Minus]21, which is expected from the 4d analysis, is an allowed value at this point. Let us further study the Jacobi identity by setting c to this value:*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]/.solution[[1]]/.c-> -21//Factor//Prune//Union*)


(* ::Text:: *)
(*So we need the following null field:*)


(* ::Input:: *)
(*null=DT[W']-NO[J,DT[W]]+4NO[DT[J],W];\[CapitalDelta][null]=7/2;F[null]=3;*)


(* ::Text:: *)
(*Note that it has the following descendant:*)


(* ::Input:: *)
(*DT[null]*)


(* ::Input:: *)
(*DBT[null]*)


(* ::Input:: *)
(*null'/5/. W''-> W''+DBT[null]//Expand*)


(* ::Text:: *)
(*So we can set up the following rule killing some of the null fields:*)


(* ::Input:: *)
(*nullRules={DT[W']-> DT[W']-null,DT[W'']-> DT[W'']-null',NO[DT[J],DT[W]]-> 0,W''-> W''+DBT[null],W'''-> W'''+DBT[null]',NO[DT[J],W']->NO[J',DT[W]]/2+NO[DDB[J],DT[W]]/2 };*)


(* ::Text:: *)
(*Then the Jacobi is now satisfied:*)


(* ::Input:: *)
(*(OPEJacobi@@#)&/@Permutations[{W,W,WB}]/.solution[[1]]/.c-> -21//.nullRules//Factor//Flatten//Union*)


(* ::Text:: *)
(*Note also that the OPE of the null field with J, W, WB contains operators with higher dimensions only:*)


(* ::Input:: *)
(*OPE[J,null]-primaryope[null]/.solution[[1]]/.c-> -21//Simplify*)


(* ::Input:: *)
(*OPE[W,null]/.solution[[1]]/.c-> -21//Simplify*)


(* ::Input:: *)
(*OPE[WB,null]/.solution[[1]]/.c-> -21//Simplify*)


(* ::Text:: *)
(*In particular, the last OPE means that the following WB descendant of the null field gives the Higgs branch relation:*)


(* ::Input:: *)
(*relation=coeffOftb12[OPEPole[2][%]]//Expand*)


(* ::Text:: *)
(*So the relation is J^2/30 =16W \!\(\*OverscriptBox[\(W\), \(_\)]\)  up to descendants.*)


(* ::Section:: *)
(*k=6*)


(* ::Text:: *)
(*The analysis for k=6 is entierly analogous, although more tedious.*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Input:: *)
(*\[CapitalDelta][W]=3;F[W]=6;*)


(* ::Input:: *)
(*OPE[J,W]=primaryope[W]//OPESimplify*)


(* ::Input:: *)
(*\[CapitalDelta][WB]=3;F[WB]=-6;*)


(* ::Input:: *)
(*OPE[J,WB]=primaryope[WB]//OPESimplify*)


(* ::Text:: *)
(*The OPE of W and WB has the general form*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(6/2)a t12**tb12**One,*)
(*a One+(5/2)t12**tb12**J,*)
(*J +  t12**DT[J]*)
(*+(4/2)t12**tb12**gL[2,A],*)
(*gL[2,A] +  t12**  DT[gL[2,A]]*)
(*+(3/2)t12**tb12**gL[3,AA],*)
(*gL[3,AA]+t12**DT[gL[3,AA]]*)
(*+(2/2) t12**tb12**gL[4,AAA],*)
(*gL[4,AAA]+t12**DT[gL[4,AAA]]*)
(*+(1/2)t12**tb12**gL[5,AAAA],*)
(*gL[5,AAAA]+t12**DT[gL[5,AAAA]]}]//OPESimplify;*)


(* ::Text:: *)
(*from which we can fix the coefficients by demanding the J-W-WB Jacobi:*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Union;*)


(* ::Input:: *)
(*solution=Solve[GetCoefficients[%]==0];*)


(* ::Text:: *)
(*We then move on to the Jacobi of W, W and WB:*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]/.solution[[1]]//Factor//Prune//Union;*)


(* ::Input:: *)
(*{remainder[[2]],remainder[[3]]}*)


(* ::Input:: *)
(*PolynomialGCD[-135-153 c+c^2,8424+12528 c-138 c^2-4 c^3+c^4]*)


(* ::Text:: *)
(*We see that c=\[Minus]33, \[Minus]15, 18 are allowed at this point. We just consider the case c=\[Minus]33. Now, to speed up the computation, we redefine the basic OPEs with the determined values:*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Input:: *)
(*OPE[J,J]=MakeOPE[{c/3 One + t12**tb12**J,tb12**DBT[J]-t12**DT[J]+t12**tb12**J'}]/.c-> -33;*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(6/2)a t12**tb12**One,*)
(*a One+(5/2)t12**tb12**J,*)
(*J +  t12**DT[J]*)
(*+(4/2)t12**tb12**gL[2,A],*)
(*gL[2,A] +  t12**  DT[gL[2,A]]*)
(*+(3/2)t12**tb12**gL[3,AA],*)
(*gL[3,AA]+t12**DT[gL[3,AA]]*)
(*+(2/2) t12**tb12**gL[4,AAA],*)
(*gL[4,AAA]+t12**DT[gL[4,AAA]]*)
(*+(1/2)t12**tb12**gL[5,AAAA],*)
(*gL[5,AAAA]+t12**DT[gL[5,AAAA]]}]/.solution[[1]]/.c-> -33//OPESimplify;*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]//Factor//Prune//Union;*)


(* ::Input:: *)
(*remainder[[2]]*)


(* ::Text:: *)
(*So we need the following null field:*)


(* ::Input:: *)
(*null=DT[W']-NO[J,DT[W]]+6NO[DT[J],W];\[CapitalDelta][null]=9/2;F[null]=5;*)


(* ::Text:: *)
(*The Jacobi identity for W, W, WB now follows if we use the following descendants. There should be a more systematic way to carry it out, but for now we just enumerate them in an ad hoc manner.*)


(* ::Text:: *)
(*First, the null field itself and its \!\(\*OverscriptBox[\(D\), \(_\)]\) leads to the following basic rules:*)


(* ::Input:: *)
(*DBT[null]//Expand*)


(* ::Input:: *)
(*basicRules={DT[W']-> DT[W']-null,*)
(*W''-> W''+DBT[null],W'''-> W'''+DBT[null]',W''''-> W''''+DBT[null]'',W'''''-> W'''''+DBT[null]'''}//Expand;*)


(* ::Text:: *)
(*Second, the D of the null field and its repeated derivatives lead to more rules:*)


(* ::Input:: *)
(*DT[null]*)


(* ::Input:: *)
(*DT[null]'//Expand*)


(* ::Input:: *)
(*DT[null]''//.basicRules//Expand*)


(* ::Input:: *)
(*DT[null]'''//.basicRules//Expand*)


(* ::Input:: *)
(*DBT[DT[null]]//.basicRules//Expand*)


(* ::Input:: *)
(*moreRules={NO[DT[J],DT[W]]-> 0,*)
(*NO[DT[J''],DT[W]]-> NO[DT[J''],DT[W]]-4/35 DT[null]'',*)
(*NO[DT[J'''],DT[W]]-> NO[DT[J'''],DT[W]]-12/119 DT[null]''',*)
(*NO[DT[J],W']-> NO[DT[J],W']+DBT[DT[null]]/7}//.basicRules//Expand;*)


(* ::Text:: *)
(*At this point the failure of the Jacobi are not very long:*)


(* ::Input:: *)
(*remainder=(OPEJacobi@@#)&/@Permutations[{W,W,WB}]//.Join[basicRules,moreRules]//Factor//Flatten//Prune//Union*)


(* ::Text:: *)
(*Finally, we need a few really ad hoc descendants to see that  the remainders above are set to zero:*)


(* ::Input:: *)
(*A=DBT[DT[null]'']//.basicRules//Simplify//Expand*)


(* ::Input:: *)
(*B=NO[J',DBT[DT[null]]]//.basicRules//Simplify//Expand*)


(* ::Input:: *)
(*evenMoreRules=Solve[{A==0,B==0},{NO[DT[J^\[Prime]\[Prime]],Derivative[1][W]], NO[\!\(\*SuperscriptBox[\(J\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\),DT[W]]}][[1]]//.basicRules//Expand*)


(* ::Input:: *)
(*S=NO[DDB[J],DT[null]']//.Join[basicRules,moreRules,evenMoreRules]//Expand*)


(* ::Input:: *)
(*T=NO[DT[DDB[J]],null']//.Join[basicRules,moreRules,evenMoreRules]//Expand*)


(* ::Input:: *)
(*furtherRules=Solve[{S==0,T==0},{ NO[DT[J],NO[DT[Derivative[1][J]],Derivative[1][W]]], NO[DT[J],NO[DDB[Derivative[1][J]],DT[W]]]}][[1]]*)


(* ::Input:: *)
(*DBT[B]12/7//.Join[basicRules,moreRules,evenMoreRules,furtherRules]//Expand*)


(* ::Input:: *)
(*finalRules={NO[\!\(\*SuperscriptBox[\(J\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\),Derivative[1][W]]-> NO[\!\(\*SuperscriptBox[\(J\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\),Derivative[1][W]]+DBT[B]12/7}//.Join[basicRules,moreRules,evenMoreRules,furtherRules]//Expand*)


(* ::Input:: *)
(*remainder//.Join[basicRules,moreRules,evenMoreRules,furtherRules,finalRules]//Factor//Flatten//Union*)


(* ::Text:: *)
(*Note also that the OPE of the null field with J, W, WB contains operators with higher dimensions only.*)


(* ::Input:: *)
(*OPE[J,null]-primaryope[null]//Simplify*)


(* ::Input:: *)
(*OPE[W,null]//Simplify*)


(* ::Text:: *)
(*The OPE with WB is very long, so we do not show it explicitly*)


(* ::Input:: *)
(*tmp=OPE[WB,null]//Simplify;*)


(* ::Text:: *)
(*But we can check that its leading pole is of order 2,*)


(* ::Input:: *)
(*MaxPole[tmp]*)


(* ::Text:: *)
(*which guarantees that indeed there is no term with higher dimension than the null operator itself.*)


(* ::Text:: *)
(*In fact the residue of this second order pole does not have a term without \[Theta] or \!\(\*OverscriptBox[\(\[Theta]\), \(_\)]\):*)


(* ::Input:: *)
(*coeffOf1[OPEPole[2][tmp]]*)


(* ::Text:: *)
(*The leading term in the OPE is the following WB descendant of the null field and it gives the Higgs branch relation:*)


(* ::Input:: *)
(*relation=coeffOftb12[OPEPole[2][tmp]]//Expand*)


(* ::Text:: *)
(*So the relation is J^2/5040 =36W \!\(\*OverscriptBox[\(W\), \(_\)]\)  up to descendants.*)


(* ::Section:: *)
(*Anti-climactic aside: k=5*)


(* ::Text:: *)
(*As a final exercise, let us show that the 2d chiral algebra for the would-be 4d N=3 k=5 theory also exists.*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Input:: *)
(*\[CapitalDelta][W]=5/2;F[W]=5;*)


(* ::Input:: *)
(*OPE[J,W]=primaryope[W]//OPESimplify*)


(* ::Input:: *)
(*\[CapitalDelta][WB]=5/2;F[WB]=-5;*)


(* ::Input:: *)
(*OPE[J,WB]=primaryope[WB]//OPESimplify*)


(* ::Input:: *)
(*OPE[J,J]=MakeOPE[{c/3 One + t12**tb12**J,tb12**DBT[J]-t12**DT[J]+t12**tb12**J'}]*)


(* ::Text:: *)
(*The OPE of W and WB has the general form*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(5/2)a t12**tb12**One,*)
(*a One+(4/2)t12**tb12**J,*)
(*J +  t12**DT[J]*)
(*+(3/2)t12**tb12**gL[2,A],*)
(*gL[2,A] +  t12**  DT[gL[2,A]]*)
(*+(2/2)t12**tb12**gL[3,AA],*)
(*gL[3,AA]+t12**DT[gL[3,AA]]*)
(*+(1/2) t12**tb12**gL[4,AAA],*)
(*gL[4,AAA]+t12**DT[gL[4,AAA]]}]//OPESimplify;*)


(* ::Text:: *)
(*from which we can fix the coefficients by demanding the J-W-WB Jacobi. Then we study the Jacobi for W, W, WB:*)


(* ::Input:: *)
(*OPEJacobi[J,W,WB,Function->Factor]//Union;*)


(* ::Input:: *)
(*solution=Solve[GetCoefficients[%]==0];*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]/.solution[[1]]//Factor//Prune//Union;*)


(* ::Input:: *)
(*{remainder[[2]],remainder[[3]]}*)


(* ::Input:: *)
(*PolynomialGCD[-1008-99 c+c^2,33768+37395 c-1469 c^2-95 c^3+9 c^4]*)


(* ::Text:: *)
(*We see that only c=\[Minus]27 is allowed at this point. Note that \[Minus]27=\[Minus]3(2k\[Minus]1). Let us redefine the basic OPEs with the determined values:*)


(* ::Input:: *)
(*ClearOPESavedValues[]*)


(* ::Input:: *)
(*OPE[J,J]=MakeOPE[{c/3 One + t12**tb12**J,tb12**DBT[J]-t12**DT[J]+t12**tb12**J'}]/.c-> -27;*)


(* ::Input:: *)
(*OPE[W,WB]=MakeOPE[{(5/2)a t12**tb12**One,*)
(*a One+(4/2)t12**tb12**J,*)
(*J +  t12**DT[J]*)
(*+(3/2)t12**tb12**gL[2,A],*)
(*gL[2,A] +  t12**  DT[gL[2,A]]*)
(*+(2/2)t12**tb12**gL[3,AA],*)
(*gL[3,AA]+t12**DT[gL[3,AA]]*)
(*+(1/2) t12**tb12**gL[4,AAA],*)
(*gL[4,AAA]+t12**DT[gL[4,AAA]]}]/.solution[[1]]/.c-> -27//OPESimplify;*)


(* ::Input:: *)
(*remainder=OPEJacobi[W,W,WB,Function->Factor]//Factor//Prune//Union;*)


(* ::Input:: *)
(*remainder[[2]]*)


(* ::Text:: *)
(*So we need the following null field:*)


(* ::Input:: *)
(*null=DT[W']-NO[J,DT[W]]+5NO[DT[J],W];\[CapitalDelta][null]=4;F[null]=4;*)


(* ::Text:: *)
(*Necessaray descendants can be constructed completely analogously as was in k=6:*)


(* ::Input:: *)
(*DBT[null]//Expand*)


(* ::Input:: *)
(*basicRules={DT[W']-> DT[W']-null,*)
(*W''-> W''+DBT[null],W'''-> W'''+DBT[null]',W''''-> W''''+DBT[null]'',W'''''-> W'''''+DBT[null]'''}//Expand;*)


(* ::Input:: *)
(*DT[null]*)


(* ::Input:: *)
(*DT[null]'//Expand*)


(* ::Input:: *)
(*DT[null]''//.basicRules//Expand*)


(* ::Input:: *)
(*DT[null]'''//.basicRules//Expand*)


(* ::Input:: *)
(*DBT[DT[null]]//.basicRules//Expand*)


(* ::Input:: *)
(*moreRules={NO[DT[J],DT[W]]-> 0,*)
(*NO[DT[J''],DT[W]]-> NO[DT[J''],DT[W]]-1/6 DT[null]'',*)
(*NO[DT[J'''],DT[W]]-> NO[DT[J'''],DT[W]]-1/7 DT[null]''',*)
(*NO[DT[J],W']-> NO[DT[J],W']+DBT[DT[null]]/6}//.basicRules//Expand;*)


(* ::Text:: *)
(*We see that the Jacobi then vanishes:*)


(* ::Input:: *)
(*remainder=(OPEJacobi@@#)&/@Permutations[{W,W,WB}]//.Join[basicRules,moreRules]//Factor//Flatten//Prune//Union*)


(* ::Text:: *)
(*Note also that the OPE of the null field with J, W, WB contains operators with higher dimensions only.*)


(* ::Input:: *)
(*OPE[J,null]-primaryope[null]//Simplify*)


(* ::Input:: *)
(*OPE[W,null]//Simplify*)


(* ::Input:: *)
(*tmp=OPE[WB,null]//Simplify;*)


(* ::Input:: *)
(*MaxPole[tmp]*)


(* ::Text:: *)
(*The residue of this second order pole does not have a term without \[Theta] or \!\(\*OverscriptBox[\(\[Theta]\), \(_\)]\):*)


(* ::Input:: *)
(*coeffOf1[OPEPole[2][tmp]]*)


(* ::Text:: *)
(*The leading term in the OPE is the following WB descendant of the null field and and it gives the Higgs branch relation:*)


(* ::Input:: *)
(*relation=coeffOftb12[OPEPole[2][tmp]]//Expand*)


(* ::Text:: *)
(*So the relation is J^2/336 =25W \!\(\*OverscriptBox[\(W\), \(_\)]\)  up to descendants.*)


(* ::Text:: *)
(*It seems that for general k we still have a corresponding 2d algebra at c=\[Minus]3(2k-1) and the relation is *)
(*               n!/(2(n-1))!J^2 =k^2W \!\(\*OverscriptBox[\(W\), \(_\)]\)  *)
(*up to descendants.*)
