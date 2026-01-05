(*^

::[	frontEndVersion = "Microsoft Windows Mathematica Notebook Front End Version 2.2";
	microsoftWindowsStandardFontEncoding;
	fontset = title, "Times New Roman", 24, L0, center, nohscroll, bold;
	fontset = subtitle, "Times New Roman", 18, L0, center, nohscroll, bold;
	fontset = subsubtitle, "Times New Roman", 14, L0, center, nohscroll, bold;
	fontset = section, "Times New Roman", 14, L0, bold, grayBox;
	fontset = subsection, "Times New Roman", 12, L0, bold, blackBox;
	fontset = subsubsection, "Times New Roman", 10, L0, bold, whiteBox;
	fontset = text, "Times New Roman", 12, L0;
	fontset = smalltext, "Times New Roman", 10, L0;
	fontset = input, "Courier New", 11, L0, nowordwrap;
	fontset = output, "Courier New", 11, L0, nowordwrap;
	fontset = message, "Courier New", 10, L0, nowordwrap, R65535;
	fontset = print, "Courier New", 10, L0, nowordwrap;
	fontset = info, "Courier New", 10, L0, nowordwrap;
	fontset = postscript, "Courier New", 8, L0, nowordwrap;
	fontset = name, "Times New Roman", 10, L0, nohscroll, italic, B65535;
	fontset = header, "Times New Roman", 10, L0, right, nohscroll;
	fontset = footer, "Times New Roman", 10, L0, right, nohscroll;
	fontset = help, "Times New Roman", 10, L0, nohscroll;
	fontset = clipboard, "Times New Roman", 12, L0, nohscroll;
	fontset = completions, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = graphics, "Courier New", 10, L0, nowordwrap, nohscroll;
	fontset = special1, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = special2, "Times New Roman", 12, L0, center, nowordwrap, nohscroll;
	fontset = special3, "Times New Roman", 12, L0, right, nowordwrap, nohscroll;
	fontset = special4, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = special5, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = leftheader, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = leftfooter, "Times New Roman", 12, L0, nowordwrap, nohscroll;
	fontset = reserved1, "Courier New", 10, L0, nowordwrap, nohscroll;]
:[font = text; inactive; preserveAspect; ]
(* :Name: SOPEN2defs` *)

(* :Title: 
  	N=2 Superfield OPEs
*)

(* :Author: 
    Sergey Krivonos and Kris Thielemans
 *)

(* :Summary:
   This package computes Operator Product Expansions of composite
   operators in meromorphic N=2 conformal field theory. All calculations
   are done in superfields. An interface is provided to construct the 
   components of the superfields (and their OPEs), such that results can 
   be checked with OPEdefs ( from KT).
*)

(* :Context: SOPEN2defs` *) 

(* :Package Version: 1.0 *)

(* :Mathematica Version: 2.0 *)
:[font = text; inactive; preserveAspect; ]
This package is written by
     Kris Thielemans
    Theoretical Physics Group
    Imperial College
    London (UK)
    Email-address: kris@tfdec1.fys.kuleuven.ac.be
and
    Sergey Krivonos
    Joint Institute for Nuclear Research
    Dubna (Russia)
    Email-address : krivonos@thsun1.jinr.dubna.su

For details, see our paper:
	A Mathematica package for Computing 
	N=2 Superfield Operator Product Expansions
Preprint Imperial/TP/95-96/13, hep-th/9512029.
Please refer to it in your publications if you used this package in 
your work. See the History section below if you used an older version
of the package before.

You're free to redistribute this package, on the only condition you
keep this header and don't distribute modified code.

Of course, we do not guarantee you get correct results (even if you use
the package correctly !), although it is tested rather extensively.

Please contact us if you want to use the package, such that we know
who uses it, and you can get any updates or related packages.  Any
comments, bugs and especially improvements are welcome.

This package modifies the definitions of Derivative and NonCommutativeMultiply.
:[font = section; inactive; preserveAspect; startGroup; Cclosed; ]
History
:[font = text; inactive; preserveAspect; endGroup; ]
first version 0.5
   based on OPEdefs 3.0
version 1.0 : (updated to OPEdefs 3.1 standards)
beta 1:
   - Replaced GeneralBosonic, GeneralFermionic with simply Bosonic, 
      Fermionic
   - (A)ChiralQ now global. Useful to be able to set it by hand (and not
      through (A)Chiral*onic).
   - Added OPEJacobi
   - N2OPEToComponents
      added factor -1/2 for the last component such that  for a superfield 
           A + theta B + thetabar C + theta thetabar D
      the components are {A,B,C,D}
   - MaxPole is now a global function.
   - GetCoefficients is now Listable. That means when acting on a list of
     operators, it will returns a list of lists. To get the same result as before,
     use Flatten.
   - Addition GetOperators.
   - OPESimplify now accepts the option Function.
   - Addition OPEParity, OPEOperator and the option ParityMethod
   - internal function SwapSign now also accepts only 1 argument
   - CallAndSave now set via SetOPEOptions. slightly faster.
   - TeXForm now works on SOPEN2defs structures. See also TeXFormDT
   - Import package Delta instead of having the definition in here.
   - Added support for Dummies` via the option EnableDummies.
   beta 2:
   - Write warning message when loading Dummies` concerning OPESaving and call      ResetDummies.

   !!!Taylor  !!!
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
BeginPackage["SOPEN2defs`", "Delta`"];
(*
:[font = section; inactive; preserveAspect; startGroup; Cclosed; ]
Exported symbols
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
SOPEN2defsVersion = 1.0;
betaVersion = 2;
Print["SOPEN2defs Version ", SOPEN2defsVersion,
      " (beta ",betaVersion,") by Kris Thielemans and Sergey Krivonos"]
Print["Type ?SOPEN2defsHelp for a primer on SOPEN2defs (try it !)."]
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
general help
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

SOPEN2defsHelp::usage = StringJoin[
"First, declare your operators using Bosonic, Fermionic,
ChiralBosonic, AChiralBosonic, ...\n
    e.g.: Bosonic[T, J[_]]\n",
"Then, define the non-zero OPEs using MakeOPE. If you declared operator
A before B, you HAVE TO define OPE[A,B] (not OPE[B,A]).\n
    e.g.: OPE[T,T] = MakeOPE[{c/4 One + t12**tb12**T,
                        tb12**DBT[T]-t12**DT[T]+t12**tb12**T'}];\n",
"Compute what you want using OPE, NO (normal ordered products),
OPEPole, OPESimplify\n",
"Note : if you want to do Poisson bracket computations, issue\n
          SetOPEOptions[OPEMethod, ClassicalOPEs]\n
before any other definitions.\n
\nType ?SOPEN2defs`* for all global symbols.
\nReloading the package (with Get or \"<<\"), clears all your definitions
and stored intermediate results concerning OPEs.\n
Type ?Basis for the definiton of the basis in superspace we use."]

Basis::usage =
"We use the following notation :\n
     DT is the spinor theta-derivative\n
     DBT is the spinor theta-bar-derivative\n
     DDB[A] =  DT[DBT[A]]-DBT[DT[A]]\n
     t12, tb12 are the theta's.\n
We use the following basis (with DT,DBT acting on theta2, z2):\n
     DT[t12] = -1,          DBT[tb12] = -1\n
     DT[z12] = -1/2 tb12,   DBT[z12] = -1/2 t12\n
     {DT,DBT} = - Derivative[1]\n
or in an obious notation :\n
     DT  = d_theta - 1/2 thetabar d_z\n
     DBT = d_thetabar - 1/2 theta d_z."
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
declarations
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
GeneralBosonic::usage =
  "GeneralBosonic is obsolete. Use Bosonic instead."
GeneralFermionic::usage =
  "GeneralFermionic is obsolete. Use Fermionic instead."

Bosonic::usage =
"Bosonic[A] declares A as a general bosonic N=2 SF.
You can use Bosonic[A,B,...]
(see also (A)ChiralBosonic,Fermionic,(A)ChiralFermionic,
OPEOperator, (A)ChiralQ)"

ChiralBosonic::usage =
"ChiralBosonic[A] declares A as a chiral bosonic N=2 SF.
You can use ChiralBosonic[A,B,...]
(see also Bosonic,AChiralBosonic,Fermionic,(A)ChiralFermionic,
OPEOperator, (A)ChiralQ)"

AChiralBosonic::usage =
"AChiralBosonic[A] declares A as a anti-chiral bosonic N=2 SF.
You can use AChiralBosonic[A,B,...]
(see also Bosonic,ChiralBosonic,Fermionic,(A)ChiralFermionic,
OPEOperator, (A)ChiralQ)"

Fermionic::usage =
"Fermionic[A] declares A as a general fermionic N=2 SF.
You can use Fermionic[A,B,...]
(see also (A)ChiralFermionic,Bosonic,(A)ChiralBosonic,
OPEOperator, (A)ChiralQ)"

ChiralFermionic::usage =
"ChiralFermionic[A] declares A as a chiral fermionic N=2 SF.
You can use ChiralFermionic[A,B,...]
(see also Bosonic,(A)ChiralBosonic,Fermionic,AChiralFermionic,
OPEOperator, (A)ChiralQ)"

AChiralFermionic::usage =
"AChiralFermionic[A] declares A as a anti-chiral fermionic N=2 SF.
You can use AChiralFermionic[A,B,...]
(see also Bosonic,(A)ChiralBosonic,Fermionic,ChiralFermionic,
OPEOperator, (A)ChiralQ)"

OPEOperator::usage =
"OPEOperator[A,p] declares A as an operator of parity p. Here, p can
be a symbolic expression. You can use OPEOperator[{A,pA},{B,pB},...]. Using
only Bosonic and Fermionic gives slightly faster calculations
(see also Bosonic,Fermionic,OPEParity,ParityMethod,(A)ChiralQ).\n
e.g.: OPEOperator[J[i_], parity[i]]"

(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
properties
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPEParity::usage =
"OPEParity[A] returns the parity of the operator A. The parity of a 
bosonic (fermionic) operator is an even (odd) integer number."

ChiralQ::usage =
"ChiralQ[A] tests if A is a chiral operator.\n
One can forse this to be true (implying DT[A]=0) by typing\n
\tChiralQ[A] = True\n
(see also AChiralQ, OPEOperator, (A)ChiralBosonic, (A)ChiralFermionic)"

AChiralQ::usage =
"AChiralQ[A] tests if A is an antichiral operator.\n
One can forse this to be true (implying DBT[A]=0) by typing\n
\tAChiralQ[A] = True\n
(See also ChiralQ, OPEOperator, (A)ChiralBosonic, (A)ChiralFermionic)"
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPE, NO ...
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPEData::usage =
        "The head of the structure returned by OPE. For future
compatibility, do not use its structure directly. See MakeOPE instead."

MakeOPE::usage =
        "Must be used to make an OPEData structure. The argument is
a list of the operators occuring in the poles of this series.
Begin with the highest order pole, and list ALL poles up to the first order
pole (i.e. also the zero ones)\n
    e.g.: MakeOPE[{c/4 One + t12**tb12**T,
                   tb12**DBT[T]-t12**DT[T]+t12**tb12**T'}]\n
Note that all terms should include an operator, use the constant
operator One for the unit operator (e.g. t12**One)."

OPE::usage =
        "OPE[A,B] gives the operator product expansion of
A[z1,theta1] B[z2,theta2] in the form of an OPEData structure which
lists the poles (with the argument z2 dropped) which occur in the
series in z12 up to order 0 in z12.\n
(see also MakeOPE)"

OPEPole::usage =
        "OPEPole[i][A,B] gives the operator occuring at the 'i'-th
order pole in the OPE of A[z] with B[w]. OPEPole[0][A,B] gives
NO[A,B] + t12**NO[DT[A],B] + tb12**NO[DBT[A],B] - 1/2 NO[DDB[A],B].
Negative 'i' can be used to get expressions for the regular part
of the OPE.\n
e.g.: (assuming T is defined as a super Virasoro operator)\n
   OPEPole[1][T,T]\n
   -->  tb12**DBT[T]-t12**DT[T]+t12**tb12**T'\n\n
OPEPole can also be used as OPEPole[i][some_ope]. In this case,
'i' must be strictly positive (otherwise, OPEPole returns 0).\n
e.g.: tt = OPE[T,T]; OPEPole[2][tt]\n
   --> c/4 One + t12**tb12**T"

NO::usage =
        "NO[A,B] denotes the normal ordered product of
A and B (point splitting convention). NO[A,B] is defined as the
theta-independent part of OPEPole[0][A,B]. Simplifications are done
to reduce NO's to a standard form : normal ordering is from the right
to the left, earlier declared operators are more to the left, and
higher derivatives are (by default) more to the left than lower ones.\n
(see also NOOrdering)."

One::usage =
        "The constant operator : derivatives of it are considered
zero, an OPE with anything else is regular, and the NO with any other
operator is that operator."

OPEMap::usage =
        "OPEMap[f, ope] maps a function 'f' to all poles of 'ope'.
OPEMap[f, ope, level] maps 'f' at a specific level (see Map).\n
e.g.: OPEMap[Expand, ope]\n
(see also OPEMapAt)."

OPEMapAt::usage =
        "OPEMapAt[f, ope, position] maps a function 'f' to a certain
point of 'ope' (see MapAt).\n
e.g.: OPEMapAt[Expand, ope, {1}]\n
  maps Expand on the highest order pole of 'ope'\n
(see also OPEMap)."

OPESimplify::usage =
        "Applied to a sum of operators, collects terms
with the same operator.\n 
Applied to an OPE collects terms with the same operator in every pole.
It effectively maps OPESimplify on all the poles of the OPE.\n
OPESimplify distributes over lists.\n
OPESimplify accepts the option Function which determines which 
function it will apply on the coefficients. By default
this is set to Expand. Use SetOptions to change the default.\n
e.g.: OPESimplify[T'' + 1/c T'' + 2 NO[T,T] + c NO[T,T] , 
           Function -> Together]\n
   --> (1+c)/c T'' + (2+c) NO[T,T]\n
The effect of the option Function
is also achieved by using OPESimplify[expr, function]."

GetCoefficients::usage =
        "GetCoefficients[expr] returns a list of the coefficients of all
operators in expr. The order in which they are given is currently
undefined.\n
(see also GetOperators)."

GetOperators::usage = 
"GetOperators[expr] returns a list of operators in expr.\n
(see also GetCoefficients)."

MaxPole::usage = "MaxPole[ope] gives the order of the highest pole."
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
saved values
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPESave::usage =
        "OPEdefs by default remembers some intermediate results
(in particular OPEs of composites, and reordering of some composites).
OPESave can be used to save these results in a file. In a next session,
you can read that file (after loading OPEdefs AND your declarations of
the operators and OPEs you use).\n
Note : do not use intermediate results saved by an older version
of OPEdefs.\n
e.g.: OPESave[\"results.m\"]\n
      Quit[]\n
      ... restart mathematica\n
      Needs[\"OPEdefs`\"]\n
      ... read your normal definitions\n
      <<results.m\n
See also OPESaving."

ClearOPESavedValues::usage = "ClearOPESavedValues[] clears all the
intermediate results that were stored. It doesn't clear the definitions
for the OPEs."
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPEJacobi
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPEJacobi::usage = "OPEJacobi[A,B,C] check the Jacobi identities
for A,B,C. It returns a double list of expressions which should
all be zero (or a null field). In general, all different orders
for the operators should be tried."
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Components
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

N2OPEToComponents::usage =
        "N2OPEToComponents[J] gives the N=0 components of a
superfield J. For J= A + theta B + thetabar C + theta thetabar D
the components are {A,B,C,D}. For a 'basic' superfield, the
components are represented as {J,DT[J],DBT[J], -1/2 DDB[J]}.
After projection on the theta independent part, this gives the
correct result.\n
N2OPEToComponents[ope,A,B] computes the 16 OPEs of the N=0
components of A and B. It is a double list composed as\n
Table[OPE[comp[i][A],comp[j][B]], {i,4},{j,4}]\n
where comp[i][A] represents N2OPEToComponents[A][[i]]."

Null
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Options
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)

SetOPEOptions::usage =
        "Function to set some options of the OPE package.
Currently available options : OPESaving, NOOrdering,
ParityMethod, OPEMethod, EnableDummies."

OPESaving::usage =
        "OPEdefs by default remembers some intermediate results
(in particular OPEs of composites, and reordering of some composites).
This of course speeds up future computations, but can require lots of
memory, which results in most of your time is spent in swapping parts of
memory to disk. A switch is provided to prevent OPEdefs to remember
these results.\n
The second argument of SetOPEOptions has to be an expression that will
evaluate to a boolean value (True of False).\n
e.g.: SetOPEOptions[OPESaving, True]\n
Setting the switch to False prevents saving of results in subsequent
calculations. The switch defaults to True.\n
A more sophisticated use is the following :\n
          SetOPEOptions[OPESaving, MaxMemoryUsed[] < 6 10^6]\n
This will save intermediate results as long as you don't have used
more than 6 Mb. (On most systems using MaxMemoryUsed[] is more appropriate
than MemoryInUse[] to avoid swapping.)"

NOOrdering::usage =
        "An option to set the order of derivatives of
fields inside normal ordered products. Specify an integer number (n) as
the second argument of SetOPEOptions. Three cases apply :\n
    n <0 : higher derivatives of a field are moved to the left (default)\n
    n==0 : no reordering (not advised, certainly when OPESaving is set
to True)\n
    n >0 : lower derivatives of a field are moved to the left.\n
Please note that resetting the switch after calculations are done (while
OPESaving is set to True) causes some intermediate results to be ignored.\n
e.g. : SetOPEOptions[NOOrdering, -1]\n
       NO[J, J']  -> NO[J', J] + ..."

OPEMethod::usage = "An option to set the way OPEs are computed.
Currently implemented methods : QuantumOPEs, ClassicalOPEs
(see e.g. ?QuantumOPEs). The first method is the default.\n
e.g.: SetOPEOptions[OPEMethod, ClassicalOPEs]."

ClassicalOPEs::usage = "A value for the OPEMethod option that
specifies that OPEs will be computed as Poisson brackets (no
multiple contractions, graded-commutative and associative normal
ordering).\n
Important : if OPESaving is set to true (the default), this option
can currently only be safely set before any OPEs are used.\n
See also OPEMethod, QuantumOPEs, OPESaving."

QuantumOPEs::usage =  "A value for the OPEMethod option that
specifies that OPEs will be computed as normal OPEs (multiple
contractions, normal ordering).\n
Important : if OPESaving is set to true (the default), this option
can currently only be safely set before any OPEs are used.\n
See also OPEMethod, ClassicalOPEs, OPESaving."

ParityMethod::usage = "An option that determines the way OPEdefs computes
signs.
When the second argument to SetOPEOptions is 0, all operators have to be
declared to be Bosonic or Fermionic (this is the default). When the argument
is 1, symbolic parities can be used.\n
Note that in the last case, powers of $-1$ are used to compute signs, which is
slightly slower than the boolean function which is used by the first
method.\n
This option is not normally needed as using OPEOperator with a 
symbolic second argument sets this option automatically.\n
e.g.: SetOPEOptions[ParityMethod,1]"
(*
:[font = input; initialization; endGroup; nowordwrap; ]
*)
EnableDummies::usage = "An option which enables use of the Dummies` package. 
Use\n
\tSetOPEOptions[EnableDummies, True]"
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
theta's
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
DT::usage =
"DT[A] gives the spinor theta derivative of A"

DBT::usage =
"DBT[A] gives the spinor theta-bar derivative of A"

DDB::usage =
"DDB[A] gives DT[DBT[A]]-DBT[DT[A]]"

t12::usage = "t12 is definition of theta_12"

tb12::usage = "tb12 is definition of theta-bar_12"
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Delta
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
(*Delta::usage =
        "Delta[i,j] is the delta symbol."*)
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Mathematica 1.2 additions
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)

If [!NumberQ[$VersionNumber] || $VersionNumber < 2.0,
StringReplace::usage =
        "Own implementation of a Mathematica Version 2.0 function.";
Fold::usage =
        "Own implementation of a Mathematica Version 2.0 function."
]
Null
(*
:[font = section; inactive; preserveAspect; startGroup; ]
Implementation
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
Begin["`Private`"];
(*
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
Nsusy=2;
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Auxiliary functions (Delta, StringReplace, Fold)
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
(* change 1.0 beta 1 : Delta moved to separate package
Clear[Delta]
Delta[a_,a_] = 1
Delta[i_,j_] := 0 /; i!=j
SetAttributes[Delta, Orderless]
*)
 (********************* StringReplace and Fold **********************)
If[!NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* 15-10-92 (3.0 beta 5)
       previous implementation changed only first occurences, unlike
       StringReplace in Mathematica 2.0 -> Use FixedPoint
    *)
    StringReplace[str_String, (str1_String)->(str2_String)] :=
        StringReplace[str,{str1->str2}];
    StringReplace[str_String, rules_List] := Block[
        { newrules =
             rules /.
                  (  (str1_String)->(str2_String) )
                  :>
                  ( {p1___,Sequence @@ Characters[str1],p2___}
                     ->{p1,{str2},p2}
                  ) ,
          chars = Characters[str]
        },
        (* Replace changes only toplevel expressions *)
        chars = FixedPoint[Replace[#,newrules]&, chars];
        StringJoin[Sequence @@ Flatten[chars]]
        ];

     Fold[func_, expr_, l_List] := Block[
         { res = expr
         },
         Scan[ (res = func[res, #])&, l ];
         res
     ]
]
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Rules for Noncommutative Multiplication
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Unprotect[NonCommutativeMultiply];
ClearAll[NonCommutativeMultiply]

Literal[NonCommutativeMultiply[a___,NonCommutativeMultiply[b__],c___]] :=
  NonCommutativeMultiply[a,b,c]

Literal[NonCommutativeMultiply[a_]] := a

NonCommutativeMultiply[] := 1

If[NumberQ[$VersionNumber],
    Literal[NonCommutativeMultiply[a___,b_Plus,c___]]:=
            Distribute[
                Lineartmp[a,b,c],
                Plus,Lineartmp,
                Plus,NonCommutativeMultiply
        ],
    (* bug in version 1.1 Distribute *)
    Literal[NonCommutativeMultiply[a___,b_+c_,d___]] :=
        NonCommutativeMultiply[a,b,d] +
        NonCommutativeMultiply[a,c,d]
]

Literal[NonCommutativeMultiply[___,0,___]]=0;

 (* rule necessary because we use thetas12[0]=1 further on *)
Literal[NonCommutativeMultiply[a___,1,b___]]:=
        NonCommutativeMultiply[a,b];

 (* rule removed as it slows down calculations almost twice *)
(*Literal[NonCommutativeMultiply[a___,b_,e___]] :=
  b NonCommutativeMultiply[a,e] /; !OperatorQ[b]
 *)

Literal[NonCommutativeMultiply[a___,s_ c_,e___]] :=
  s NonCommutativeMultiply[a,c,e] /; OperatorQ[c]

 (* note: we can have only theta's more than once,
    these are always fermionic *)
Literal[NonCommutativeMultiply[___,b_,___,b_,___]] :=
         0 (* /; !BosonQ[b]*)

Literal[NonCommutativeMultiply[c___,b_,a_,e___]] :=
SwapSign[a,b] NonCommutativeMultiply[c,a,b,e] /; ord[a,b] >0

Clear[ord]
ord[tb12,t12]=-1;
ord[t12,tb12]=1;
ord[_,_]=0;
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Functions to manipulate OPEdatas
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPEData
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[OPEData]
OPEData /: n_ * OPEData[A_List] := OPEData[n* A]

extend[A_List,n_Integer] :=
        Join[Table[0,{n-Length[A]}], A]
 (* 01-11-92 (3.0 beta 6)
  Previous implementation :
   OPEData /: Plus[A__OPEData] :=
        Block[{maxP = Max[MaxPole /@ {A}]},
            OPEData[ Plus @@ (extend[#[[1]],maxP]& /@ {A}) ]
        ]
  When given Plus[OPEData[a1], OPEData[a2]], this rule applies to give
  Plus[OPEData[a1+a2]]. Then the rule matches again (and does nothing).
  Also, when entering an (invalid) sum of an OPEData with something else,
  this rule causes infinite recursion.
  Change : make sure the rule matches only when there are at least 2
           OPEDatas. Also use a pure-function instead of "extend", this
           turns out to be more efficient.
  Algorithm : use First to get the list of the poles for each OPE.
           Then "extend" the lists with zeroes such that they have equal length.
           Add these lists, and bracket with OPEData.
*)
A1_OPEData + A2__OPEData ^:=
        Block[{maxP = Max[MaxPole /@ {A1,A2}]},
            OPEData[ Plus @@
                (Join[ Table[0,{maxP-Length[#]}], # ]& /@
                    (First /@ {A1,A2})
                )
            ]
        ]
 (* 13-05-93 (3.0 beta 8) added .., although it doesn't change
    anything upto Mathematica 2.1 (and further ?)
 *)
OPEData[{(0).., A___}] := OPEData[{A}]
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
MakeOPE
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[MakeOPE]
Ord[z_,w_,i_Integer] := SeriesData[z,w,{},i,i,1]
SeriesDataToOPEData =
        Literal[SeriesData[z_,w_,A_List,highest_Integer,0,1]] :>
        OPEData[Join[A /. (OP_[w] -> OP),
                                 Table[0, {-highest-Length[A]}]]];
checkOne[ope_OPEData] :=
     OPEMap[If[OperatorQ[#], #, # One]& , ope]

MakeOPE[l_List] := OPEData[l]

MakeOPE[l_SeriesData] := checkOne[(l /. SeriesDataToOPEData)]

MakeOPE[ope_OPEData] := ope

MakeOPE[l_] := (Message[MakeOPE::other, l]; MakeOPE[{}])
MakeOPE::other = "Error : `1` is not a SeriesData or a List object. Assuming
it is equivalent to MakeOPE[{}].";
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
simplification
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[OperatorPattern, PoleSimplify, OPESimplify]

OperatorPattern := Alternatives @@ OperatorList

PoleSimplify[term_] :=
   Block[{expterm = Expand[term], var},
        var = ExtractOperators[expterm];
        Apply[Plus, (# Coefficient[expterm, #])& /@ var]
   ]
(* change 1.0 beta 1 : added option *)
PoleSimplify[term_, Function -> func_] :=
   PoleSimplify[term, func] 

PoleSimplify[term_, func_] :=
   Block[{expterm = Expand[term], var},
        var = ExtractOperators[expterm];
        Apply[Plus, (# func[Coefficient[expterm, #]])& /@ var]
   ]

Clear[ExtractOperators]
 (* ExtractOperators[expr_] returns a list of operators in expr
    !! It assumes that expr is Expanded
 *)
ExtractOperators[a_Plus] := Union[ ExtractOperators[#][[1]]& /@ (List @@ a)]
If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* old version of Mathematica *)
    ExtractOperators[a_Times] := Select[List@@a, OperatorQ]
    ,
    (* new version : use Select[ , ,1] to extract the first operator
       you find (there is only one !)
    *)
    ExtractOperators[a_Times] := Select[List@@a, OperatorQ, 1]
]
ExtractOperators[0] = {};
ExtractOperators[a_] := {a}

Clear[ExtractOperatorsNoTheta]
 (* ExtractOperatorsNoTheta[expr_] returns a list of operators in expr
    !! It assumes that expr is Expanded
 *)
ExtractOperatorsNoTheta[a_Plus] := Union[ ExtractOperatorsNoTheta[#][[1]]& /@ (List @@ a)]
If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* old version of Mathematica *)
    ExtractOperatorsNoTheta[a_Times] :=
        Select[List@@a/.{t12->1,tb12->1}, OperatorQ]
    ,
    (* new version : use Select[ , ,1] to extract the first operator
       you find (there is only one !)
    *)
    ExtractOperatorsNoTheta[a_Times] :=
        Select[List@@a/.{t12->1,tb12->1}, OperatorQ, 1]
]
ExtractOperatorsNoTheta[0] = {};
ExtractOperatorsNoTheta[a_] := {a/.{t12->1,tb12->1}}

(* change 1.0 beta 1: added options *)
Options[OPESimplify] = {Function -> Expand};

OPESimplify[A_OPEData, func___] :=
        OPEMap[PoleSimplify[#,func]& , A]
 (* 13-05-93 (3.0 beta 8) added two following definitions *)
OPESimplify[A_List, func___] :=
        OPESimplify[#,func]& /@ A
OPESimplify[A_, args___] :=
        PoleSimplify[A, args]

If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* bug in old versions : OPEData[{a}] - OPEData[{a}] -> 0 *)
    OPESimplify[0, ___] = 0;
    OPEMap[_, 0,___] = 0;
    OPEMapAt[_, 0,_] = 0;
]
(*
:[font = subsubsection; inactive; startGroup; Cclosed; ]
Mapping
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; backColorRed = 65535; backColorGreen = 65535; backColorBlue = 65535; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; plain; fontName = "Courier New"; fontSize = 11; ]
*)
Clear[OPEMap, OPEMapAt]

OPEMap[f_, ope_OPEData, levelspec___] :=
    MapAt[Map[f, #, levelspec] & , ope, {1}]
 (*
  OPEMap[f_, ope_OPEData] :=
        MapAt[Map[If[Head[#]=!=Plus, f[#], (f /@ #1)] &,#]& , ope, {1}]
 *)
OPEMapAt[f_, ope_OPEData, position_List] :=
    MapAt[MapAt[f, #, position] & , ope, {1}]
(*
:[font = subsubsection; inactive; startGroup; Cclosed; backColorRed = 65535; backColorGreen = 65535; backColorBlue = 65535; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; bold; fontName = "Times New Roman"; fontSize = 10; ]
GetCoefficients, GetOperators
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; backColorRed = 65535; backColorGreen = 65535; backColorBlue = 65535; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; plain; fontName = "Courier New"; fontSize = 11; ]
*)
ClearAll[GetCoefficients]
SetAttributes[GetCoefficients,Listable]
GetCoefficients[something_] := Block[{ii,term},
        ii=1;
        OPESimplify[something,If[#===0,0,term[ii++] = #]&];
        Array[term,{ii-1}]
]

ClearAll[GetOperators]
SetAttributes[GetOperators, Listable]

GetOperators[A_OPEData] := Join @@ Map[GetOperators, First[A]]

GetOperators[a_Plus] := Union[Join@@ (GetOperators /@ (List @@ a))]
If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* old version of Mathematica *)
    GetOperators[a_Times] := 
           GetOperators@@Select[List@@a, OperatorQ]
    ,
    (* new version : use Select[ , ,1] to Get the first operator
       you find (there is only one !)
    *)
    GetOperators[a_Times] := 
           GetOperators@@Select[List@@a, OperatorQ, 1]
]
GetOperators[0] = {};
GetOperators[a_] := {a}
(*
:[font = subsubsection; inactive; startGroup; Cclosed; ]
try using LinOperate (does not work yet, too slow)
:[font = input; inactive; preserveAspect; endGroup; endGroup; nowordwrap; backColorRed = 65535; backColorGreen = 65535; backColorBlue = 65535; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; plain; fontName = "Courier New"; fontSize = 11; ]
PoleSimplify[A_] :=
        PoleSimplify[A, Sequence@@ Options[OPESimplify]]
PoleSimplify[term_, Function -> func_] :=
( (*tmp = *)
    LinCollect[term, GetOperators[term], func,
      Pattern -> _NonCommutativeMultiply](*;
  If[ Expand[tmp-term/. Global`cc->8] =!= 0,
      Print[tmp//InputForm]];
  tmp*)
)
PoleSimplify[A_,func_] :=
        PoleSimplify[A, Function->func]

GetCoefficients[A_OPEData] :=
        Flatten[GetCoefficients[First[A]]]
GetCoefficients[pole_] :=
   LinOperate[pole, ExtractOperators[pole], #2&, List,
      Pattern -> _NonCommutativeMultiply]
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Rules for Derivative
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[Derivative]
Derivative[i__][A_Plus] := Derivative[i] /@ A
Derivative[i__][A_ s_]  := s Derivative[i][A] /;  OperatorQ[A]
Derivative[_][0] = 0
 (* A comment on efficiency of the rule for a product given above.
    In Mathematica order, s follows A. Because Times is Orderless,
    this means that during pattern matching the A_ pattern will be
    matched to the factors in the product. So, for each factor
    OperatorQ is called once, until an operator is found.
        (a b OP)'     ?b OP (a')    OperatorQ[a]
                      ?a OP (b')    OperatorQ[b]
                      ?a b (OP')    OperatorQ[OP]
                  --> a b OP'
    With the rule
    Derivative[i__][B_ a_]  := a Derivative[i][B] /; OperatorQ[B]
    (a is ordered before B) you get
        (a b OP)'     ?a (b OP)'    OperatorQ[b OP]  (needs testing b and OP)
                  --> a (b OP)'
                      ? a b OP'     OperatorQ[OP]
                  --> a b OP'
    When more factors are present, this is obviously much less efficient
   ( n evaluations of the rules, (n-1)! or so OperatorQ )
 *)
CommutativeQ[Derivative[_][a_]]:= CommutativeQ[a]
ChiralQ[Derivative[_][a_]]:= ChiralQ[a]
AChiralQ[Derivative[_][a_]]:= AChiralQ[a]

Literal[Derivative[i_][NonCommutativeMultiply[a__,B_]]] :=
        NonCommutativeMultiply[a,Derivative[i][B]]
Derivative[n_][t12|tb12|t|tb] := 0 /; n>0
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Operator Handling
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
Clear[Bosonic,Fermionic,OPEOperator,OPEParity,
      OPEPosition,OPEOrder,OperatorQ, BosonQ,SwapSign]
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
declaration of SF
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
BosonicHelp[A_] := (If[OperatorQ[A], Message[Bosonic::operator, A]];
              OPEParity[A]=0; BosonQ[A]=True; OperatorQ[A]=True;
              ChiralQ[A]=False; AChiralQ[A]=False;
              OPEPosition[A]=OPEPositionCounter++;AppendTo[OperatorList, A])
FermionicHelp[A_] := (If[OperatorQ[A], Message[Fermionic::operator, A]];
              OPEParity[A]=1; BosonQ[A]=False; OperatorQ[A]=True;
              ChiralQ[A]=False; AChiralQ[A]=False;
              OPEPosition[A]=OPEPositionCounter++;AppendTo[OperatorList, A])
ChiralBosonicHelp[A_] := (BosonicHelp[A]; ChiralQ[A]=True)
ChiralFermionicHelp[A_] := (FermionicHelp[A]; ChiralQ[A]=True)
AChiralBosonicHelp[A_] := (BosonicHelp[A]; AChiralQ[A]=True)
AChiralFermionicHelp[A_] := (FermionicHelp[A]; AChiralQ[A]=True)

OPEOperatorHelp[A_,n_Integer] :=
     If[Mod[n,2]==0,BosonicHelp[A],FermionicHelp[A]]
OPEOperatorHelp[A_,p_] :=
     (If[OperatorQ[A], Message[OPEOperator::operator, A]];
      OPEParity[A]=p; OperatorQ[A]=True;
      ChiralQ[A]=False; AChiralQ[A]=False;
      OPEposition[A] = OPEpositionCounter++;
      AppendTo[OperatorList, A];
      SetOPEOptions[ParityMethod,1])

Bosonic::operator =
Fermionic::operator = 
OPEOperator::operator = "Warning : `1` is already declared operator.
Defining other OPEs as before for `1` may give wrong results."

OperatorList = {}
OPEPositionCounter = 0

GeneralBosonic=Bosonic;
GeneralFermionic=Fermionic;

Bosonic[A___] := Scan[BosonicHelp,{A}]
Fermionic[A___] := Scan[FermionicHelp,{A}]
ChiralBosonic[A___] := Scan[ChiralBosonicHelp,{A}]
ChiralFermionic[A___] := Scan[ChiralFermionicHelp,{A}]
AChiralBosonic[A___] := Scan[AChiralBosonicHelp,{A}]
AChiralFermionic[A___] := Scan[AChiralFermionicHelp,{A}]

OPEOperator[A:({_,_}..)] := Scan[OPEOperator, {A}]
OPEOperator[A_,p_] := (OPEOperatorHelp[A,p];Null)
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Order
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

 (******************** ordering of operators ***********************)
 (* OPEOrder[A,B] returns > 0 when A and B are ordered,
  *                       = 0 when A == B,
  *                       < 0 otherwise.
  * Ordering is fixed by the order in which the operators are
  * declared. If A and B match the same operator-pattern, use
  * standard Mathematica order.
  *)
(*
order!!
OPEOrder[a_, b_] :=
    Block[{res = OPEPosition[b] - OPEPosition[a]},
        If[res == 0, Order[a, b], res]
    ]
*)
OPEOrder[a_, b_] := OPEPosition[b] - OPEPosition[a]
 (*
 OPEOrder[A_,B_NO] := (Print["NOR ",A," ",B,Stack[_RuleCondition]];1)
 OPEOrder[A_NO,B_] := (Print["NOL ",A," ",B,Stack[_RuleCondition]];-1)
 *)
OPEOrder[A_,B_NO] =1
OPEOrder[A_NO,B_] =-1

(*OPEOrder[a_,a_]=0*)
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
OperatorQ
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

 (* OperatorQ[A] tests if A is an operator *)

OperatorQ[A_+B_] := OperatorQ[A]
 (*
OperatorQ[A_Plus] := (
        Or@@ (OperatorQ /@ (List@@A))
) *)
OperatorQ[A_Times] :=
        Or@@ (OperatorQ /@ (List@@A))

OperatorQ[A_NonCommutativeMultiply]=True

OperatorQ[Derivative[_][A_]] := OperatorQ[A]

 (* all other things are no operators *)
OperatorQ[_] = False
OperatorQ[0] = True
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
BosonQ, OPEParity
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
BosonQ[Derivative[_][A_]] := BosonQ[A]
OPEParity[Derivative[_][A_]] := OPEParity[A]
OPEParity[A_+B_]:= OPEParity[A]
OPEParity[a_ A_]:= OPEParity[A]/;!OperatorQ[a]
Literal[OPEParity[NonCommutativeMultiply[A__]]] :=
    Plus@@ (OPEParity /@ {A})
OPEParity[1] = 0;
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
theta handling
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
theta properties
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
(* TODO : weg ? *)
AChiralFermionic[t,t12]
ChiralFermionic[tb,tb12]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
DT1 DBT1 properties
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[DT1,DBT1]
 (* We will use that in **, thetas are always ordered first. Moreover,
    DT[theta] is always a number. This provides a slight simplification
    of Leibniz rule.
 *)
DT1[A_Plus] := DT1 /@ A
DT1[A_ s_]  := s DT1[A] /; OperatorQ[A]
DT1[0] = 0
Literal[DT1[NonCommutativeMultiply[a_,b__]]]:=
        DT1[a] NonCommutativeMultiply[b]+
        SwapSign[a] a**DT1[NonCommutativeMultiply[b]]
DT1[t12]=1
DT1[_]=0

DBT1[A_Plus] := DBT1 /@ A
DBT1[A_ s_]  := s DBT1[A] /; OperatorQ[A]
DBT1[0] = 0
Literal[DBT1[NonCommutativeMultiply[a_,b__]]]:=
        DBT1[a] NonCommutativeMultiply[b]+
        SwapSign[a] a**DBT1[NonCommutativeMultiply[b]]
DBT1[tb12]=1
DBT1[_]=0
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
DT DBT DDB properties
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
Clear[DT,DBT,DDB]
 (* We will use that in **, thetas are always ordered first. Moreover,
    DT[theta] is always a number. This provides a slight simplification
    of the Leibniz rule.
 *)

BosonQ[DT[a_]] := !BosonQ[a]
OPEParity[DT[a_]] := OPEParity[a]+1
OperatorQ[DT[a_]] := OperatorQ[a]
OPEPosition[DT[a_]] :=OPEPosition[a]
BosonQ[DBT[a_]] := !BosonQ[a]
OPEParity[DBT[a_]] := OPEParity[a]+1
OperatorQ[DBT[a_]] := OperatorQ[a]
OPEPosition[DBT[a_]] :=OPEPosition[a]
BosonQ[DDB[a_]] := BosonQ[a]
OPEParity[DDB[a_]] := OPEParity[a]+1
OperatorQ[DDB[a_]] := OperatorQ[a]
OPEPosition[DDB[a_]] :=OPEPosition[a]

DT[A_Plus] := DT /@ A
DT[A_ s_]  := s DT[A] /; OperatorQ[A]
DT[0] = 0
Literal[DT[NonCommutativeMultiply[a_,b__]]]:=
        DT[a] NonCommutativeMultiply[b]+
        SwapSign[a] a**DT[NonCommutativeMultiply[b]]
(*Literal[DT[DT[a_]]]=0*)
ChiralQ[DT[_]]=True;
DT[a_]:=0 /; (*!OperatorQ[a] ||*) ChiralQ[a]
DT[t]=1
DT[t12]=-1

DBT[A_Plus] := DBT /@ A
DBT[A_ s_]  := s DBT[A] /; OperatorQ[A]
DBT[0] = 0
Literal[DBT[NonCommutativeMultiply[a_,b__]]]:=
        DBT[a] NonCommutativeMultiply[b]+
        SwapSign[a] a**DBT[NonCommutativeMultiply[b]]
(*Literal[DBT[DBT[a_]]]=0*)
AChiralQ[DBT[_]]=True;
DBT[a_]:=0 /; (*!OperatorQ[a] || *) AChiralQ[a]
DBT[tb]=1
DBT[tb12]=-1

DDB[A_Plus] := DDB /@ A
DDB[A_ s_]  := s DDB[A] /; OperatorQ[A]
DDB[0] = 0
Literal[DDB[NonCommutativeMultiply[a_,b__]]]:=
        (*DDB[a] NonCommutativeMultiply[b]+*)
        SwapSign[a] (
            -2 DBT[a] DT[NonCommutativeMultiply[b]] +
            2 DT[a] DBT[NonCommutativeMultiply[b]]) +
        a**DDB[NonCommutativeMultiply[b]]
(*DDB[a_]:=0 /; !OperatorQ[a]*)
DDB[tb]=0
(*DDB[t12|tb12]=0*)

Literal[DT[DBT[a_]]]:=-a' /;  ChiralQ[a] ;
Literal[DBT[DT[a_]]]:=-a' /; AChiralQ[a] ;

(* ddb *)
(*Literal[ DDB[DT[a_]]]:=-DT[a'] ;
Literal[ DDB[DBT[a_]]]:=DBT[a'] ;*)
Literal[ DT[DDB[a_]]]:=DT[a'] ;
Literal[ DBT[DDB[a_]]]:=-DBT[a'];
Literal[ DDB[DDB[a_]]]:=a'' ;

Literal[DDB[a_]]:=-a' /;  ChiralQ[a] ;
Literal[DDB[a_]]:=a' /;  AChiralQ[a] ;

Literal[ DT[DBT[a_]]]:=1/2 DDB[a] - 1/2 a';
Literal[ DBT[DT[a_]]]:=-1/2 DDB[a] - 1/2 a';
(* end ddb *)
(* TODO omgekeerd ? *)
Literal[Derivative[n_][DT[a_]]] := DT[Derivative[n][a]]
Literal[Derivative[n_][DBT[a_]]] := DBT[Derivative[n][a]]
Literal[Derivative[n_][DDB[a_]]] := DDB[Derivative[n][a]]

(* DDB[x_] := DT[DBT[x]]-DBT[DT[x]] *)
(* 2 following rules not necessary when DDB is used *)
(* Literal[ DT[DBT[DT[a_]]]]:=-DT[a'] ;
   Literal[ DBT[DT[DBT[a_]]]]:=-DBT[a'] ;*)
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPE handling
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
Clear values
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[OPE,
      OPEDerivativeHelpR, OPEDerivativeHelpL,
      OPECompositeHelpR, OPECompositeHelpL,
      OPECommuteHelp]
OPE::trace = "`1`";
OPEPole::trace = "`1`";
OPEDerivativeHelpL::trace = "`1`";
OPEDerivativeHelpR::trace = "`1`";
OPECompositeHelpL::trace = "`1`";
OPECompositeHelpR::trace = "`1`";
OPECommuteHelp::trace = "`1`";
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
MaxPole
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (* MaxPole[ope] gives the order of the highest pole *)

MaxPole::other = "`1` should be a OPEData object, something wrong
here, assuming there are no poles in its OPE."

MaxPole[OPEData[A_List]] := Length[A]
MaxPole[A_] := (Message[MaxPole::other, A]; 0)

(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPE rules
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (********************** OPE properties ****************************)
 (* 11-06-93 (3.0 beta 10) added following rule *)
Literal[OPE[___,0,___]] = OPEData[{}];

 (*Linear*)
If[NumberQ[$VersionNumber],
    (*
    Literal[OPE[a___,b_Plus,c___]]:=
        Plus @@
            Distribute[
                Lineartmp[a,b,c],
                Plus,Lineartmp,
                Lineartmp2,OPE
        ]
    *)
    (* change 5-10-92 *)
    Literal[OPE[a___,b_Plus,c___]]:=
            Distribute[
                Lineartmp[a,b,c],
                Plus,Lineartmp,
                Plus,OPE
        ],
    (* bug in version 1.1 Distribute *)
    Literal[OPE[a___,b_+c_,d___]] :=
        OPE[a,b,d] + OPE[a,c,d]
]

Literal[OPE[A___,s_ B_,C___]] :=
         s OPE[A,B,C] /; OperatorQ[B]

Literal[OPE[NonCommutativeMultiply[theta__,A_],B_]]:=
        OPEMap[NonCommutativeMultiply[theta,#]&,OPE[A,B]]

Literal[OPE[A_,NonCommutativeMultiply[theta__,B_]]]:=
        Apply[Times,SwapSign[A,#]& /@{theta}]*
        OPEMap[NonCommutativeMultiply[theta,#]&,OPE[A,B]]

Literal[OPE[Derivative[i_][A_],B_]]:=
        OPEDerivativeHelpL[A,B,i]

Literal[OPE[DT[A_],B_]]:=
        OPEDTHelpL[OPESimplify[OPE[A,B]]]

Literal[OPE[DBT[A_],B_]]:=
        OPEDBTHelpL[OPESimplify[OPE[A,B]]]

Literal[OPE[DDB[A_],B_]]:=
        OPEDDBHelpL[OPESimplify[OPE[A,B]]]

Literal[OPE[A_,Derivative[i_][B_]]]:=
        OPEDerivativeHelpR[A,B,i]

Literal[OPE[A_,DT[B_]]]:=
        SwapSign[A] OPEDTHelpR[OPESimplify[OPE[A,B]]]

Literal[OPE[A_,DBT[B_]]]:=
        SwapSign[A] OPEDBTHelpR[OPESimplify[OPE[A,B]]]

Literal[OPE[A_,DDB[B_]]]:=
        OPEDDBHelpR[OPESimplify[OPE[A,B]]]

Literal[OPE[A_,NO[B_,C_]]] :=
        CallAndSave[OPECompositeHelpR,A,B,C]

Literal[OPE[NO[A_,B_],C_]] :=
        CallAndSave[OPECompositeHelpL,A,B,C] /;
                 Not[SameQ[Head[B], NO]]

Literal[OPE[B_,A_]] :=
        OPECommuteHelp[B,A] /;
                 SameQ[Head[B],NO] || OPEOrder[A,B]>0

 (* All non-defined OPEs are assumed to be zero *)
Literal[OPE[_,_]]= OPEData[{}];
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
common auxiliary functions
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
derivatives
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (* OPE[Derivative[i_][A_],B_]] *)
 (* 08-11-92 (3.0 beta 6)*)
OPEDerivativeHelpL[A_,B_,i_] :=
        OPEData[
            Block[{j, AB = OPE[A,B]},
                Join[(-1)^i *
                        Table[Pochhammer[j,i]  opepole[j][AB],
                            {j,MaxPole[AB],1,-1}
                        ],
                     Table[0, {i}]
                ]
            ]
        ]

OPEDTHelpL[AB_] :=
        OPEData[
            Join[
                Range[MaxPole[AB],1,-1]/2 *(tb12**#&/@AB[[1]]),
                {0}
            ] +
            Join[{0}, DT1 /@ AB[[1]]]
        ]

OPEDBTHelpL[AB_] :=
        OPEData[
            Join[
                Range[MaxPole[AB],1,-1]/2 *(t12**#&/@AB[[1]]),
                {0}
            ] +
            Join[{0}, DBT1 /@ AB[[1]]]
        ]

OPEDDBHelpL[AB_] :=
        OPEDTHelpL[OPEDBTHelpL[AB]] -
        OPEDBTHelpL[OPEDTHelpL[AB]]

 (* OPE[A_,Derivative[i_][B_]] *)
 (* 08-11-92 (3.0 beta 6)*)
OPEDerivativeHelpR[A_,B_,i_] :=
    Block[{der, j, k, AB = OPE[A,B], maxAB},
        maxAB = MaxPole[AB];
        der[0] = Reverse[ AB[[1]] ];
        Do[der[j] = Map[Derivative[1], der[j-1]], {j,i}];
        OPEData[
            Table[
                Sum[der[i-k][[j-k]] Binomial[i,k] Pochhammer[j-k,k],
                    {k, Max[0, j-maxAB], Min[i, j-1]}],
                {j,maxAB+i,1,-1}]
        ]
    ]

OPEDTHelpR[AB_] :=
        OPEData[
            Join[
                Range[MaxPole[AB],1,-1]/2 *(tb12**#&/@AB[[1]]),
                {0}
            ] +
            Join[{0}, DT /@ AB[[1]]]
        ]

OPEDBTHelpR[AB_] :=
        OPEData[
            Join[
                Range[MaxPole[AB],1,-1]/2 *(t12**#&/@AB[[1]]),
                {0}
            ] +
            Join[{0}, DBT /@ AB[[1]]]
        ]


OPEDDBHelpR[AB_] :=
        OPEDTHelpR[OPEDBTHelpR[AB]] -
        OPEDBTHelpR[OPEDTHelpR[AB]]
(*
:[font = text; inactive; preserveAspect; startGroup; Cclosed; ]
commutation
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
 (* OPE[B,A] in terms of OPE[A,B]
  * Formula:
  * [BA](q) = Sum[(-1)^l / (l-q)! D[[AB](l),{.,l-q}],
  *               {l,q,MaxPole[AB]}]
  *)
Taylor[pole_] :=
    PoleSimplify[
        PoleSimplify[pole] /.
            {t12->-t12,tb12->-tb12} /.
            x:Alternatives@@ExtractOperatorsNoTheta[pole] :>
                   x + t12**DT[x] + tb12**DBT[x] -
                   1/2 t12**tb12**DDB[x]
    ]

OPECommuteHelp[B_,A_] :=
    Block[{q,l,term,AB = OPE[A,B], max},
        max = MaxPole[AB];
        OPEData[
            SwapSign[A,B] *
            Table[
                (term[q] = (-1)^q Taylor[opepole[q][AB]]
                ) +
                Sum[term[l] = (term[l]') / (l-q),
                      {l,q+1,max}
                ],
                {q,max,1,-1}
            ]
        ]
    ]
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
auxiliary functions quantum case
:[font = text; inactive; preserveAspect; startGroup; Cclosed; ]
composites at the right
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
Clear[proj]
thetas12[0]=1;
thetas12[1]=t12;
thetas12[2]=tb12;
thetas12[3]=t12**tb12;
Dops[0]=Identity;
Dops[1]=DT;
Dops[2]=DTB;
Dops[3]=DDB;
proj[0][x_] := x/.t12|tb12->0
proj[1][x_Plus] := Select[x,!FreeQ[#,t12]&&FreeQ[#,tb12]&]/.t12->1
proj[2][x_Plus] := Select[x,FreeQ[#,t12]&&!FreeQ[#,tb12]&]/.tb12->1
proj[3][x_Plus] := Select[x,!FreeQ[#,t12]&&!FreeQ[#,tb12]&]/.t12|tb12->1
proj[1][x_] := If[!FreeQ[x,t12]&&FreeQ[x,tb12],x/.t12->1,0]
proj[2][x_] := If[FreeQ[x,t12]&&!FreeQ[x,tb12],x/.tb12->1,0]
proj[3][x_] := If[!FreeQ[x,t12]&&!FreeQ[x,tb12],x/.t12|tb12->1,0]
(*
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

 (* OPE[A,NO[B,C]] *)
ClearOPECompositeHelpRQ[] := (OPECompositeHelpRQ[A_,B_,C_] =.)

RedefineOPECompositeHelpRQ[] := (
    OPECompositeHelpRQ[A_,B_,C_] :=
        Block[{q,l,sign = SwapSign[A,B], ABC, AB, AC,
                 maxAB, maxABC, maxq},
            AB = OPESimplify[OPE[A,B]];
            AC = If[ SameQ[B,C], AB, OPE[A,C]];
            maxAB = MaxPole[AB];
            ABC = Table[
                      Sum[OPEMap[thetas12[l]**proj[0][#]&,
                                OPE[proj[l][opepole[q][AB]], C]
                          ],{l,0,3}
                      ],
                      {q,maxAB}
                  ];
            (* change 6-4-92
               maxABC = Max[MaxPole /@ ABC];
               maxq = Max[maxABC + maxAB, MaxPole[AC]];
            *)
            maxABC = MaxPole /@ ABC;
            maxq = Max[maxABC + Range[maxAB], MaxPole[AC]];
            maxABC = Max[maxABC,0];
            OPEData[
                Table[
                    PoleSimplify[
                        sign * NO[B,OPEPole[q][AC]] +
                        NO[OPEPole[q][AB],C] +
                        Sum[Binomial[q-1,l] *
                            OPEPole[l][ ABC[[q-l]]] ,
                            {l,Max[1,q-maxAB], Min[q-1, maxABC]}
                        ]
                    ],
                    {q,maxq,1,-1}
                ]
            ]
       ]
)
RedefineOPECompositeHelpRQ[]
(*
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites at the left
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
(* OPE[NO[A,B],C] *)
ClearOPECompositeHelpLQ[] := (OPECompositeHelpLQ[A_,B_,C_] =.)

RedefineOPECompositeHelpLQ[] := (
    OPECompositeHelpLQ[A_,B_,C_] :=
          Block[{AC, BC, q,l,sign = SwapSign[A,B], BAC,
                 maxAC, maxBC, maxBAC, maxq, der, res},
            AC = OPESimplify[OPE[A,C]];
            maxAC = MaxPole[AC];
            BAC = Table[
                      OPE[B, opepole[q][AC]],
                      {q,maxAC}
                  ];
            (* der[l][A] = Derivative[l][A//Taylor] *)
            der[0][B] = B + t12**DT[B] + tb12**DBT[B] -
                   1/2 t12**tb12**DDB[B];
            der[l_][X_] := der[l][X] =
                PoleSimplify[der[l-1][X]', Together];
            res = OPEData[
                    Table[
                        Sum[ NO[der[l][B], opepole[l+q][AC]] / l!,
                            {l,0, maxAC-q}
                        ],
                        {q,maxAC,1,-1}
                    ]
                  ];
            If[ SameQ[A,B],
                res = (1+sign)*res;
                , (* else *)
                BC= OPE[B,C];
                maxBC = MaxPole[BC];
                der[0][A] = A + t12**DT[A] + tb12**DBT[A] -
                   1/2 t12**tb12**DDB[A];
                res = sign*res +
                        OPEData[
                            Table[
                                Sum[ NO[der[l][A], opepole[l+q][BC]] / l!,
                                    {l,0, maxBC-q}
                                ],
                                {q,maxBC,1,-1}
                            ]
                        ];
            ];
            (* change 6-4-92
               maxBAC = Max[MaxPole /@ BAC];
               maxq = Max[maxBAC + maxAC, 0];
            *)
            maxBAC = MaxPole /@ BAC;
            (* Special case : maxAC==0,
             * then Max[maxBAC]==-Infinity. However,
             * Sum[i,{i,1,-Infinity}] does not evaluate to 0 as expected.
             * So we need the Max[...,0].
             *)
            maxq = Max[maxBAC + Range[maxAC], 0];
            maxBAC = Max[maxBAC, 0];
            OPESimplify[
                res  +
                OPEData[sign*
                    Table[
                        Sum[ OPEPole[l][ BAC[[q-l]]] ,
                            {l,Max[1,q-maxAC], Min[q-1, maxBAC]}
                        ],
                        {q,maxq,1,-1}
                    ]
                ], Together
            ]
       ]
)
RedefineOPECompositeHelpLQ[]
(*
:[font = subsubsection; inactive; preserveAspect; startGroup; Cclosed; ]
auxiliary functions Poisson case
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites at the right
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

 (* OPE[A,NO[B,C]] *)
ClearOPECompositeHelpRPB[] := (OPECompositeHelpRPB[A_,B_,C_] =.)

RedefineOPECompositeHelpRPB[] := (
    OPECompositeHelpRPB[A_,B_,C_] :=
        Block[{q,l,sign = SwapSign[A,B], AB, AC,
                 maxAB, maxq},
            AB = OPE[A,B];
            AC = If[ SameQ[B,C], AB, OPE[A,C]];
            maxAB = MaxPole[AB];
            maxq = Max[MaxPole[AC],MaxPole[AB]];
            OPEData[
                Table[
                    PoleSimplify[
                        sign * NO[B,OPEPole[q][AC]] +
                        NO[OPEPole[q][AB],C]
                    ],
                    {q,maxq,1,-1}
                ]
            ]
       ]
)
RedefineOPECompositeHelpRPB[]
(*
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites at the left
:[font = input; initialization; preserveAspect; endGroup; endGroup; endGroup; nowordwrap; ]
*)
(* OPE[NO[A,B],C] *)
ClearOPECompositeHelpLPB[] := (OPECompositeHelpLPB[A_,B_,C_] =.)

RedefineOPECompositeHelpLPB[] := (
    OPECompositeHelpLPB[A_,B_,C_] :=
          Block[{AC, BC, q,l,sign = SwapSign[A,B],
                 maxAC, maxBC, maxq, der},
            AC = OPESimplify[OPE[A,C]];
            maxAC = MaxPole[AC];
            (* der[l][A] = Derivative[l][A//Taylor] *)
            der[0][B] = B + t12**DT[B] + tb12**DBT[B] -
                   1/2 t12**tb12**DDB[B];
            der[l_][X_] := der[l][X] =
                PoleSimplify[der[l-1][X]', Together];
            res = OPEData[
                    Table[
                        Sum[ NO[der[l][B], opepole[l+q][AC]] / l!,
                            {l,0, maxAC-q}
                        ],
                        {q,maxAC,1,-1}
                    ]
                  ];
            If[ SameQ[A,B],
                res = (1+sign)*res;
                , (* else *)
                BC= OPE[B,C];
                maxBC = MaxPole[BC];
                der[0][A] = A + t12**DT[A] + tb12**DBT[A] -
                   1/2 t12**tb12**DDB[A];
                res = sign*res +
                        OPEData[
                            Table[
                                Sum[ NO[der[l][A], opepole[l+q][BC]] / l!,
                                    {l,0, maxBC-q}
                                ],
                                {q,maxBC,1,-1}
                            ]
                        ];
            ];
            OPESimplify[res, Together]
       ]
)
RedefineOPECompositeHelpLPB[]
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPEPole
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
Clear values
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[OPEPole,OPEPoleHelpR,OPEPoleHelpL]
OPEPoleHelpL::trace = "`1`";
OPEPoleHelpR::trace = "`1`";
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
OPEPole rules
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (* 11-06-93 (3.0 beta 10) added following rule *)
Literal[OPEPole[_][___,0,___]] = 0;

 (* change 25-9-92 *)
 (*
OPEPole[n_Integer][A_OPEData] :=
        If [1 <= n <= Length[A[[1]]],
             A[[1, Length[A[[1]]]+1-n]],
             0
        ]
*)
OPEPole[n_Integer][OPEData[A_List]] :=
        If [1 <= n <= Length[A],
             A[[-n]],
             0
        ]

 (* opepole introduced 6-6-93 (3.0 beta 11) *)
opepole[n_Integer][OPEData[A_List]] :=
        If [1 <= n <= Length[A],
             A[[-n]],
             Message[opepole::range,n,OPEData[A]];0
        ]
opepole::range = "WARNING : pole `1` out of range in `2`";

If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0,
    (* bug in old versions : OPEData[{a}] - OPEData[{a}] -> 0 *)
    OPEPole[_][0] = 0;
    opepole[_][0] = 0
]

Literal[OPEPole[0][A_,B_]] :=
        NO[A,B] + t12**NO[DT[A],B] + tb12**NO[DBT[A],B] -
        1/2 t12**tb12**NO[DDB[A],B]

Literal[OPEPole[i_][A_,B_]] :=
        1/(-i)! OPEPole[0][Derivative[-i][A],B] /; i<0

 (* Linear *)
If[!NumberQ[$VersionNumber],
    (* other implementation : Distribute doesn't want OPEPole[i] as
       fifth argument.
    *)
    Literal[OPEPole[i_][a___,b_Plus,c___]]:=
        Distribute[
            Lineartmp[a,b,c]
        ] /. Lineartmp -> OPEPole[i],
    Literal[OPEPole[i_][a___,b_+c_,d___]] :=
        OPEPole[i][a,b,d] + OPEPole[i][a,c,d]
]
Literal[OPEPole[i_][A___,B_ s_,C___]] :=
        s OPEPole[i][A,B,C] /; OperatorQ[B]

Literal[OPEPole[_][NonCommutativeMultiply[___,th_,___,_],
        NonCommutativeMultiply[___,th_,___,_]]] = 0;

Literal[OPEPole[i_][NonCommutativeMultiply[theta__,A_],B_]]:=
        NonCommutativeMultiply[theta,OPEPole[i][A,B]]

Literal[OPEPole[i_][A_,NonCommutativeMultiply[theta__,B_]]]:=
        Apply[Times,SwapSign[A,#]& /@{theta}]*
        NonCommutativeMultiply[theta,OPEPole[i][A,B]]

Literal[OPEPole[i_][Derivative[j_][A_],B_]] :=
        If[ i>j,
                Pochhammer[1-i,j] * OPEPole[i-j][A,B],
                0
        ]

Literal[OPEPole[i_][DT[A_],B_]] :=
        If[i==1, 0,
          (i-1)/2  tb12**OPEPole[i-1][A,B]
        ] +
        DT1[OPEPole[i][A,B]]

Literal[OPEPole[i_][DBT[A_],B_]] :=
        If[i==1, 0,
          (i-1)/2  t12**OPEPole[i-1][A,B]
        ] +
        DBT1[OPEPole[i][A,B]]

Literal[OPEPole[i_][DDB[A_],B_]] :=
        If[i==1, 0,
          (i-1)/2 (tb12**OPEPole[i-1][DBT[A],B] -
                   t12**OPEPole[i-1][DT[A],B])
        ] +
        DT1[OPEPole[i][DBT[A],B]] -
        DBT1[OPEPole[i][DT[A],B]]

Literal[OPEPole[j_][A_, Derivative[i_][B_]]] :=
        Block[{k},
                Sum[Pochhammer[j-k,k] Binomial[i,k]*
                                Derivative[i-k][OPEPole[j-k][A,B]],
                        {k, 0, Min[i,j-1]}
                ]
        ] /; i>0
Literal[OPEPole[i_][DT[A_],B_]] :=
        If[i==1, 0,
          (i-1)/2  tb12**OPEPole[i-1][A,B]
        ] +
        DT1[OPEPole[i][A,B]]

Literal[OPEPole[i_][A_,DT[B_]]] :=
        SwapSign[A] *
        (   If[i==1, 0,
              (i-1)/2  tb12**OPEPole[i-1][A,B]
            ] +
            DT[OPEPole[i][A,B]]
        )

Literal[OPEPole[i_][A_,DBT[B_]]] :=
        SwapSign[A] *
        (   If[i==1, 0,
              (i-1)/2  t12**OPEPole[i-1][A,B]
            ] +
            DBT[OPEPole[i][A,B]]
        )

Literal[OPEPole[i_][A_,DDB[B_]]] :=
        SwapSign[A] *
        (   If[i==1, 0,
              (i-1)/2 (tb12**OPEPole[i-1][A,DBT[B]] -
                   t12**OPEPole[i-1][A,DT[B]])
            ] +
            DT[OPEPole[i][A,DBT[B]]] -
            DBT[OPEPole[i][A,DT[B]]]
        )

Literal[OPEPole[i_][A_,NO[B_,C_]]] :=
        OPEPoleHelpR[A,B,C,i]

Literal[OPEPole[i_][NO[A_,B_],C_]] :=
        OPEPoleHelpL[A,B,C,i]

 (* 7-07-93 (3.0 beta 11) added following rule.*)
Literal[OPEPole[i_][B_,A_]] :=
        OPEPoleCommuteHelp[i][B,A] /; OPEOrder[A,B]>0

Literal[OPEPole[i_][A_,C_]] :=
        OPEPole[i][OPE[A,C]]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
common auxiliary functions
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPEPoleCommuteHelp[q_][B_,A_] :=
    Block[{l,AB = OPE[A,B], sign = SwapSign[A,B] },
        Sum[ sign (-1)^l/(l-q)!*
              Derivative[l-q][Taylor[opepole[l][AB]]],
           {l,q,MaxPole[AB]}
        ]
    ]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
auxiliary functions quantum case
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites at the right
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)

 (* OPEPole[q][A,NO[B,C]] *)
(*
OPEPoleHelpRQ[A_,B_,C_,q_]:=
    Block[{BA, l,sign = SwapSign[A,B], maxBA},
        BA = OPE[B,A];
        maxBA = MaxPole[BA];
        PoleSimplify[
            sign * NO[B, OPEPole[q][A,C]] + sign*
            Sum[(-1)^l *
                OPEPole[q-l][ opepole[l][BA], C ],
                {l,maxBA}
            ]
        ]
    ]/; SameQ[Head[A],NO] || OPEOrder[A,B]<0
*)
(* Assumes q>0 *)
OPEPoleHelpRQ[A_,B_,C_,q_]:=
    Block[{AB, l,sign = SwapSign[A,B], maxAB},
        AB = OPESimplify[OPE[A,B]];
        maxAB = MaxPole[AB];
        PoleSimplify[
            sign * NO[B,
                If[SameQ[B,C],OPEPole[q][AB], OPEPole[q][A,C]]
            ] +
            NO[OPEPole[q][AB],C] +
            Sum[Binomial[q-1,l] *
                Sum[thetas12[p]**proj[0][
                           OPEPole[l][proj[p][OPEPole[q-l][AB]], C]
                    ],{p,0,3}
                ],
                {l,Max[1,q-maxAB], q-1}
            ]
        ]
    ]
(*
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites at the left
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
(* OPEPole[q][NO[A,B],C] *)
OPEPoleHelpLQ[A_,B_,C_,q_] :=
      Block[{AC, BC, l,sign = SwapSign[A,B],
             maxAC, maxBC, der, res},
            AC = OPESimplify[OPE[A,C]];
            maxAC = MaxPole[AC];
            (* der[l][A] = Derivative[l][A//Taylor] *)
            der[0][B] = B + t12**DT[B] + tb12**DBT[B] -
                   1/2 t12**tb12**DDB[B];
            der[l_][X_] := der[l][X] =
                PoleSimplify[der[l-1][X]', Together];
            res = Sum[ NO[der[l][B], opepole[l+q][AC]] / l!,
                       {l,0, maxAC-q}
                  ];
            If[ SameQ[A,B],
                res = (1+sign)*res;
                , (* else *)
                BC = OPESimplify[OPE[B,C]];
                maxBC = MaxPole[BC];
                der[0][A] = A + t12**DT[A] + tb12**DBT[A] -
                   1/2 t12**tb12**DDB[A];
                res = sign*res +
                    Sum[
                        NO[der[l][A],opepole[l+q][BC]] /l!,
                        {l,0, maxBC-q}
                    ]
            ];
            PoleSimplify[
                res  +
                sign*
                    Sum[ OPEPole[l][B, opepole[q-l][AC] ] ,
                        {l,Max[1,q-maxAC], q-1}
                    ],
                Together
            ]
       ]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
auxiliary functions Poisson case
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
(*
OPEPoleHelpRPB[A_,B_,C_,q_]:=
    Block[{BA, l,sign = SwapSign[A,B], maxBA},
        BA = OPE[B,A];
        maxBA = MaxPole[BA];
        PoleSimplify[
            sign * NO[B, OPEPole[q][A,C]] + sign*
            Sum[(-1)^l *
                OPEPole[q-l][ opepole[l][BA], C ],
                {l,maxBA}
            ]
        ]
    ]/; SameQ[Head[A],NO] || OPEOrder[A,B]<0
*)
OPEPoleHelpRPB[A_,B_,C_,q_]:=
    PoleSimplify[
            SwapSign[A,B] * NO[B, OPEPole[q][A,C]] +
            NO[OPEPole[q][A,B],C]
    ]

OPEPoleHelpLPB[A_,B_,C_,q_] :=
      Block[{AC, BC, l,sign = SwapSign[A,B],
        maxAC, maxBC, tA, tB},
        AC = OPESimplify[OPE[A,C]];
        maxAC = MaxPole[AC];
        (* der[l][A] = Derivative[l][A//Taylor] *)
        der[0][B] = B + t12**DT[B] + tb12**DBT[B] -
                   1/2 t12**tb12**DDB[B];
        der[l_][X_] := der[l][X] =
                PoleSimplify[der[l-1][X]', Together];
        res = Sum[ NO[der[l][B], opepole[l+q][AC]] / l!,
                       {l,0, maxAC-q}
        ];
        PoleSimplify[
            If[ SameQ[A,B],
                res = (1+sign)*res;
                , (* else *)
                BC = OPESimplify[OPE[B,C]];
                maxBC = MaxPole[BC];
                der[0][A] = A + t12**DT[A] + tb12**DBT[A] -
                   1/2 t12**tb12**DDB[A];
                res = sign*res +
                    Sum[
                        NO[der[l][A],opepole[l+q][BC]] /l!,
                        {l,0, maxBC-q}
                    ]
            ]
        ]
    ]
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
NO handling
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
Clear values
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[NO,NOCommuteHelp,NOCompositeHelpR,NOCompositeHelpL,
      NODerivativeHelp, NOOrder,NOOrderHelp,NOOrderHelp2]
NOCompositeHelpL::trace = "`1`"
NOCompositeHelpR::trace = "`1`"
NOCommuteHelp::trace = "`1`"

(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
Relation with Operators
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OperatorQ[_NO]=True
Literal[BosonQ[NO[A_,B_]]]:= !Xor[BosonQ[A],BosonQ[B]]
Literal[OPEParity[NO[A_,B_]]]:= OPEParity[A] + OPEParity[B]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
NOOrder
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (* NOOrder is the same as OPEOrder, except that derivatives
  * (which can not occur as argument of OPEOrder) are ignored.
  * Only in the case of comparing derivatives of the same operator.
  * Then higher derivatives come after lower if NOOrderingValue > 0.
  *)
 (*!!
NOOrderHelp[Derivative[_][A_],B_] := NOOrderHelp[A,B]
NOOrderHelp[A_,Derivative[_][B_]] := NOOrderHelp[A,B]
 *)
NOOrderHelp[A_,B_] :=
        OPEOrder[A/.Derivative[_]->Identity,
                B/.Derivative[_]->Identity]
 (* include these when NOCompositeHelpL is not used *)
 (*
NOOrder[A_,B_NO] =1
NOOrder[A_NO,B_] =-1
 *)
NOOrder[A_,B_] :=
    Block[{res = NOOrderHelp[A,B]},
        If[res == 0,
           res=NOOrderHelp2[B]-NOOrderHelp2[A];
           If[res == 0,res = Order[A,B]]
        ];
        res
    ]
NOOrderHelp2[(DT|DBT|DDB)[Derivative[i_][_]]] :=
        NOOrderingValue * i
NOOrderHelp2[Derivative[i_][_]] := NOOrderingValue * i
NOOrderHelp2[_] = 0;
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
NOs of thetas
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Literal[NO[A_,t12]] ^:= A**t12;
Literal[NO[t12,A_]] ^:= t12**A;
Literal[NO[A_,tb12]] ^:= A**tb12;
Literal[NO[tb12,A_]] ^:= tb12**A;
Literal[NO[A_,t12**tb12]] := t12**tb12**A;
Literal[NO[t12**tb12,A_]] := t12**tb12**A;
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
NO rules
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Literal[NO[A_,B_,C__]] := NO[A,NO[B,C]]

Literal[NO[0,_]] = 0;
Literal[NO[_,0]] = 0;

 (* Linear *)
If[NumberQ[$VersionNumber],
    Literal[NO[a___,b_Plus,c___]]:=
        Distribute[
            Lineartmp[a,b,c],
            Plus,Lineartmp,
            Plus,NO
        ],
    Literal[NO[a___,b_+c_,d___]] :=
        NO[a,b,d] + NO[a,c,d]
]

Literal[NO[A___,B_ s_,C___]] := s NO[A,B,C] /; OperatorQ[B]

Literal[NO[NonCommutativeMultiply[___,th_,___,_],
        NonCommutativeMultiply[___,th_,___,_]]] = 0;
Literal[NO[NonCommutativeMultiply[theta__,A_],B_]] :=
        NonCommutativeMultiply[theta,NO[A,B]]
Literal[NO[A_,NonCommutativeMultiply[theta__,B_]]] :=
        Apply[Times,SwapSign[A,#]& /@{theta}]*
        NonCommutativeMultiply[theta,NO[A,B]]

(* 08-11-92 beta 7
  * Use If statements instead of relying on NOOrder[_,_NO]=1. This
  * is safer, more clear, and skips an additional evaluation of NOOrder.
  *)

Literal[NO[NO[A_,B_],NO[C_,D_]]] :=
        If[NOOrder[A,C]<=0,
           CallAndSave[NOCompositeHelpR,NO[A,B],C,D] (* -> {C,{{A,B},D}}*),
           CallAndSave[NOCompositeHelpL,A,B,NO[C,D]] (* -> {A,{B,{C,D}}}*)
        ]

Literal[NO[NO[A_,B_],C_]] :=
        If[NOOrder[A,C]>0,
           CallAndSave[NOCompositeHelpL,A,B,C]       (* -> {A,{B,C}}*),
           SwapSign[C,NO[A,B]] (NO[C,NO[A,B]] - NOCommuteHelp[C,NO[A,B]])
        ]

Literal[NO[B_,NO[A_,C_]]] :=
        CallAndSave[NOCompositeHelpR,B,A,C] /;
        NOOrder[A,B]>0

 (* Special case: A=B fermionic *)
Literal[NO[A_,NO[A_,C_]]] := 1/2 NO[NOCommuteHelp[A,A],C] /;
        !BosonQ[A]

 (* Fermionic case :
    NO[B,A] = -NO[A,B] + NOCommuteHelp[A,B]
    Replacing B with A, we get the following rule
 *)
Literal[NO[A_,A_]] :=
        NOCommuteHelp[A,A]/2 /; !BosonQ[A]

Literal[NO[B_,A_]] :=
        SwapSign[A,B] (NO[A,B] - NOCommuteHelp[A,B]) /;
        !SameQ[Head[A],NO] && NOOrder[A,B]>0
 (* Rewrite NO[A,B]' in terms of NO[A,B'] + ...*)
NO /: Literal[Derivative[i_Integer][NO[A_,B_]]] :=
        NODerivativeHelp[i,A,B] /; i>0
Literal[DT[NO[A_,B_]]] :=
        NO[DT[A],B] + SwapSign[A] NO[A,DT[B]]
Literal[DBT[NO[A_,B_]]] :=
        NO[DBT[A],B] + SwapSign[A] NO[A,DBT[B]]
Literal[DDB[NO[A_,B_]]]:=
        NO[DDB[A],B]+
        SwapSign[A] (
            -2 NO[DBT[A],DT[B]] +
            2 NO[DT[A],DBT[B]]) +
        NO[A,DDB[B]]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
common auxiliary functions
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
derivatives
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
NODerivativeHelp[i_,A_,B_] := Block[{j},
    PoleSimplify[
        Sum[Binomial[i,j] *
            NO[Derivative[j][A],Derivative[i-j][B]],
          {j,0,i}
        ]
    ]
]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
auxiliary functions quantum case
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
commutation
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
 (* definition:
  * NOCommuteHelp[A,B] = NO[A,B]-SwapSign[A,B] NO[B,A]
  * see for actual computation comment on OPECommuteHelp
  *)

NOCommuteHelpQ[A_,B_] :=
    Block[{m, AB = OPE[A,B]},
        Sum[-(-1)^m /m! Derivative[m][
                opepole[m][AB] /.{t12|tb12->0}
                ],
                {m,MaxPole[AB]}
        ]
    ]
(*
:[font = text; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
composites
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
 (* NO[B, NO[A,C]] = (-) NO[A,NO[B,C]] +
  *                  NO[NO[B,A] - (-)NO[A,B], C]
  *)
ClearNOCompositeHelpRQ[] := (NOCompositeHelpRQ[B_,A_,C_] =.)

RedefineNOCompositeHelpRQ[] := (
    NOCompositeHelpRQ[B_,A_,C_] :=
            PoleSimplify[
                SwapSign[A,B] NO[A, NO[B, C]] +
                NO[NOCommuteHelp[B, A], C]
            ]
)
RedefineNOCompositeHelpRQ[]

 (* NO[NO[A, B], C] *)
ClearNOCompositeHelpLQ[] := (NOCompositeHelpLQ[A_,B_,C_] =.)

(* change 1.0 beta 1
   Moved projections on theta-free part to the definitions
   of AC,BC. They used to be after the opepole[l][AC]
*)
RedefineNOCompositeHelpLQ[] := (
NOCompositeHelpLQ[A_,B_,C_] :=
    Block[{AC, BC, l,sign = SwapSign[A,B],
             maxAC, maxBC, res},
        AC = OPE[A,C]/.t12|tb12->0;
        maxAC = MaxPole[AC];
        res = Sum[ NO[Derivative[l][B], opepole[l][AC]] /l!,
                {l, maxAC}
            ];
        PoleSimplify[
            NO[A,NO[B,C]] +
            If[SameQ[A,B], (* A=B will always be bosonic because of NO[A,A] *)
                (1+sign)res 
                , (* else *)
                BC = OPE[B,C]/.t12|tb12->0;
                maxBC = MaxPole[BC];
                sign res +
                Sum[ NO[Derivative[l][A],
                         opepole[l][BC]] /l!,
                    {l, maxBC}
                ]
            ]
        ]
    ]
)
RedefineNOCompositeHelpLQ[]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
auxiliary functions Poisson case
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)
NOCommuteHelpPB[B_,A_] := 0

ClearNOCompositeHelpRPB[] := (NOCompositeHelpRPB[B_,A_,C_] =.)

(* NO[B, NO[A,C]] = (-) (NO[A,NO[B,C]]  *)
RedefineNOCompositeHelpRPB[] := (
    NOCompositeHelpRPB[B_,A_,C_] := SwapSign[A,B] NO[A,NO[B,C]]
)
RedefineNOCompositeHelpRPB[]

(* NO[NO[A, B], C] *)
ClearNOCompositeHelpLPB[] := (NOCompositeHelpLPB[A_,B_,C_] =.)

RedefineNOCompositeHelpLPB[] := (
NOCompositeHelpLPB[A_,B_,C_] := NO[A,NO[B,C]]
)
RedefineNOCompositeHelpLPB[]
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
The operator One
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
DefineConstantOperator[op_Symbol] := (
    OPEParity[op]^=0; BosonQ[op]^=True; OperatorQ[op]^=True;
    Derivative[n_][op] ^:= 0 /; n>0;
    DT[op] ^= 0;
    DBT[op] ^= 0;
    DDB[op] ^= 0;
    AppendTo[OperatorList, op];
    Literal[OPE[_,op]] ^:= OPEData[{}];
    Literal[OPE[op,_]] ^:= OPEData[{}];
    Literal[NO[A_,op]] ^:= A;
    Literal[NO[op,A_]] ^:= A;
)
Clear[One]
DefineConstantOperator[One]
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Components
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
N2OPEToComponents[A_, B_] :=
    N2OPEToComponents[OPESimplify[OPE[A,B]],A, B]
(* change version 1.0 : 
   added factor -1/2 for the last component
*)
(* provide default values : A and B bosonic *)
N2OPEToComponents[ope_OPEData,A_:One, B_:One] :=
    Block[{tmp = OPESimplify[
                {ope, OPEDTHelpL[ope],
                 OPEDBTHelpL[ope],-1/2 OPEDDBHelpL[ope]}],
           signs = If[BosonQ[A],{1,-1,-1,1},{-1,1,1,-1}]
                (*SwapSign[#]&/@{A,DT[A],DT[A],A}*)
          },
        {proj[0][tmp],
         signs proj[0][OPEDTHelpR/@tmp],
         signs proj[0][OPEDBTHelpR/@tmp],
         -1/2 proj[0][OPEDDBHelpR/@tmp]}//Transpose
    ]//OPESimplify
N2OPEToComponents[A_] :=
    {A,DT[A],DBT[A],-1/2 DDB[A]}//OPESimplify
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
OPEJacobi
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
derLop[0]=Identity;
derLop[1]=OPEDTHelpL;
derLop[2]=OPEDBTHelpL;
derLop[3]=-1/2 OPEDDBHelpL[#]&;

Clear[OPEJacobi]

Options[OPEJacobi] = {Function -> Expand};

OPEJacobi[A_,B_,C_,opts___Rule] :=
    Block[{AB, AC, BC,
          AnBC, BnAC, ABnC, m,n,p, sign = SwapSign[A,B],
          maxAB, maxAC, maxBC, maxAnBC, maxBnAC, maxABnC,
          maxn, maxm, x, y, dLAB,
          simopts= Function -> 
              (Function /. {opts} /. Options[OPEJacobi])},
        AB = OPESimplify[OPE[A,B], simopts];
        AC = OPESimplify[OPE[A,C], simopts];
        maxAB = MaxPole[AB];
        maxAC = MaxPole[AC];
        dLAB[y_]:=dLAB[y]=OPESimplify[derLop[y][AB], simopts];
        AnBC[y_,n_] := (AnBC[y,n]=
            OPESimplify[OPE[A, proj[y][OPEPole[n][BC]]], simopts])/;n<=maxBC;
        ABnC[y_,n_] :=
    		 (ABnC[y,n]=OPESimplify[OPE[OPEPole[n][
    				dLAB[y]/.{t12->-t12,tb12->-tb12}
    			], C], simopts])/;n<=(maxAB + Ceiling[y/Nsusy]);
        AnBC[__] = ABnC[__] = OPEData[{}];
        If[SameQ[A,B],
            BnAC=AnBC; BC=AC; maxBC=maxAC
    	, (*else*)
            BC= OPESimplify[OPE[B,C], simopts];
            maxBC = MaxPole[BC];
            BnAC[y_,n_] :=
    		 (BnAC[y,n]=OPESimplify[OPE[B, proj[y][OPEPole[n][AC]]], simopts])/;n<=maxAC;
                      BnAC[__] = OPEData[{}];
        ];
        Print["AnBC"];
        maxAnBC = Max[MaxPole /@ 
		Flatten[Table[AnBC[y,n],{n,maxBC},{y,0,3}]]];
        Print["BnAC"];
        maxBnAC = Max[MaxPole /@
		Flatten[Table[BnAC[y,n],{n,maxAC},{y,0,3}]]];
        Print["ABnC"];
        maxABnC = Max[MaxPole /@
		Flatten[Table[ABnC[y,n],{n,maxAB},{y,0,3}]]];
        maxn = Max[ maxAnBC, maxAC, maxAB ];
        maxm = Max[ maxBC, maxBnAC, maxABnC ];
        Print["equations"];
        Table[
            PoleSimplify[
               proj[x][OPEPole[n][AnBC[y,m]]] -
               (-1)^( OPEParity[thetas12[x]**A] OPEParity[thetas12[y]**B]) *
                   proj[y][OPEPole[m][BnAC[x,n]]] -
               (-1)^( OPEParity[thetas12[x]**A] OPEParity[thetas12[y]])*
                   Sum[ Binomial[n-1,p-1] proj[y][OPEPole[m+n-p][ABnC[x,p]]],
                     {p,1,n}
                   ],
              simopts],
            {m, maxm},{n,maxn},{x,0,3},{y,0,3}
        ]
]//Flatten
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Saving of intermediate results
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
ClearOPESavedValues[] := (
    Clear[
        OPECompositeHelpRQ,OPECompositeHelpLQ,
        NOCompositeHelpRQ,NOCompositeHelpLQ,
        OPECompositeHelpRPB,OPECompositeHelpLPB,
        NOCompositeHelpRPB,NOCompositeHelpLPB];
    RedefineOPECompositeHelpRQ[];
    RedefineOPECompositeHelpLQ[];
    RedefineNOCompositeHelpRQ[];
    RedefineNOCompositeHelpLQ[];
    RedefineOPECompositeHelpRPB[];
    RedefineOPECompositeHelpLPB[];
    RedefineNOCompositeHelpRPB[];
    RedefineNOCompositeHelpLPB[];
)
ClearOPESavedValues[]
(*
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
OPESave[name_String] := (
    If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0, Off[Unset::norep]];
    ClearOPECompositeHelpR[];
    ClearOPECompositeHelpL[];
    ClearNOCompositeHelpR[];
    ClearNOCompositeHelpL[];
    If[ !NumberQ[$VersionNumber] || $VersionNumber < 2.0, On[Unset::norep]];
    Put[ Definition[
            NOOrderingValue,
            OPECompositeHelpR, OPECompositeHelpL,
            NOCompositeHelpR, NOCompositeHelpL
         ],
         name];
    RedefineOPECompositeHelpR[];
    RedefineOPECompositeHelpL[];
    RedefineNOCompositeHelpR[];
    RedefineNOCompositeHelpL[]
)

(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Global options
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
SetAttributes[SetOPEOptions, HoldRest]

 (* 21-07-94 (OPEdefs 3.1 beta 1) changed definition of CallAndSave according
    to which option is set for OPESaving. Previously, all cases where
    handled through the last of the 3 definitions which follow.
    The first two are shortcuts for faster respons.
  *)
SetOPEOptions[OPESaving, True] :=
    ( Clear[CallAndSave]; CallAndSave[f_, arg__] := (f[arg] = f[arg]) )
SetOPEOptions[OPESaving, False] :=
    ( Clear[CallAndSave]; CallAndSave[f_, arg__] := f[arg] )
SetOPEOptions[OPESaving, expr_] :=
    If[!SameQ[expr,True] && !SameQ[expr, False],
        Message[OPESaving::nobool, HoldForm[expr]],
        (* else *)
        Clear[CallAndSave];
        CallAndSave[f_, arg__] :=
           If[ expr, f[arg] = f[arg], f[arg]]
    ]
OPESaving::nobool =
    "Error : `` is not an expression that evaluates to a boolean value.
OPESaving not changed.";

(*default*)
SetOPEOptions[OPESaving, True]
(*
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)

SetOPEOptions[NOOrdering, n_Integer]  := (NOOrderingValue = n)
NOOrderingValue = -1;
(*
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
SetOPEOptions[ParityMethod,0] :=
    (Message[ParityMethod::warning];
     SwapSign[A_,B_] := If[BosonQ[A] || BosonQ[B],1,-1];
     SwapSign[A_] := If[BosonQ[A],1,-1])
SetOPEOptions[ParityMethod,1] :=
    (SwapSign[A_,B_] := (-1)^(OPEParity[A] OPEParity[B]);
     SwapSign[A_] := (-1)^OPEParity[A])
ParityMethod::warning = "Warning : wrong results will occur if you use operators
which are declared with OPEOperator with a symbolic parity."

(*default SetOPEOptions[ParityMethod,0] is not called because of the message.*)
SwapSign[A_,B_] := If[BosonQ[A] || BosonQ[B],1,-1]
SwapSign[A_] := If[BosonQ[A],1,-1]
(*
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
SetOPEOptions[OPEMethod, QuantumOPEs] :=
    setQuantum[]
SetOPEOptions[OPEMethod, ClassicalOPEs] :=
    setPB[]

setQuantum[] := (
    ClearOPECompositeHelpR=ClearOPECompositeHelpRQ;
    ClearOPECompositeHelpL=ClearOPECompositeHelpLQ;
    ClearNOCompositeHelpR=ClearNOCompositeHelpRQ;
    ClearNOCompositeHelpL=ClearNOCompositeHelpLQ;
    RedefineOPECompositeHelpR=RedefineOPECompositeHelpRQ;
    RedefineOPECompositeHelpL=RedefineOPECompositeHelpLQ;
    RedefineNOCompositeHelpR=RedefineNOCompositeHelpRQ;
    RedefineNOCompositeHelpL=RedefineNOCompositeHelpLQ;
    NOCommuteHelp=NOCommuteHelpQ;
    OPEPoleHelpL=OPEPoleHelpLQ;
    OPEPoleHelpR=OPEPoleHelpRQ;
    OPECompositeHelpR=OPECompositeHelpRQ;
    OPECompositeHelpL=OPECompositeHelpLQ;
    NOCompositeHelpR=NOCompositeHelpRQ;
    NOCompositeHelpL=NOCompositeHelpLQ;
)

setPB[] := (
    ClearOPECompositeHelpR=ClearOPECompositeHelpRPB;
    ClearOPECompositeHelpL=ClearOPECompositeHelpLPB;
    ClearNOCompositeHelpR=ClearNOCompositeHelpRPB;
    ClearNOCompositeHelpL=ClearNOCompositeHelpLPB;
    RedefineOPECompositeHelpR=RedefineOPECompositeHelpRPB;
    RedefineOPECompositeHelpL=RedefineOPECompositeHelpLPB;
    RedefineNOCompositeHelpR=RedefineNOCompositeHelpRPB;
    RedefineNOCompositeHelpL=RedefineNOCompositeHelpLPB;
    NOCommuteHelp=NOCommuteHelpPB;
    OPEPoleHelpL=OPEPoleHelpLPB;
    OPEPoleHelpR=OPEPoleHelpRPB;
    OPECompositeHelpR=OPECompositeHelpRPB;
    OPECompositeHelpL=OPECompositeHelpLPB;
    NOCompositeHelpR=NOCompositeHelpRPB;
    NOCompositeHelpL=NOCompositeHelpLPB;
)
setPB[]
setQuantum[]
(*
:[font = input; initialization; nowordwrap; backColorRed = 65535; backColorGreen = 65535; backColorBlue = 65535; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; plain; fontName = "Courier New"; fontSize = 11; ]
*)
(* change 1.0 beta 2: added warning message *)
SetOPEOptions[EnableDummies,True] :=
(   Needs["Dummies`"];
    Dummies`Renumber[expr_OPEData, arg_] := 
        OPEMap[Dummies`Renumber[#,arg]&, expr];
    Dummies`NewDummies[expr_OPEData, arg_] := 
        OPEMap[Dummies`NewDummies[#,arg]&, expr];
    Dummies`DummySimplify[expr_OPEData, arg_] := 
        OPEMap[Dummies`DummySimplify[#,arg]&, expr];
    Dummies`SumDummy[expr_OPEData, arg__] := 
        OPEMap[Dummies`SumDummy[#,arg]&, expr];
    Dummies`SetDummiesOptions[Dummies`FunctionPattern, NO];
    SetOPEOptions[OPESaving,False];
    Message[SetOPEOptions::Dummies];
    Dummies`ResetDummies[]
)

SetOPEOptions::Dummies = "OPESaving switched off. Keep it off while doing
calculations with dummy indices.";
(*
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Protect[OPESaving, NOOrdering, SeriesArguments, OPEMethod, EnableDummies];

SetOPEOptions[op_, val_] := Message[SetOPEOptions::novalopt,op,val]
SetOPEOptions::novalopt =
    "Error : option `` with value `` not a valid combination.";
(*
:[font = subsection; inactive; preserveAspect; startGroup; Cclosed; ]
Output routines
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
OutputForm
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
Clear[formatstring]

 (* formatstring[n_,start_] makes a formatstring like
         "start|| `` ||start-1|| `` ... `` ||start-n+1|| ``"
    There are n ``'s to plug in n items.
 *)
formatstring[0,_] = " "
formatstring[start_, n_] := formatstring[start, n] =
    Block[{i},
        StringJoin[
             ToString[start],
             "|| ``",
             Sequence  @@ Flatten[
                 Table[{" ||",ToString[start+1-i],"|| ``"},{i,2,n}]
             ]
        ]
    ]

Format[OPEData[{a___}]] :=
    StringForm[
        StringJoin["<< ", formatstring[Length[{a}], Length[{a}]], " >>"],
        a
    ]

 (* 09/09/92
    OPEData written to file will now be in appropriate MakeOPE form.
    Note : using HoldForm leads to infinite recursion,
           SequenceForm works nice, but OutputForm of the strings is needed
           to discard the quotes in the output.
 *)
Format[OPEData[A_List],InputForm] :=
        SequenceForm[ OutputForm["MakeOPE["], A, OutputForm["]"] ]
(*
:[font = subsubsection; inactive; initialization; preserveAspect; startGroup; Cclosed; ]
TeXForm
:[font = input; initialization; preserveAspect; endGroup; endGroup; nowordwrap; ]
*)

(*Z2,Z12 are symbols used only for TeXForm *)
Clear[Z2,Z12]

Format[OPEData[{}], TeXForm] :=
   SequenceForm["{\\rm O}[",Z12,"]^0 "]

Format[OPEData[A_List], TeXForm] := 
  Block[{listop, listz, n, texlist},
     listop = (#[Z2]& /@ A) //.
        {(f_ +  g_) [x_] :> f[x] + g[x],
         (f_ *  g_) [x_] :> f * g[x] /; OperatorQ[g],
         Literal[NonCommutativeMultiply[theta__, g_] [x_]] :> 
             NonCommutativeMultiply[theta, g[x]],
         0[_] -> 0};
     listz = Table[{Z12, "^", -n}, {n,Length[A],1,-1}];
     texlist = {"\\left(", listop[[1]], "\\right)\,", listz[[1]]};
     Do[texlist = Join[texlist,
            {"\\ + \\ \\left(", listop[[n]], "\\right)\,", listz[[n]]}
        ], {n,2,Length[A]}
     ];
     texlist = Join[texlist,{"\ + \ {\\rm O}[",Z12,"]^0 "}]//Flatten;
     SequenceForm@@texlist
  ] 

Format[A_NonCommutativeMultiply,TeXForm] :=
   SequenceForm@@A
Format[t12, TeXForm] = "\\theta_{12} ";
Format[tb12, TeXForm] = "\\bar\\theta_{12} ";
Format[Z2, TeXForm] = "Z_2 ";
Format[Z12, TeXForm] = "Z_{12} ";
Format[DT[A_],TeXForm] := SequenceForm[TeXFormDT,A]
Format[DBT[A_],TeXForm] := SequenceForm[TeXFormDBT,A]
Format[DDB[A_],TeXForm] := SequenceForm[TeXFormDDB,A]
TeXFormDT = "D ";
TeXFormDBT = StringJoin["\\bar ",TeXFormDT];
TeXFormDDB = StringJoin["{[}", TeXFormDT,",",TeXFormDBT,"{]}"];
(*
:[font = input; initialization; preserveAspect; endGroup; nowordwrap; ]
*)
End[];
(*
:[font = input; initialization; preserveAspect; nowordwrap; ]
*)
EndPackage[]
(*
^*)