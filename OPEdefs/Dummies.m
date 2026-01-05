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
	fontset = message, "Courier New", 10, L0, nowordwrap, R32896;
	fontset = print, "Courier New", 10, L0, nowordwrap;
	fontset = info, "Courier New", 10, L0, nowordwrap;
	fontset = postscript, "Courier New", 8, L0, nowordwrap;
	fontset = name, "Times New Roman", 10, L0, nohscroll, italic, B32896;
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
(* :Name: Dummies` *)

(* :Title:  Summation over dummy indices *)

(* :Author: 
        Kris Thielemans
        Theoretical Physics Group
        Imperial College
        London SW7 2BZ
        kris@tfdec1.fys.kuleuven.ac.be
 *)

(* :Summary:
       This package is an attempt to implement summation over dummy indices.
       Some simple rules for the delta symbol are defined, and functions are 
       available which try to simplify expressions with dummy indices. 
       However, these rules don't take symmetries into account. 
*)

(* :Context: Dummies` *) 

(* :Package Version: 1.0 *)

(* :Mathematica Version: 2.0 *)
:[font = section; inactive; startGroup; Cclosed; ]
History
:[font = text; inactive; endGroup; ]
changes in beta 2 :
- added extra argument to DefineDummy
- changed Renumber and NewDummies to take into account that
  indices that occur only once (not summed over) could overlap
  with the new indices these functions assign. Change for Renumber 
  is efficient in time, but doesn't admit "holes" in the indices that are
  used.
beta 4:
- added Unravel to avoid a complete Expand in DummySimplify.
!!!
!!! Currently, NewDummies and Renumber are only safe on Expanded 
!!! expressions if indices are involved where no summation is implied.
!!!
beta 5:
- bug removed in SumDummy which prevented summing expressions without
  dummy indices.
beta 6:
- small bug removed in Unravel preventing spurious warning messages of Select.
- SumDummy still summed over dummy indices occuring only once.
beta 7: 
- removed support for OPEdefs (is now in OPEdefs itself)
- added SetDummiesOption[FunctionPattern,_] 
- replaced Unravel by Expand in Mathematica 2.3
beta 8:
- added ResetDummies[]
beta 9:
- Renumber, NewDummies and SumDummy were unsafe when
  expressions which are not fully expanded are used.
  Used Unravel now to prevent this problem. This is done by an
  option Unravel->True to these functions.
- Renumber now includes also the option Min -> number.
:[font = input; initialization; nowordwrap; ]
*)
BeginPackage["Dummies`", "Delta`"]
(*
:[font = input; initialization; nowordwrap; ]
*)
DummiesVersion = 1.0;
DummiesbetaVersion = 9;
Print["Dummies Version ", DummiesVersion,
      " (beta ",DummiesbetaVersion,") by Kris Thielemans"]

If[!NumberQ[$VersionNumber] || $VersionNumber < 2.0,
   Dummies::Version = "Dummies requires Mathematica 2.0 or later. Sorry.";
   Message[Dummies:Version]
]
(*
:[font = section; inactive; startGroup; Cclosed; ]
Usages
:[font = input; initialization; endGroup; nowordwrap; ]
*)
DefineDummy::usage = "DefineDummy[a] prepares the package for 
using dummy indices of the form a[something]. DefineDummy[a,int]
makes sure that new dummy indices will start from a[int].";

dummy::usage = "dummy[x] makes a new dummyindex of the form x[1].
DefineDummy[x] should be called first.";

Renumber::usage = "Renumber[expr, {x}] renumbers the \"dummy\"
arguments (which are of the form x[something]) in expr.
Only indices which occur at least twice are renumbered.\n
The second argument can also be a list. If it is not present, 
the list of declared dummies is used.\n
Renumber takes two options:\n
\t-Min -> value gives the starting value of the renumbered indices.
\t-Min -> NewDummies gives all indices unique values
\t-Unravel -> True specifies that Unravel[expr,x] should be called
first. Do not set to False if expr is not unraveled !\n
The second argument can also be a list. If it is not present, 
the list of declared dummies is used.";

dimension::usage = "dimension[k] is the length of the range over
which k[something] indices are summed.";

NewDummies::usage = "NewDummies[expr, {x}] gives the \"dummy\"
arguments (which are of the form x[something]) in expr a new
(unique) number. This is useful when e.g. multiplying identical
expressions.\n See Renumber for more details.";

DummySimplify::usage = "DummySimplify[expr, {x}] renumbers the \"dummy\"
arguments (which are of the form x[something]) in expr and applies
some rules for Delta.\n
The second argument can also be a list. If it is not present, 
the list of declared dummies is used.";

SumDummy::usage= "SumDummy[expr,{x,start,end}] explicitly sums the
dummy indices x[_] from start to end (where start defaults to 1).";

Unravel::usage = "Unravel[expr,pat]  is like Expand[expr,pat], but
much faster for big expressions.\n
It is also an option for Renumber, NewDummies, SumDummy";

SetDummiesOptions::usage = "A function to set global options for the 
Dummies` package. Current options : FunctionPattern."

FunctionPattern::usage = "SetDummiesOptions[FunctionPattern, pat_] can be
used to tell Dummies` which objects it should Hold in the process of
renumbering. The default\n
\tSetDummiesOptions[FunctionPattern,NonCommutativeMultiply]\n
for example makes sure that rules for NonCommutativeMultiply are not
applied in intermediate steps."

ResetDummies::usage = "ResetDummies[] clears previous definitions of
dummy indices."
(*
:[font = section; inactive; startGroup; ]
Implementation
:[font = input; initialization; nowordwrap; ]
*)
Begin["Dummies`Private`"];
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
SetDummiesOptions
:[font = input; initialization; endGroup; nowordwrap; ]
*)
SetDummiesOptions[FunctionPattern, pat_] :=
   (holdpattern = pat;)

SetDummiesOptions[FunctionPattern, NonCommutativeMultiply]
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
ResetDummies, DefineDummy, dummy
:[font = input; initialization; endGroup; nowordwrap; ]
*)
dummylist = {};
drule = {};

ResetDummies[] := 
  If[ !SameQ[dummylist, {}],
     Message[ResetDummies::reset];
     Delta[#[i_],#[i_]]=. & /@ dummylist;
     dummylist = {};
     drule = {};
  ]
ResetDummies::reset = "Clearing previous definitions for dummy indices.";

ResetDummies[] 

DefineDummy[a_,start_Integer:100] :=
  Block[{i},
     If[ MemberQ[dummylist, a], Return[]];	
	 dummyindex[a] = start;
	 AppendTo[dummylist,a];
     (* strange way of doing this because we have to be able to clear
        this rule in ResetDummies. Otherwise problems with i$_
     *)
	 (Delta[#[i_],#[i_]] = dimension[#]) & /@ {a};
	 AppendTo[drule, Delta[a[i_],aa_] x_ :> (x/.a[i]->aa) /;
		!FreeQ[x,a[i]]];
	 AppendTo[drule, Delta[a[i_],aa_]^2 :> Delta[aa,aa]];
  ]
 (* drule is to be applied to Expanded expressions only !
    Otherwise  a (x Delta[] + y) /. drule
    will not replace dummy indices in a ! 
  *)

dummy[a_] := a[dummyindex[a]++]
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
Unravel
:[font = input; initialization; endGroup; nowordwrap; ]
*)
Clear[Unravel,mult]
mult[a_,b_Plus] := a #&/@ b
mult[a_,b_] := a b

(* change 3/8/95 1.0 beta 7: forget about Unravel in new versions *) 
If[NumberQ[$VersionNumber] && $VersionNumber > 2.2,
   Unravel[expr_,pat_] := Expand[expr,pat],
   (* else *)
   Unravel[a_ b_,pat_] := mult[a, Unravel[b,pat]] /; 
	   FreeQ[a,pat];
   Unravel[a_Plus,pat_] := Unravel[#,pat]&/@a;

   (* change 3/8/95 1.0 beta 6 
      added the Return[] statement. In some cases, tmp is
      already simplified to a single term (for example 0).
      Without the Return, Select complained.
   *)
   Unravel[a_Times,pat_]:=
     Block[{tmp=Unravel[#,pat]&/@ a,tmp2},
        If[!SameQ[Head[tmp], Times], Return[tmp]];
        tmp2=Select[tmp,SameQ[Head[#],Plus]&&!FreeQ[#,pat]&];
        mult[tmp/tmp2,
   	       If[SameQ[Head[tmp2],Times],
   	          Outer[Times,Sequence@@tmp2],
   	          tmp2
   	       ]
        ] 
     ];
   Unravel[a_,_]:=a;
]
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
getDummies , getUnique
:[font = input; initialization; endGroup; nowordwrap; ]
*)
Clear[extractindices, getUnique] 

 (* extractindices[expr,arg] returns a list of all the 
    subexpressions of the form arg[_]
  *)
(* old, slow version
extractindices[arg_[i_],{___,arg_,___}]:={arg[i]}
 (*extractindices[a_,arg_]:= {} /;
	FreeQ[a,Alternatives@@arg] 
 *)
 (* change 28-06-93 : Flatten -> Join*)
extractindices[a_[b__],arg_] :=
    Join[extractindices[a,arg],
	 Sequence@@(extractindices[#,arg]& /@ {b})]
extractindices[_,_] := {}
*)
(* faster version, but not general
extractindices[expr_,arg_List] := Block[{ll={}},
   Scan[If[MatchQ[#,(Alternatives@@arg)[_]],
		AppendTo[ll,#]]&,expr,{-2}];
   ll
]
end of commented out version *)

(* actual implementation (fastest and general) *)
(* change 1.0 beta 7:
   removed next rule. 
   Purpose was to have NO's treated first in renumbering indices
   It didn't apply always.
   TODO: I should handle sums, lists, OPEData ,... also.
         Use of OperatorQ is problematic.

	testit=True;
	extractindices[a_ b_,arg_] :=
	    Block[{testit=False},
		Join[extractindices[a,arg],
		     extractindices[b,arg]]
	    ] /;
	    testit && !SameQ[Head[a],Times]&&OperatorQ[a];

*)
(* change 1.0 beta 7:
   introduced holdpattern instead of listing NO and others by hand.
   Hold is used to prevent reevaluation of complicated normal orderings, 
   noncommutative multiplications...
*)
extractindices[expr_,arg_List]:=
    Block[{l={}},
		expr /. {
		  x:holdpattern :> Hold[x],
		  x:(Alternatives@@arg)[_] :> (AppendTo[l,x];x)
		};
		l
   	]

 (* getUnique is like Union, but it preserves order 
  *)
getUnique[l_List] :=
    Block[{l2={}},
       If[!MemberQ[l2,#],AppendTo[l2,#]]& /@ l;
       l2
    ]
 (* getDummies[expr, arg] returns a list with 2 elements : 
    first the dummies that occur more than once,
    second the dummies that occur only once.
    It is important that the first list contains indices
    occuring in the same order as in expr.
    !! getDummies works only correctly on Unraveled expressions
  *)
getDummies[expr_,arg_] :=
    Block[{ind = extractindices[expr,arg], unind},
	unind = getUnique[ind];
	ind = Select[unind, Count[ind,#]>1&];
	{ind, Complement[unind,ind]}
    ]
(* an attempt to define getDummies on 'raveled' expressions.
   not debuged.

getDummies[a_Times,arg_] :=
   Block[{indices,onlyOnce},
      indices = getDummies[List@@a,arg];
      If[MatchQ[indices, {_,{}}],
         Return[{getUnique[First/@indices],{}}]];
      
      onlyOnce = Join@@Transpose[indices][[2]];
      indices = Apply[
         Function[{summed, nosum},
            Join[summed, Select[nosum, Count[onlyOnce,#]>1&]]
         ],
         indices,
         1
      ];
      {indices, Complement[onlyOnce, indices]}
   ]
*)
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
Renumber, NewDummies
:[font = input; initialization; endGroup; nowordwrap; ]
*)
Clear[Renumber, NewDummies]

 (* Renumber renumbers all summation indices starting 
    from 1.
    change 1.0 beta 9:
    Now always require a List as second argument.
  *)
Options[Renumber] = {Min -> 1, Unravel->True};

Renumber[expr_,opts___Rule] :=
    Renumber[expr,dummylist,opts]
Renumber[expr_, {},opts___Rule] :=
    expr

Renumber[expr_List, arg_List,opts___Rule] :=
    Renumber[#,arg,opts]& /@ expr
Renumber[expr_Plus, arg_List,opts___Rule] :=
    Renumber[#,arg,opts]& /@ expr
Renumber[expr_, {arg1_,args___},opts___Rule] :=
	Renumber[
	   Renumber[
	      Unravel[expr,Alternatives[arg1,args]],
	      {arg1},Unravel->False,opts
	   ], {args},Unravel->False,opts
	] /; (Unravel /. {opts}/. Options[Renumber])

Renumber[expr_,{arg_},opts___Rule]:=
    expr /; FreeQ[expr,arg]
Renumber[expr_,{arg_},opts1___Rule,Unravel->False,opts2___Rule] :=
  Block[{ind,nosum,newind,max, intersect, start},
    start = Min /. {opts1,opts2} /. Options[Renumber];
    {ind,nosum}=getDummies[expr,{arg}];
    max=Length[ind];
	If[ start === NewDummies,
	    start = dummyindex[arg];
	    dummyindex[arg] += max
	];
	newind=arg/@ Range[start,start+max-1];
	
	While[Length[intersect=Intersection[newind,nosum]]>0,
	    (* if Min->NewDummies one should never get here *)
	    newind = newind /. 
		Thread[intersect->Map[max+#&,intersect,{2}]];
		(* make sure the indices of newind are
		   smaller than dummyindex[arg] *)
		dummyindex[arg] += max;
	];
    expr /. Thread[ind->Sort[newind]]
 ]

NewDummies[expr_] := NewDummies[expr,dummylist]
NewDummies[expr_,arg_] := Renumber[expr,arg,Min->NewDummies]
(* change 1.0 beta 9:
   modified Renumber to include NewDummies

NewDummies[expr_, {}] := expr
NewDummies[expr_List, arg_] := 
     Map[NewDummies[#,arg]&, expr]
NewDummies[expr_Plus,arg_] := NewDummies[#,arg]& /@ expr

NewDummies[expr_, {arg1_,args__}] :=
	NewDummies[ NewDummies[expr,{arg1}], {args}]

NewDummies[expr_,{arg_}] :=
     Block[{ind, nosum, newind, max},
        {ind,nosum}=getDummies[expr,{arg}];
	dummyindex[arg] += (max = Length[ind]);
	newind=arg/@ 
		Range[dummyindex[arg]-max,dummyindex[arg]-1];
	While[Length[intersect=Intersection[newind,nosum]]>0,
	     Print["NewDummies:shouldn't be here. expr :",expr];
	     newind = newind /. 
		Thread[intersect->Map[max+#&,intersect,{2}]];
		(* make sure the indices of newind are
		   smaller than dummyindex[arg] *)
		dummyindex[arg] += max;

	];
        expr /. Thread[ind->Sort[newind]]
     ]
*)
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
DummySimplify
:[font = input; initialization; endGroup; nowordwrap; ]
*)
DummySimplify[expr_] := DummySimplify[expr,dummylist]

(* change 1.0 beta 9: 
   - now always requires a List as second argument
   - added next rule
 *)
DummySimplify[expr_, {}] := expr
DummySimplify[expr_List, arg_List] := 
     Map[DummySimplify[#,arg]&, expr]

(* change 1.0 beta 9: added Unravel->False *)
DummySimplify[expr_,arg_] :=
     Renumber[
        Unravel[expr,Alternatives@@dummylist]//.drule,
        arg, Unravel-> False
     ]
(*
:[font = subsection; inactive; startGroup; Cclosed; ]
SumDummy
:[font = input; initialization; endGroup; nowordwrap; ]
*)
Clear[SumDummy]
(* change 1.0 beta 9 : added Unravel option*)
Options[SumDummy] = {Unravel -> True};

SumDummy[expr_List, arg__] := 
     Map[SumDummy[#,arg]&, expr]
SumDummy[expr_Plus, arg__] := 
     Map[SumDummy[#,arg]&, expr]
(* change 1.0 beta 9 : added this trivial case *)
SumDummy[expr_, opts___Rule] := expr
SumDummy[expr_,ss:({_,__Integer}..)] :=
   SumDummy[expr, ss, Sequence@@Options[SumDummy]]
SumDummy[expr_,ss:({_,__Integer}..), Unravel -> True ] :=
     SumDummy[
         Unravel[expr,Apply[Alternatives,First/@ {ss}]],
         ss,
         Unravel -> False
     ]
     
makeDimRule[{k_,start_:1,end_}] := dimension[k]->end-start+1
makeSumRule[indices_List] :=
	{Head[#]->#}& /@ indices
(* change 1.0 beta 5 :
   added test to see if no dummies are present
   change 1.0 beta 6 :
   replace extractindices by First[getDummies]. This to prevent
   summing over indices which occur only once.
 *)
SumDummy[exp_,ss:({_,__Integer}..), Unravel-> False] :=
  Block[
	{ind = Union[First[getDummies[exp,{First[#]}]]]& /@ {ss}},
    If[ SameQ[ind, {{}}],
       Return[exp //. Map[makeDimRule,{ss}]]
    ];
	Sum[Release[exp //. Map[makeDimRule,{ss}]], 
	    Release[Sequence @@
		    Flatten[MapThread[ReplaceAll,
				{{ss},makeSumRule/@ ind}
		    ],1]
	    ]
	]
  ]
(*
:[font = input; initialization; endGroup; nowordwrap; ]
*)
End[];
(*
:[font = input; initialization; nowordwrap; ]
*)
EndPackage[]
(*
^*)