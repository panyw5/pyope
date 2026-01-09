#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter::Closed:: *)
(*READ ME*)


(* ::Text:: *)
(*In order for this notebook to evaluate properly, please download the packages specified in the "Load Packages" section.*)


(* ::Text:: *)
(*We begin by giving a quick run through on how to generate vacuum characters of certain VOAs with this code. We will also give advice on how to generalise the code to other VOAs after some sample calculations.*)


(* ::Text:: *)
(*NOTE: We recommend restarting the kernel when swapping between different VOA calculations. This is due to the fact that the codes share definitions of variables which may cause errors.*)


(* ::Subsection:: *)
(*User Guide*)


(* ::Text:: *)
(*Here we give a list of steps and some examples on how to run the code.*)


(* ::Text:: *)
(*1. Open the "Load Packages" chapter and make sure that the package requirements listed within are loaded correctly by Mathematica. We included a set directory command just in case the packages need to be navigated to. *)
(**)
(*NOTE: The first package should be downloaded and saved as "ope.math" but called with "<<ope.math.m" as seen below.*)
(**)
(*This is what it should look like when they are loaded correctly:*)


(* ::Input:: *)
(*<<ope.math.m*)


(* ::Input:: *)
(*<<OPEDefs`*)


(* ::Text:: *)
(*2. Open up the "Settings" chapter. Here one can change a pair of parameters the notebook takes when running the calculations.*)


(* ::Input:: *)
(*fcal=Checkbox[False];*)


(* ::Text:: *)
(*Ticking the above box  speeds up the code by not searching for the degeneracies of the negatively charged-states (under the VOA symmetries). Calculating that part of the spectrum is very intensive, but it can be gotten from the positively-charged states by conjugation."Checkbox[False]" means it is turned off, while "Checkbox[True]" means it is on. The conservative reader can leave the box unticked.*)


(* ::Input:: *)
(*lev=;*)


(* ::Text:: *)
(*This parameter controls the weight at which the expansion will be calculated. This has to be specified for the code to be able to run. To get the vacuum character up to a given order, the code needs to be evaluated successively for different values of lev.*)
(**)
(*NOTE: This number can be integer or half-integer only.*)


(* ::Text:: *)
(*3. Pick a specific VOA. There are 4 copies of very similar code in this notebook; everything which will be said from here on is relevant for all cases. *)
(**)
(*Expand the title related to the relevant VOA and run the entirety of the "Definitions" section. This will define the variables used for the OPEs (bosonic and fermionic), as well as the definitions of the abstract VOA generators in terms of the \[Beta]\[Gamma]bc system. There is no user input required at this point.*)


(* ::Text:: *)
(*4. Open the "Character Calculation" chapter and run the entirety of the "Field Generation" section. This will generate all possible configurations built out of the strong generators and their derivatives (given in the earlier sections in terms of free fields) at the specific weight specified by the user in lev. No user input is required at this stage.*)


(* ::Text:: *)
(*5. Run the entirety of the "Selection via Quantum Numbers" section. This section of the code grades the previously generated operators by the symmetries of the N=2 VOA. No user input is required at this stage.*)


(* ::Text:: *)
(*6. Run the entirety of the "Pre-Processing" section. No user input is required at this stage.*)


(* ::Text:: *)
(*7. Choose whether to run the most refined calculation, or save computational power and do an unrefined calculation. The former includes an R-filtration grading, while the latter one does not. The codes for calculating both are given in the "Calculating Physical States" section, given by "R-Filtration Level-by-Level Character Calculator" and "Level-by-Level Character Calculator" subsections respectively. Evaluating the required subsubsection computes the vacuum character at the requested level.*)
(**)
(*NOTE: At this stage the code determines the degeneracy of null states for fixed quantum numbers. The result is then subtracted from the naive degeneracy of states.*)


(* ::Text:: *)
(*8. An example for the output of the above procedure in the refined case, and for the VOA corresponding to the complex reflection group Subscript[Z, 3] at lev = 3, is the following:*)


(* ::Input:: *)
(*rcharacter=Quiet[RFullCharCalculator[If[fcal===Checkbox[False],fullgrouping,sfullgrouping]]/.{MatrixRank[{}]->0}]*)


(* ::Text:: *)
(*NOTE: If the fcal option was selected, the above degeneracies will only include non-negatively charged states (under the red gl(1) symmetry). This result shows the list of all the states along with their degeneracies and all the quantum numbers. The bold number is the degeneracy; inside the square brackets are the quantum numbers: red is for m (gl(1) in osp(2|2)), while blue is for r (gl(1) outer automorphism).*)


(* ::Input:: *)
(*rfcharacter=If[fcal===Checkbox[False],Quiet[Sum[fermlist[[i,j]]* Head[rcharacter[[i,j,k]]][[1]] z^rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^rcharacter[[i,j,k,1,2,1]],{i,Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]],Quiet[Sum[sfermlist[[i,j]]*Head[rcharacter[[i,j,k]]][[1]] z^rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^rcharacter[[i,j,k,1,2,1]],{i,Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]+Sum[sfermlist[[i,j]]*Head[rcharacter[[i,j,k]]][[1]] z^-rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^-rcharacter[[i,j,k,1,2,1]],{i,If[(IntegerQ[lev]||MemberQ[mparam,0]),2,1],Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]]]*)


(* ::Text:: *)
(*In the above line we are expressing the refined vacuum character at a certain order in a q expansion (in this case q^3). This is evaluated by a signed sum of monomials in the fugacities, raised to the power of the corresponding charges of the states and multiplied by the degeneracies. The v and z are fugacities  associated with the gl(1) outer automorphism and gl(1) subalgebra of the 2D N = 2 SCA respectively, while y with the weight under the R-filtration.*)


(* ::Text:: *)
(*9. The output in the unrefined case is  similar to that of the previous step, up to the R-filtration grading, and given here for the Subscript[Z, 3]case at lev =3 :*)


(* ::Input:: *)
(*character=Quiet[FullCharCalculator[If[fcal===Checkbox[False],fullgrouping,sfullgrouping]]/. MatrixRank[{}]->0]*)


(* ::Input:: *)
(*fcharacter=If[fcal===Checkbox[False],Quiet[Sum[fermlist[[i,j]]*Head[character[[i,j]]][[1]] z^character[[i,j,1,1,1]] v^character[[i,j,1,2,1]],{i,Length[character]},{j,Length[character[[i]]]}]/. MatrixRank[{}]->0],Quiet[Sum[sfermlist[[i,j]]*Head[character[[i,j]]][[1]] z^character[[i,j,1,1,1]] v^character[[i,j,1,2,1]],{i,Length[character]},{j,Length[character[[i]]]}]+Sum[sfermlist[[i,j]]*Head[character[[i,j]]][[1]] z^-character[[i,j,1,1,1]] v^-character[[i,j,1,2,1]],{i,If[(IntegerQ[lev]||MemberQ[mparam,0]),2,1],Length[character]},{j,Length[character[[i]]]}]/. MatrixRank[{}]->0]]*)


(* ::Text:: *)
(*NOTE: We recommend restarting the kernel when swapping between different VOA. This is due to the fact that the codes share definitions of variables which may cause errors.*)


(* ::Subsection::Closed:: *)
(*Advice on Generalisations*)


(* ::Text:: *)
(*In this section we will be discussing how the reader can generalise the code for applications to VOAs not discussed in the accompanying paper. Below we will list the relevant points where modifications would have to be made:*)


(* ::Text:: *)
(*Definitions:*)
(**)
(*opfields - (add the new operators)*)
(*fields - (add the new operators)*)
(*Bosonic[...] - (add as many bc\[Beta]\[Gamma] fields as ranks + add the new operators)*)
(*Fermionic[...] - (add as many bc\[Beta]\[Gamma] fields as ranks + add the new operators)*)
(*OPE[...] - (add additional OPEs between new bc\[Beta]\[Gamma] fields)*)
(*the list of nested "Derivativers" - (add the new operators)*)
(*list of operators at the end - (add the new operators in free field representation, REQUIRES KNOWLEDGE OF REPRESENTATION)*)
(**)
(*Field Generation//Definitions:*)
(**)
(*gen[...] - (add the new operators + change weights of operators)*)
(**)
(*Selection via Quantum Numbers//Initial Quantum Grouping:*)
(**)
(*list of operators right under title each beginning with d (eg: dww[...]) - (add new operators)*)
(*fulllist - (add the new operators)*)
(*mnumber - (add the new operators + change quantum numbers)*)
(*rnumber - (add the new operators + change quantum numbers)*)
(*rfilnumber - (add the new operators + change quantum numbers)*)
(*fnumber - (add the new operators + change quantum numbers)*)
(**)
(*Calculating the Physical States//R-Filtration Level-by-Level Character Calculator:*)
(**)
(*rcounter - (add the new operators in the substitution rule)*)
(*RSingleCharCalculator - (add the new operators in the substitution rule)*)
(*rfcharacter - (change the grading or representation of the character in any way)*)
(**)
(*Calculating the Physical States//Level-by-Level Character Calculator:*)
(**)
(*SingleCharCalculator - (add the new operators in the substitution rule)*)
(*fcharacter - (change the grading or representation of the character in any way)*)


(* ::Chapter:: *)
(*Load Packages*)


(* ::Text:: *)
(*In this section the code will import the necessary packages. The two packages required are:*)
(*- OPEdefs (download + documentation: https://people.maths.ox.ac.uk/beem/voaschool.html )*)
(*- ope.math (download + documentation: https://www.u-aizu.ac.jp/~a-fujitu/mathematica/ )*)


(* ::Input:: *)
(*SetDirectory[""]*)


(* ::Input:: *)
(*<<ope.math.m*)


(* ::Input:: *)
(*<<OPEDefs`*)


(* ::Chapter:: *)
(*Settings*)


(* ::Text:: *)
(*For speeding up the calculation, by cutting out the negative quantum numbers. (Due to symmetry of the spectrum under conjugation)*)


(* ::Input:: *)
(*fcal=Checkbox[True];*)


(* ::Text:: *)
(*Enter the level/weight of the desired calculation.*)


(* ::Input:: *)
(*lev=4;*)


(* ::Title:: *)
(*Subscript[Z, 3] VOA *)


(* ::Chapter::Closed:: *)
(*Subscript[Z, 3] - Definitions*)


(* ::Text:: *)
(*In this section we define all the operators being used in our code. Both the free fields (b, c, \[Beta], \[Gamma]) and the abstract algebra strong generators. (In order to generalise the code to higher ranks or other algebras, one must change and/or add the generators appearing in the new algebra in the definitions below. We will highlight where change is needed throughout the code with comments CNH = Change Needed Here. The comment should be understood for multiple lines. There is more information about this under the "READ ME//Advice on Generalisation" section.)*)


(* ::Input:: *)
(*opfields={ww[0,1,z],jj0[0,1,z],wwb[0,1,z],tt[0,1,z],gg[0,1,z],ggt[0,1,z],ggw[0,1,z],ggwb[0,1,z]};(*CNH*)*)
(*fields={ww,jj0,wwb,tt,gg,ggt,ggw,ggwb};*)


(* ::Input:: *)
(*Bosonic[\[Beta]1,\[Gamma]1,ww[i_,j_,z],jj0[i_,j_,z],wwb[i_,j_,z],tt[i_,j_,z]];(*CNH*)*)
(*Fermionic[b1,c1,gg[i_,j_,z],ggt[i_,j_,z],ggw[i_,j_,z],ggwb[i_,j_,z]];*)


(* ::Input:: *)
(*OPE[\[Beta]1,\[Gamma]1]=MakeOPE[-One (z-w)^-1+Ord[z,w,0]];*)
(*OPE[b1,c1]=MakeOPE[One (z-w)^-1+Ord[z,w,0]];*)


(* ::Input:: *)
(*Derivativer[x_]:=x';*)


(* ::Input:: *)
(*wwww[i_,j_,z_]:=Nest[Derivativer,w,i];*)
(*j0j0[i_,j_,z_]:=Nest[Derivativer,j0,i];*)
(*wbwb[i_,j_,z_]:=Nest[Derivativer,wb,i];*)
(*tttt[i_,j_,z_]:=Nest[Derivativer,t,i];(*CNH*)*)
(*gggg[i_,j_,z_]:=Nest[Derivativer,g,i];*)
(*gtgt[i_,j_,z_]:=Nest[Derivativer,gt,i];*)
(*gwgw[i_,j_,z_]:=Nest[Derivativer,gw,i];*)
(*gwbgwb[i_,j_,z_]:=Nest[Derivativer,gwb,i];*)


(* ::Text:: *)
(*Below is where the abstract generators are written in terms of the free-field realisation.*)


(* ::Input:: *)
(*w=\[Beta]1;*)
(*j0=2 NO[b1,c1]+3 NO[\[Beta]1,\[Gamma]1];*)
(*wb=NO[\[Beta]1,NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,\[Gamma]1]]]]+2 NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,NO[b1,c1]]]]-4 NO[\[Beta]1,NO[Derivative[1][\[Gamma]1],\[Gamma]1]]-4/3 NO[\[Gamma]1,NO[b1,Derivative[1][c1]]]+2/3 NO[\[Gamma]1,NO[Derivative[1][b1],c1]]+2/3 NO[Derivative[1][\[Beta]1],NO[\[Gamma]1,\[Gamma]1]]-8/3 NO[Derivative[1][\[Gamma]1],NO[b1,c1]]+(10 \[Gamma]1^\[Prime]\[Prime])/9;*)
(*t=-2 NO[b1,Derivative[1][c1]]-3/2 NO[\[Beta]1,Derivative[1][\[Gamma]1]]-NO[Derivative[1][b1],c1]-1/2 NO[Derivative[1][\[Beta]1],\[Gamma]1];*)
(*g=NO[\[Gamma]1,b1];*)
(*gt=2NO[\[Beta]1',c1]+3NO[\[Beta]1,c1'];*)
(*gw=b1;*)
(*gwb=8/3 NO[b1,NO[c1^\[Prime]\[Prime],c1]]+3 NO[\[Beta]1,NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,Derivative[1][c1]]]]]-4 NO[\[Beta]1,NO[\[Gamma]1,NO[b1,NO[Derivative[1][c1],c1]]]]-4 NO[\[Beta]1,NO[\[Gamma]1,c1^\[Prime]\[Prime]]]-4 NO[\[Beta]1,NO[Derivative[1][\[Gamma]1],Derivative[1][c1]]]-2/3 NO[Derivative[1][b1],NO[Derivative[1][c1],c1]]+2 NO[Derivative[1][\[Beta]1],NO[\[Beta]1,NO[\[Gamma]1,NO[\[Gamma]1,c1]]]]-8/3 NO[Derivative[1][\[Beta]1],NO[Derivative[1][\[Gamma]1],c1]]+2/3 NO[\[Beta]1^\[Prime]\[Prime],NO[\[Gamma]1,c1]]+(10 \!\(\*SuperscriptBox[\(c1\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\))/9;*)


(* ::Input:: *)
(*(*CNH*)*)


(* ::Subsection:: *)
(*Functions*)


(* ::Input:: *)
(*NO[x_]:=x;*)


(* ::Input:: *)
(*pwrexpand[x_]:=Module[*)
(*{y=x,list={},i,j,p=Variables[x]},*)
(*For[i=1,i<Length[p]+1,i++,For[j=1,j<(Exponent[y,p[[i]]]+1),j++,list=Append[list,Head[p[[i]]]]]];*)
(*list*)
(*]*)


(* ::Input:: *)
(*parorder[x_]:=Block[*)
(*{y=x,l1,r},*)
(*l1=Table[NO@@y[[i]],{i,Length[y]}];*)
(*r=Total[l1];*)
(*r*)
(*]*)


(* ::Input:: *)
(*parordermod[x_]:=Block[*)
(*{y=x,r},*)
(*r=NO@@y;*)
(*r*)
(*]*)


(* ::Chapter:: *)
(*Subscript[Z, 3] - Character Calculation*)


(* ::Section::Closed:: *)
(*Field Generation*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*In this section of the code we are defining the rules on how to generate the available fields at each level. This is done by taking into account the weight of each field and whether said weight is integer or odd half-integer.*)


(* ::Input:: *)
(*gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===ww[0,1,z]):=ww[gr-3/2,1,z];*)
(*gen[op_, gr_]/;(gr>=1&&IntegerQ[gr]&&op===jj0[0,1,z]):=jj0[gr-1,1,z];*)
(*gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===wwb[0,1,z]):=wwb[gr-3/2,1,z];*)
(*gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===tt[0,1,z]):=tt[gr-2,1,z];*)
(*gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===gg[0,1,z]):=gg[gr-3/2,1,z];*)
(*gen[op_, gr_]/;(gr>=3/2&&OddQ[2gr]&&op===ggt[0,1,z]):=ggt[gr-3/2,1,z];*)
(*gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===ggw[0,1,z]):=ggw[gr-2,1,z];*)
(*gen[op_, gr_]/;(gr>=2&&IntegerQ[gr]&&op===ggwb[0,1,z]):=ggwb[gr-2,1,z];(*CNH*)*)
(*gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===ww[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;(!IntegerQ[gr]&&op===jj0[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===wwb[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===tt[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===gg[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!OddQ[2gr]||gr<3/2)&&op===ggt[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===ggw[0,1,z]):=Nothing;*)
(*gen[op_, gr_]/;((!IntegerQ[gr]||gr<2)&&op===ggwb[0,1,z]):=Nothing;*)


(* ::Input:: *)
(*OperatorList[oplis_,grlis_]:=Table[gen[op, gr],{gr,grlis},{op,oplis}]*)


(* ::Input:: *)
(*OperatorMult[a_,b_,right__]:=NO[a,b,right];*)
(*OperatorMult[a_,b_]:=NO[a,b];*)
(*OperatorMult[a_]:=NO[a];*)


(* ::Input:: *)
(*GenOp[partlis_]:=Flatten@Table[Outer[OperatorMult,Sequence@@OperatorList[opfields,partition]],{partition,partlis}]*)


(* ::Subsection::Closed:: *)
(*Level*)


(* ::Text:: *)
(*In this section we employ the definitions set in the previous sections and generate all the possible field combinations at a given level. This is done by generating all integer partitions and then allocating them to specific fields. The answer we end up with is a list of normal ordered products of fields and their derivatives.*)


(* ::Input:: *)
(*part=1/2IntegerPartitions[2 lev]*)


(* ::Input:: *)
(*generated=Select[GenOp[part], !NumberQ[#]&]/.-1*x_:>x;*)
(*DeleteDuplicates[%];*)
(*finalgen=Table[If[Head[%[[i]]]===NO,ReplaceAll[%[[i]],NO->Nl],Nl[%[[i]]]],{i,Length[%]}]*)


(* ::Section::Closed:: *)
(*Selection via Quantum Numbers*)


(* ::Text:: *)
(*In this section we will be solely focusing on grouping the previously generated field combinations with respect to their quantum numbers. We get the quantum numbers from the prescription in Table 1 of the accompanying paper.*)


(* ::Subsection:: *)
(*Initial Quantum Grouping*)


(* ::Input:: *)
(*dww[x_,y_,z_]:=ww;*)
(*djj0[x_,y_,z_]:=jj0;*)
(*dwwb[x_,y_,z_]:=wwb;*)
(*dtt[x_,y_,z_]:=tt;*)
(*dgg[x_,y_,z_]:=gg;(*CNH*)*)
(*dggt[x_,y_,z_]:=ggt;*)
(*dggw[x_,y_,z_]:=ggw;*)
(*dggwb[x_,y_,z_]:=ggwb;*)


(* ::Input:: *)
(*finalgen/.Nl->List;*)
(*Times@@#&/@%;*)
(*fulllist=%/.{ww->dww,jj0->djj0,wwb->dwwb,tt->dtt,gg->dgg,ggt->dggt,ggw->dggw,ggwb->dggwb}(*CNH*)*)


(* ::Input:: *)
(*mnumber[x_]:=Association[{ww->3/2,jj0->0,wwb->-(3/2),tt->0,gg->-(1/2),ggt->1/2,ggw->1,ggwb->-1}]@x*)
(*rnumber[x_]:=Association[{ww->0,jj0->0,wwb->0,tt->0,gg->1/2,ggt->-(1/2),ggw->1/2,ggwb->-(1/2)}]@x(*CNH*)*)
(*rfilnumber[x_]:=Association[{ww->3/2,jj0->1,wwb->3/2,tt->1,gg->1,ggt->1,ggw->3/2,ggwb->3/2}]@x*)
(*fnumber[x_]:=Association[{ww->0,jj0->0,wwb->0,tt->0,gg->1,ggt->1,ggw->1,ggwb->1}]@x*)


(* ::Input:: *)
(*mtotal[x_]:=Module[*)
(*{y=x,m=0,i=1},*)
(*For[i=1,i<Length[fields]+1,i++,m=m+Exponent[y,fields[[i]]]mnumber[fields[[i]]]];*)
(*m*)
(*]*)


(* ::Input:: *)
(*rtotal[x_]:=Module[*)
(*{y=x,r=0,i=1},*)
(*For[i=1,i<Length[fields]+1,i++,r=r+Exponent[y,fields[[i]]]rnumber[fields[[i]]]];*)
(*r*)
(*]*)


(* ::Input:: *)
(*rfiltotal[x_]:=Module[*)
(*{y=x,r=0,i=1},*)
(*For[i=1,i<Length[fields]+1,i++,r=r+Exponent[y,fields[[i]]]rfilnumber[fields[[i]]]];*)
(*r*)
(*]*)


(* ::Input:: *)
(*ftotal[x_]:=Module[*)
(*{y=x,r=0,i=1},*)
(*For[i=1,i<Length[fields]+1,i++,r=r+Exponent[y,fields[[i]]]fnumber[fields[[i]]]];*)
(*r*)
(*]*)


(* ::Input:: *)
(*mtset=mtotal[#]&/@fulllist;*)


(* ::Input:: *)
(*(*Here we will begin by grouping by the quantum number m first.*)*)


(* ::Input:: *)
(*mgroupstates=Table[Append[{},finalgen[[#]]]&/@Position[mtset,i]//Flatten,{i,Min[mtset],Max[mtset],1/2}]/.{}->Nothing*)


(* ::Input:: *)
(*mgfulllist=Table[Append[{},fulllist[[#]]]&/@Position[mtset,i]//Flatten,{i,Min[mtset],Max[mtset],1/2}]/.{}->Nothing*)


(* ::Input:: *)
(*rtset=Table[Table[rtotal[mgfulllist[[i,j]]],{j,Length[mgfulllist[[i]]]}],{i,Length[mgfulllist]}]*)


(* ::Input:: *)
(*(*Then we will introduce a new grading via the quantum number r.*)*)


(* ::Input:: *)
(*fullgrouping=Table[Table[Append[{},mgroupstates[[j,#]]]&/@Position[rtset[[j]],i]//Flatten,{i,Min[rtset[[j]]],Max[rtset[[j]]],1/2}],{j,Length[rtset]}]/.{}->Nothing*)


(* ::Input:: *)
(*fulllistgrouping=Table[Table[Append[{},mgfulllist[[j,#]]]&/@Position[rtset[[j]],i]//Flatten,{i,Min[rtset[[j]]],Max[rtset[[j]]],1/2}],{j,Length[rtset]}]/.{}->Nothing*)


(* ::Input:: *)
(*fermlist=Table[Table[(-1)^ftotal[fulllistgrouping[[i,j,1]]],{j,Length[fulllistgrouping[[i]]]}],{i,Length[fulllistgrouping]}]*)


(* ::Subsection:: *)
(*Complicated Quantum Grouping*)


(* ::Input:: *)
(*mlbl=Table[If[finalgen[[#]]&/@Flatten[Position[mtset,i]]=={},Nothing,i],{i,Min[mtset],Max[mtset],1/2}]*)


(* ::Input:: *)
(*ttt=Table[Append[{},{finalgen[[#]],i}]&/@Position[mtset,i]//Flatten,{i,Min[mtset],Max[mtset],1/2}]/.{}->Nothing*)


(* ::Input:: *)
(*Table[Table[Append[{},{mgroupstates[[j,#]],{mlbl[[j]],i}}]&/@Flatten[Position[rtset[[j]],i]],{i,Min[rtset[[j]]],Max[rtset[[j]]],1/2}],{j,Length[rtset]}]/.{}->Nothing;*)
(*y1=%*)


(* ::Input:: *)
(*flt1[x_]:=*)
(*Module[*)
(*{uuu=x,k,l,j,i},*)
(*For[l=1,l<Length[uuu]+1,l++,For[k=1,k<Length[uuu[[l]]]+1,k++,For[j=1,j<Length[uuu[[l,k]]]+1,j++,For[i=1,Length[uuu[[l,k,j]]]==1,i++,uuu[[l,k,j]]=Flatten[uuu[[l,k,j]],1]]]]];*)
(*uuu*)
(*]*)


(* ::Input:: *)
(*nothinger[x_,right___]:=Nothing;*)


(* ::Input:: *)
(*ioi=flt1[y1]*)


(* ::Input:: *)
(*ioi/.Nl->nothinger;*)
(*Table[Table[DeleteDuplicates[%[[i,j]]],{j,Length[%[[i]]]}],{i,Length[%]}];*)
(*Table[Table[Flatten[%[[i,j]],2],{j,Length[%[[i]]]}],{i,Length[%]}];*)
(*fcollectioncoloured=Table[Table[Table[If[k==1,Style[%[[i,j,k]],Red],Style[%[[i,j,k]],Blue]],{k,1,2}],{j,Length[%[[i]]]}],{i,Length[%]}]*)


(* ::Section::Closed:: *)
(*Pre-Processing*)


(* ::Text:: *)
(*In this section we set up certain calculational parameters to ease the character calculation later on. This is especially important if the faster calculation (fcal) option has been ticked in the settings.*)


(* ::Input:: *)
(*mparam=mtset//Union*)


(* ::Input:: *)
(*modstart=If[EvenQ[Length[mparam]],Length[mparam]/2,(Length[mparam]-1)/2+1]*)


(* ::Input:: *)
(*ansmodstart=If[EvenQ[Length[mparam]],Length[mparam]/2+1,(Length[mparam]-1)/2+1]*)


(* ::Input:: *)
(*sfullgrouping=Take[fullgrouping,-modstart];*)
(*sfulllistgrouping=Take[fulllistgrouping,-modstart];*)
(*sfcollectioncoloured=Take[fcollectioncoloured,-modstart];*)
(*sfermlist=Take[fermlist,-modstart];*)


(* ::Section::Closed:: *)
(*Calculating the Physical States*)


(* ::Text:: *)
(*In this section we give two separate calculations for the vacuum character of the VOA. The first gives a more refined character with 3 fugacities, while the second only 2 fugacities (their content is elaborated upon in the individual sections). It is important to note that the character is calculated on a level-by-level basis which has been set in the settings tab.*)


(* ::Subsection::Closed:: *)
(*R-Filtration Level-by-Level Character Calculator*)


(* ::Text:: *)
(*In this calculation the character comprises of 3 separate fugacities z, y and v. These correspond to gl(1) (osp(2|2)), R-filtration and gl(1) (outer automorphism) respectively. The fugacity corresponding to the sl(2) of osp(2|2) is neglected here since we are at a specific level in its expansion.  Note that due to this extra refinement this method of calculation of the character takes longer than the subsequent section.*)


(* ::Input:: *)
(*RedRemover[g_,i_,j_]:=Module[*)
(*{y=g,char,dta,dta1,sts,sts1,sts2,sm,css,smple,sml,sml1,sml2,sml3,sml4},*)
(*dta=DeleteCases[y[[i,j]],0]//ExpandAll;*)
(*dta1=Fsort[#,z]&/@dta;*)
(*sts=Flatten[Table[If[Head[dta1[[\[Mu]]]]===Plus,List@@dta1[[\[Mu]]],dta1[[\[Mu]]]],{\[Mu],Length[dta1]}]];*)
(*sts1=DeleteCases[Flatten[Table[If[Head[sts[[iii]]]===Times,List@@sts[[iii]],sts[[iii]]],{iii,Length[sts]}]],_Integer];*)
(*sts2=DeleteDuplicates[sts1];*)
(*sts2*)
(*]*)


(* ::Input:: *)
(*rcounter[g_]:=Module[*)
(*{y=g,r1,r2,r3},*)
(*r1=y/.Nl->Times;*)
(*r2=r1/.{ww->dww,jj0->djj0,wwb->dwwb,tt->dtt,gg->dgg,ggt->dggt,ggw->dggw,ggwb->dggwb};(*CNH*)*)
(*r3=rfiltotal[#]&/@r2;*)
(*r3*)
(*]*)


(* ::Input:: *)
(*rlistcreator[n_,list_]:=Module[*)
(*{y=list,strt={},i},*)
(*For[i=Min[y],i<n+1,i++,strt=Join[strt,Position[y,i]]];*)
(*strt*)
(*]*)


(* ::Input:: *)
(*RSingleCharCalculator[g_,i_,n_]:=Module[*)
(*{y=g,w,wrk,char,smb,sm,css,smple,sml,sml1,sml2,sml3,sml4},*)
(*wrk=Extract[y,rlistcreator[n,rcounter[y]]];*)
(*smb=wrk/.Nl->NO;*)
(*sm=Sum[\[Alpha][jj]smb[[jj]],{jj,Length[wrk]}]/.{ww->wwww,jj0->j0j0,wwb->wbwb,tt->tttt,gg->gggg,ggt->gtgt,ggw->gwgw,ggwb->gwbgwb};(*CNH*)*)
(*smple=GetCoefficients[sm];*)
(*char=MatrixRank[Table[Table[Coefficient[smple[[ll]],\[Alpha][kk]],{kk,Length[wrk]}],{ll,Length[smple]}]];*)
(*char*)
(*]*)


(* ::Input:: *)
(*redsubtractor[list_,flist_,ii_,jj_]:=Module[*)
(*{y=list,i,j,res={Style[list[[1]],Bold][Append[If[fcal===Checkbox[False],fcollectioncoloured,sfcollectioncoloured][[ii,jj]],NumericalSort[DeleteDuplicates[rcounter[flist[[ii,jj]]]]][[1]]]]}},*)
(*For[i=1,i<Length[y],i++,res=Append[res,Style[y[[i+1]]-y[[i]],Bold][Append[If[fcal===Checkbox[False],fcollectioncoloured,sfcollectioncoloured][[ii,jj]],NumericalSort[DeleteDuplicates[rcounter[flist[[ii,jj]]]]][[i+1]]]]]];*)
(*res*)
(*]*)


(* ::Input:: *)
(*RFullCharCalculator[x_]:=Module[*)
(*{y=x,chart,chartable},*)
(*chart=Table[Table[Table[If[MemberQ[rcounter[y[[i,j]]],k],RSingleCharCalculator[RedRemover[y,i,j],i,k],Nothing],{k,Min[rcounter[y[[i,j]]]],Max[rcounter[y[[i,j]]]],1/2}],{j,Length[y[[i]]]}],{i,Length[y]}];*)
(*chartable=Table[Table[redsubtractor[chart[[i,j]],y,i,j],{j,Length[chart[[i]]]}],{i,Length[chart]}];*)
(*chartable*)
(*]*)


(* ::Input:: *)
(*rcharacter=Quiet[RFullCharCalculator[If[fcal===Checkbox[False],fullgrouping,sfullgrouping]]/.{MatrixRank[{}]->0}]*)


(* ::Text:: *)
(*The bold number is the degeneracy; inside the square brackets are the quantum numbers: red is for m (gl(1) in osp(2|2)), while blue is for r (outer automorphism) and black is for R (R-filtration).*)


(* ::Input:: *)
(*rfcharacter=If[fcal===Checkbox[False],Quiet[Sum[fermlist[[i,j]]* Head[rcharacter[[i,j,k]]][[1]] z^rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^rcharacter[[i,j,k,1,2,1]],{i,Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]],Quiet[Sum[sfermlist[[i,j]]*Head[rcharacter[[i,j,k]]][[1]] z^rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^rcharacter[[i,j,k,1,2,1]],{i,Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]+Sum[sfermlist[[i,j]]*Head[rcharacter[[i,j,k]]][[1]] z^-rcharacter[[i,j,k,1,1,1]] y^rcharacter[[i,j,k,1,3]] v^-rcharacter[[i,j,k,1,2,1]],{i,If[(IntegerQ[lev]||MemberQ[mparam,0]),2,1],Length[rcharacter]},{j,Length[rcharacter[[i]]]},{k,Length[rcharacter[[i,j]]]}]]]*)


(* ::Text:: *)
(*This is using a specific grading for the character: fugacity z is for the gl(1) (in osp(2|2)), y is for the R (from the R-filtration) and v is  for the outer automorphism of the VOA.*)


(* ::Subsection:: *)
(*Level-by-Level Character Calculator*)


(* ::Text:: *)
(*In this calculation the character comprises of 2 separate fugacities z and v. These correspond to gl(1) (osp(2|2)) and gl(1) (outer automorphism) respectively. The fugacity corresponding to the sl(2) of osp(2|2) is neglected here since we are at a specific level in its expansion. *)


(* ::Input:: *)
(*SingleCharCalculator[g_,i_,j_]:=Module[*)
(*{y=g,char,dta,dta1,sts,sts1,sts2,smb,sm,css,smple,sml,sml1,sml2,sml3,sml4,bsm,aftope},*)
(*dta=DeleteCases[y[[i,j]],0]//ExpandAll;*)
(*dta1=Fsort[#,z]&/@dta;*)
(*sts=Flatten[Table[If[Head[dta1[[\[Mu]]]]===Plus,List@@dta1[[\[Mu]]],dta1[[\[Mu]]]],{\[Mu],Length[dta1]}]];*)
(*sts1=DeleteCases[Flatten[Table[If[Head[sts[[iii]]]===Times,List@@sts[[iii]],sts[[iii]]],{iii,Length[sts]}]],_Integer];*)
(*sts2=DeleteDuplicates[sts1];*)
(*smb=sts2/.Nl->NO;*)
(*sm=Sum[\[Alpha][jj]smb[[jj]],{jj,Length[sts2]}]/.{ww->wwww,jj0->j0j0,wwb->wbwb,tt->tttt,gg->gggg,ggt->gtgt,ggw->gwgw,ggwb->gwbgwb};(*CNH*)*)
(*smple=GetCoefficients[sm];*)
(*char=MatrixRank[Table[Table[Coefficient[smple[[ll]],\[Alpha][kk]],{kk,Length[sts1]}],{ll,Length[smple]}]];*)
(*char*)
(*]*)


(* ::Input:: *)
(*FullCharCalculator[x_]:=Module[*)
(*{y=x,chartable},*)
(*chartable=Table[Table[Style[SingleCharCalculator[y,i,j],Bold][If[fcal===Checkbox[False],fcollectioncoloured,sfcollectioncoloured][[i,j]]],{j,Length[y[[i]]]}],{i,Length[y]}];*)
(*chartable*)
(*]*)


(* ::Input:: *)
(*character=Quiet[FullCharCalculator[If[fcal===Checkbox[False],fullgrouping,sfullgrouping]]/. MatrixRank[{}]->0]*)


(* ::Text:: *)
(*The bold number is the degeneracy; inside the square brackets are the quantum numbers: red is for m (gl(1) in osp(2|2)), while blue is for r (gl(1) outer automorphism).*)


(* ::Input:: *)
(*fcharacter=If[fcal===Checkbox[False],Quiet[Sum[fermlist[[i,j]]*Head[character[[i,j]]][[1]] z^character[[i,j,1,1,1]] v^character[[i,j,1,2,1]],{i,Length[character]},{j,Length[character[[i]]]}]/. MatrixRank[{}]->0],Quiet[Sum[sfermlist[[i,j]]*Head[character[[i,j]]][[1]] z^character[[i,j,1,1,1]] v^character[[i,j,1,2,1]],{i,Length[character]},{j,Length[character[[i]]]}]+Sum[sfermlist[[i,j]]*Head[character[[i,j]]][[1]] z^-character[[i,j,1,1,1]] v^-character[[i,j,1,2,1]],{i,If[(IntegerQ[lev]||MemberQ[mparam,0]),2,1],Length[character]},{j,Length[character[[i]]]}]/. MatrixRank[{}]->0]]*)


(* ::Text:: *)
(*This is using a specific grading for the character: fugacity z is for the gl(1) (in osp(2|2)), v is for the outer automorphism of the VOA.*)
