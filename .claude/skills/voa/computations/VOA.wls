#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
Import["freeFieldVOA.m"]
SetOptions[EvaluationNotebook[],CommonDefaultFormatTypes->{"Output"->StandardForm}]


(* ::Title::Closed:: *)
(*Virasoro*)


(* ::Input:: *)
(*(* \:6ce8\:518c\:73bb\:8272\:7b97\:7b26 T *)*)
(*Bosonic[T];*)


(* ::Input:: *)
(*(* \:786e\:5b9a OPE *)*)
(*OPE[T,T]=MakeOPE[((c/2)One[w])/(z-w)^4+(2T[w])/(z-w)^2+Derivative[1][T][w]/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*(* \:8ba1\:7b97 NO \:548c\:590d\:5408\:7b97\:7b26 OPE *)*)
(*(* OPEdefs \:7684\:5d4c\:5957\:6b63\:89c4\:4e58\:79ef NO \:7684 \[OpenCurlyDoubleQuote]\:6807\:51c6\:6392\:5e8f\[CloseCurlyDoubleQuote] \:662f NO[a, NO[b, NO[c, NO[d,, ...]]]]\:ff0c\:5373 NO \:9760\:53f3\:653e\:7f6e *)*)
(*(* NO[NO[a, b], c] \:8981\:5229\:7528 OPE \:7684\:6781\:70b9\:4fe1\:606f\:83b7\:5f97 *)*)
(*{NO[Derivative[1][T],Derivative[1][T]],NO[NO[Derivative[1][Derivative[1][Derivative[1][T]]],Derivative[1][Derivative[1][Derivative[1][T]]]],Derivative[1][T]],NO[T,NO[NO[Derivative[1][Derivative[1][Derivative[1][T]]],Derivative[1][Derivative[1][Derivative[1][T]]]],Derivative[1][T]]]}//Expand*)
(*{OPE[Derivative[1][T],Derivative[1][T]]}*)
(*(* \:5229\:7528 OPE \:8ba1\:7b97 conformal weight = 4 \:7684 null states \:7684\:5b58\:5728\:6761\:4ef6 *)*)
(*{OPEPole[6][T,NO[T,T]+ppp Derivative[1][Derivative[1][T]]],OPEPole[4][T,NO[T,T]+ppp Derivative[1][Derivative[1][T]]]}==0//Solve*)


(* ::Title::Closed:: *)
(*bc\[Beta]\[Gamma]*)


(* ::Input:: *)
(*(* \:6ce8\:518c\:57fa\:672c\:7b97\:7b26 *)*)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)


(* ::Input:: *)
(*OPE[b,c]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]];*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[-(One[w]/(z-w))+Ord[z,w,0]];*)


(* ::Input:: *)
(*{NO[\[Beta],\[Gamma]],NO[\[Gamma],\[Beta]],NO[\[Beta],\[Gamma]]-NO[\[Gamma],\[Beta]],NO[\[Beta],\[Beta]],NO[\[Gamma],\[Gamma]]}*)
(*{NO[NO[Derivative[1][\[Beta]],\[Gamma]],\[Gamma]],NO[\[Gamma],NO[Derivative[1][\[Beta]],\[Gamma]]]}*)
(*{NO[NO[Derivative[1][\[Beta]],Derivative[1][\[Gamma]]],NO[Derivative[1][\[Beta]],Derivative[1][\[Gamma]]]]}//Expand*)
(*OPE[NO[Derivative[1][\[Beta]],Derivative[1][\[Gamma]]],NO[Derivative[1][\[Beta]],Derivative[1][\[Gamma]]]]*)


(* ::Input:: *)
(*T=(\[Lambda]-1)NO[\[Beta],Derivative[1][\[Gamma]]]+\[Lambda] NO[\[Gamma],Derivative[1][\[Beta]]];*)
(*{NO[T,\[Beta]]//Expand,OPE[T,\[Beta]]}//Simplify*)
(*OPEPole[2][T,\[Beta]]-(1-\[Lambda])\[Beta]//Simplify*)
(*OPEPole[1][T,\[Beta]]-Derivative[1][\[Beta]]//Simplify*)


(* ::Input:: *)
(*{NO[b,b],NO[c,c],NO[Derivative[1][b],Derivative[1][b]],NO[Derivative[1][b],Derivative[1][c]],NO[Derivative[1][b],Derivative[1][Derivative[1][b]]],NO[Derivative[1][c],Derivative[1][Derivative[1][c]]]}*)
(*{NO[NO[b,c],NO[b,c]],NO[NO[Derivative[1][b],c],NO[b,Derivative[1][c]]]}//Simplify*)
(*{Derivative[1][NO[b,c]],Derivative[1][NO[Derivative[1][b],c]],NO[NO[Derivative[1][b],Derivative[1][c]],NO[b,c]]}//Simplify*)
(*{NO[NO[Derivative[1][Derivative[1][b]],Derivative[1][c]],NO[b,Derivative[1][c]]]}//Simplify*)
(**)
(*T=(\[Lambda]-1)NO[b,Derivative[1][c]]+\[Lambda] NO[Derivative[1][b],c];*)
(*{NO[T,b]//Expand,OPE[T,b]}//Simplify*)
(*OPEPole[2][T,b]-(1-\[Lambda])b//Simplify*)
(*OPEPole[1][T,b]-Derivative[1][b]//Simplify*)


(* ::Input:: *)
(**)


(* ::Title::Closed:: *)
(*Subscript[U(1), k]*)


(* ::Input:: *)
(*(* \:6ce8\:518c T, J *)*)
(*Bosonic[T,J]*)


(* ::Input:: *)
(*(* \:672a\:6ce8\:518c OPE \:65f6\:7684 NO *)*)
(*{OPE[T,J],NO[T,J],NO[J,NO[T,J]],NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][T]]],NO[Derivative[1][T],NO[Derivative[1][T],J]],NO[NO[Derivative[1][T],NO[Derivative[1][T],J]],NO[Derivative[1][T],NO[Derivative[1][T],J]]]}*)


(* ::Input:: *)
(*(* \:6ce8\:518c OPE *)*)
(*OPE[J,J]:=MakeOPE[One[w]/(z-w)^2+Ord[z,w,0]];*)
(*OPE[T,J]:=MakeOPE[J[w]/(z-w)^2+Derivative[1][J][w]/(z-w)+Ord[z,w,0]];*)
(*OPE[T,T]:=MakeOPE[(c/2 One[w])/(z-w)^4+(2T[w])/(z-w)^2+Derivative[1][T][w]/(z-w)+Ord[z,w,0]];*)


(* ::Input:: *)
(*(* \:6ce8\:518c OPE \:4e4b\:540e\:7684 NO *)*)
(**)
(*{OPE[T,J],NO[T,J],NO[J,NO[T,J]],NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][T]]],NO[Derivative[1][T],NO[Derivative[1][T],J]],*)
(*NO[NO[T,J],NO[T,NO[Derivative[1][T],J]]],*)
(*NO[NO[Derivative[1][T],NO[T,J]],NO[T,NO[Derivative[1][T],J]]]*)
(*}//Simplify*)


(* ::Input:: *)
(*OPE[T,Derivative[1][J]]//OPEToSeries*)
(*OPE[T,T]//OPEToSeries*)
(*OPE[Derivative[1][T],T]//OPEToSeries*)
(*OPE[NO[J,J],J]//OPEToSeries*)


(* ::Title:: *)
(*Subscript[sl(2), k]*)


(* ::Chapter::Closed:: *)
(*J^(1,2,3) basis*)


(* ::Input:: *)
(*Bosonic[Jp[_]]*)
(*OPE[Jp[a_],Jp[b_]]:=MakeOPE[(k 1/2 Delta[a,b]One[w])/(z-w)^2+I Sum[LeviCivitaTensor[3][[a,b,c]]*1/(z-w) Jp[c][w],{c,3}]+Ord[z,w,0]]*)


(* ::Input:: *)
(*Table[OPEPole[q][Jp[1],Derivative[1][Jp[2]]],{q,0,3}]*)
(*Table[(q-1)OPEPole[q-1][Jp[1],Jp[2]]+Derivative[1][OPEPole[q][Jp[1],Jp[2]]],{q,0,3}]*)
(**)
(*Table[OPEPole[q][Jp[1],Derivative[1][Jp[1]]],{q,0,3}]*)
(*Table[(q-1)OPEPole[q-1][Jp[1],Jp[1]]+Derivative[1][OPEPole[q][Jp[1],Jp[1]]],{q,0,3}]*)


(* ::Input:: *)
(*NO[Jp[1],NO[Jp[2],Jp[3]]]-NO[Jp[1],NO[Jp[3],Jp[2]]]-NO[Jp[2],NO[Jp[1],Jp[3]]]+NO[Jp[2],NO[Jp[3],Jp[1]]]+NO[Jp[3],NO[Jp[1],Jp[2]]]-NO[Jp[3],NO[Jp[2],Jp[1]]]//Simplify*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*TT=1/(2(k+2)) Sum[2KroneckerDelta[a,b]NO[Jp[a],Jp[b]],{a,3},{b,3}]*)
(*OPE[TT,Jp[1]]//OPEToSeries//Simplify*)
(*OPE[TT,Jp[1]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*Sum[(NNO[J[bb],NNO[J[cc],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2)+(NNO[J[cc],NNO[J[bb],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2)-(NNO[J[cc],NNO[J[dd],J[bb]]] f[a][aa,bb] f[aa][cc,dd])/(hcheck+k)^2-(NNO[J[aa],NNO[J[bb],J[cc]]] f[][aa,bb,dd] f[a,dd][cc])/(4 (hcheck+k)^2),{aa,3},{bb,3},{cc,3},{dd,3}]/.{a->3,J->Jp}/.{f[a_][b_,c_]:>2LeviCivitaTensor[3][[a,b,c]],f[][a_,b_,c_]:>4LeviCivitaTensor[3][[a,b,c]],f[a_,b_][c_]:>LeviCivitaTensor[3][[a,b,c]]}/.NNO->NO//Simplify*)


(* ::Input:: *)
(*\[Chi]=NO[TT,TT]+\[Alpha]/2 Derivative[1][Derivative[1][TT]]+Sum[\[Beta]/2 KKD[ap,bp]NO[Derivative[1][Derivative[1][Jp[ap]]],Jp[bp]],{ap,3},{bp,3}]+\[Gamma] Sum[KKD[ap,bp]NO[Derivative[1][Jp[ap]],Derivative[1][Jp[bp]]],{ap,3},{bp,3}];*)
(**)
(*(OPE[Jp[3],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-4}]&//Simplify)/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}//Simplify*)
(**)
(*(OPE[Jp[3],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&//Simplify)/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}//Simplify*)
(**)
(*(OPE[Jp[3],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-2}]&//Simplify)/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}//Simplify*)


(* ::Input:: *)
(*Sum[(k KKD[aa,bb] NNO[J[a],NNO[J[aa],J[bb]]])/(4 (hcheck+k)^2)+(k KKD[aa,bb] NNO[J[aa],NNO[J[a],J[bb]]])/(2 (hcheck+k)^2)+((4 hcheck+k) KKD[aa,bb] NNO[J[aa],NNO[J[bb],J[a]]])/(4 (hcheck+k)^2),{aa,3},{bb,3}]+(1/(2 (hcheck+k)^2))(k^2 (1+\[Alpha]+k \[Beta])+2 hcheck^3 (\[Beta]+\[Gamma])+hcheck k (3+4 \[Alpha]+4 k \[Beta]+2 k \[Gamma])+hcheck^2 (3+3 \[Alpha]+5 k \[Beta]+4 k \[Gamma])) J[a]^\[Prime]\[Prime]+Sum[1/(4 (hcheck+k)^2) I ((4 hcheck^2 \[Beta]+4 hcheck (3+\[Alpha]+2 k \[Beta])+k (3+4 \[Alpha]+4 k \[Beta])) NNO[J[aa],Derivative[1][J[bb]]]+(8 hcheck^2 \[Gamma]+k (-1+4 \[Alpha]+8 k \[Gamma])+hcheck (3+4 \[Alpha]+16 k \[Gamma])) NNO[Derivative[1][J[aa]],J[bb]]) f[a][aa,bb],{aa,3},{bb,3}]+Sum[(NNO[J[bb],NNO[J[cc],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2)+(NNO[J[cc],NNO[J[bb],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2)-(NNO[J[cc],NNO[J[dd],J[bb]]] f[a][aa,bb] f[aa][cc,dd])/(hcheck+k)^2-(NNO[J[aa],NNO[J[bb],J[dd]]] f[][aa,bb,cc] f[a,cc][dd])/(4 (hcheck+k)^2),{aa,3},{bb,3},{cc,3},{dd,3}]/.{a->3,J->Jp}/.{NNO->NO}/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}/.{f[a_][b_,c_]:>2LeviCivitaTensor[3][[a,b,c]],f[][a_,b_,c_]:>4LeviCivitaTensor[3][[a,b,c]],f[a_,b_][c_]:>LeviCivitaTensor[3][[a,b,c]]}//Simplify*)
(*%/.{hcheck->2}//Simplify*)


(* ::Input:: *)
(*(Sum[+((NNO[J[bb],NNO[J[cc],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2))+(NNO[J[cc],NNO[J[bb],J[dd]]] f[a][aa,bb] f[aa][cc,dd])/(2 (hcheck+k)^2)-(NNO[J[cc],NNO[J[dd],J[bb]]] f[a][aa,bb] f[aa][cc,dd])/(hcheck+k)^2-(NNO[J[aa],NNO[J[bb],J[dd]]] f[][aa,bb,cc] f[a,cc][dd])/(4 (hcheck+k)^2),{aa,3},{bb,3},{cc,3},{dd,3}])/.{a->3,J->Jp}/.{NNO->NO}/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}/.{f[a_][b_,c_]:>2LeviCivitaTensor[3][[a,b,c]],f[][a_,b_,c_]:>4LeviCivitaTensor[3][[a,b,c]],f[a_,b_][c_]:>LeviCivitaTensor[3][[a,b,c]]}//Simplify*)
(**)
(**)
(*(Sum[NNO[J[bb],NNO[J[cc],J[dd]]]/(hcheck+k)^2 (+(( f[a][aa,bb] f[aa][cc,dd])/2 )+( f[a][aa,cc] f[aa][bb,dd])/2 - f[a][aa,dd] f[aa][bb,cc]-( f[a][aa,dd]f[aa][bb,cc] )/4),{aa,3},{bb,3},{cc,3},{dd,3}])/.{a->3,J->Jp}/.{NNO->NO}/.{KKD[a_,b_]:>2KroneckerDelta[a,b]}/.{f[a_][b_,c_]:>2LeviCivitaTensor[3][[a,b,c]],f[][a_,b_,c_]:>4LeviCivitaTensor[3][[a,b,c]],f[a_,b_][c_]:>LeviCivitaTensor[3][[a,b,c]]}//Simplify*)


(* ::Input:: *)
(*+(( f[a][aa,bb] f[aa][cc,dd])/2 )+( f[a][aa,cc] f[aa][bb,dd])/2 - f[a][aa,dd] f[aa][bb,cc]-( f[a][aa,dd]f[aa][bb,cc] )/4*)


(* ::Input:: *)
(*Table[Sum[-f[c][b,d]f[d][a,e]+f[c][a,d]f[d][b,e]+f[c][e,d]f[d][a,b],{e,3}],{a,3},{b,3},{c,3},{d,3}]/.{f[a_][b_,c_]:>2LeviCivitaTensor[3][[a,b,c]]}*)


(* ::Chapter::Closed:: *)
(*J3, SuperPlus[J], SuperMinus[J]*)


(* ::Input:: *)
(**)
(*Bosonic[J3,Jplus,Jminus];*)
(**)
(*OPE[J3,J3]=MakeOPE[((k/2) One[w])/(z-w)^2+Ord[z,w,0]]*)
(*OPE[J3,Jplus]=MakeOPE[Jplus[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[J3,Jminus]=MakeOPE[-Jminus[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[Jplus,Jminus]=MakeOPE[(k One[w])/(z-w)^2+(2J3[w])/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*(2NO[Derivative[1][J3],J3]+NO[Derivative[1][Jplus],Jminus]+NO[Derivative[1][Jminus],Jplus])-(2NO[J3,Derivative[1][J3]]+NO[Jplus,Derivative[1][Jminus]]+NO[Jminus,Derivative[1][Jplus]])*)
(*(2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])-(2NO[J3,Derivative[1][Derivative[1][J3]]]+NO[Jplus,Derivative[1][Derivative[1][Jminus]]]+NO[Jminus,Derivative[1][Derivative[1][Jplus]]])*)


(* ::Input:: *)
(*T=1/(2(k+2)) (2NO[J3,J3]+NO[Jplus,Jminus]+NO[Jminus,Jplus])*)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-4}]&*)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-3}]&*)
(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-2}]&)-(2T//MakeVertexField[#][w]&)//Simplify*)
(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-1}]&)-(T'//MakeVertexField[#][w]&)//Simplify*)


(* ::Input:: *)
(*OPE[T,J3]//OPEToSeries//Simplify*)


(* ::Subsection:: *)
(*level-4 null descendant*)


(* ::Input:: *)
(*Derivative[1][Derivative[1][T]]-(1/(2+k) (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+1/(2+k) (2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]]))//Simplify*)


(* ::Input:: *)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])(*+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])*)//Simplify*)
(**)
(*OPE[J3,\[Chi]]//OPEToSeries//Simplify*)
(**)


(* ::Input:: *)
(*T\[Chi]OPE=OPE[T,\[Chi]]//OPEToSeries//Simplify;*)
(*T\[Chi]OPE//SeriesCoefficient[#,{z,w,-3}]&*)
(*T\[Chi]OPE//SeriesCoefficient[#,{z,w,-4}]&*)
(*T\[Chi]OPE//SeriesCoefficient[#,{z,w,-5}]&*)


(* ::Input:: *)
(*\[Beta]\[Gamma]sol=Solve[{3+3 \[Alpha]+(4+3 k) \[Beta]+4 \[Gamma]==0,\[Alpha]+\[Beta]+2 (1+k) \[Gamma]==0},{\[Beta],\[Gamma]}][[1]]*)
(*\[Beta]-2\[Gamma]/.\[Beta]\[Gamma]sol//Simplify*)


(* ::Input:: *)
(*\[Chi]p=(OPE[J3,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*NO[Jplus,T]+1/6 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Expand*)
(*OPE[J3,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)
(*OPE[Jplus,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)
(*OPE[Jminus,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*\[Chi]p=(OPE[J3,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*(* this is actually the Subscript[v, sing] in Adamovic https://arxiv.org/abs/math/0401023 *)*)
(*\[Chi]p=(OPE[Jplus,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*\[Chi]p=(OPE[Jminus,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*OPE[J3,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(*OPE[Jplus,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(*OPE[Jminus,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(**)
(**)
(*OPE[J3,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(*OPE[Jplus,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(*OPE[Jminus,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(**)
(**)
(*OPE[J3,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)
(*OPE[Jplus,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)
(*OPE[Jminus,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*Solve[{-2+k+k^2==0}]*)


(* ::Input:: *)
(*(* k = -4/3 *)*)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])(*+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])*)//Simplify*)
(**)
(*OPE[J3,\[Chi]]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*\[Chi]p=(%//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}*)
(**)
(*OPE[J3,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(**)


(* ::Subsection::Closed:: *)
(*level-4 zero vector*)


(* ::Input:: *)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])/.{\[Beta]->-(3/2)-(3 \[Alpha])/2,\[Gamma]->-((3 \[Alpha])/4),\[Delta]->-(3/4)}//Simplify*)


(* ::Title::Closed:: *)
(*Subscript[sl(3), k]*)


(* ::Input:: *)
(*Bosonic[Jp[_]]*)
(*OPE[Jp[a_],Jp[b_]]:=MakeOPE[(k 1/2 Delta[a,b]One[w])/(z-w)^2+I Sum[LeviCivitaTensor[3][[a,b,c]]*1/(z-w) Jp[c][w],{c,3}]+Ord[z,w,0]]*)


(* ::Input:: *)
(*NO[Jp[1],NO[Jp[2],Jp[3]]]-NO[Jp[1],NO[Jp[3],Jp[2]]]-NO[Jp[2],NO[Jp[1],Jp[3]]]+NO[Jp[2],NO[Jp[3],Jp[1]]]+NO[Jp[3],NO[Jp[1],Jp[2]]]-NO[Jp[3],NO[Jp[2],Jp[1]]]//Simplify*)


(* ::Input:: *)
(*TT=1/(2(k+2)) Sum[2KroneckerDelta[a,b]NO[Jp[a],Jp[b]],{a,3},{b,3}]*)
(*OPE[Jp[1],NO[TT,TT]]//OPEToSeries//SeriesCoefficient[#,{z,w,-4}]&*)


(* ::Input:: *)
(**)
(*Bosonic[J3,Jplus,Jminus];*)
(**)
(*OPE[J3,J3]=MakeOPE[((k/2) One[w])/(z-w)^2+Ord[z,w,0]]*)
(*OPE[J3,Jplus]=MakeOPE[Jplus[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[J3,Jminus]=MakeOPE[-Jminus[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[Jplus,Jminus]=MakeOPE[(k One[w])/(z-w)^2+(2J3[w])/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*T=1/(2(k+2)) (2NO[J3,J3]+NO[Jplus,Jminus]+NO[Jminus,Jplus])*)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-4}]&*)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-3}]&*)
(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-2}]&)-(2T//MakeVertexField[#][w]&)//Simplify*)
(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-1}]&)-(T'//MakeVertexField[#][w]&)//Simplify*)


(* ::Subsection:: *)
(*level-4 null descendant*)


(* ::Input:: *)
(*(\[Chi]/.{\[Beta]->-((3+3 k+\[Alpha]+3 k \[Alpha])/(2+7 k+3 k^2)),\[Gamma]->-((-3+\[Alpha]+3 k \[Alpha])/(2 (2+7 k+3 k^2)))}//Simplify)//SeriesCoefficient[#,{\[Alpha],0,1}]&*)


(* ::Input:: *)
(*Derivative[1][Derivative[1][T]]-(1/(2+k) (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+1/(2+k) (2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]]))//Simplify*)


(* ::Input:: *)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])(*+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])*)//Simplify*)
(**)
(*OPE[J3,\[Chi]]//OPEToSeries//Simplify*)
(**)


(* ::Input:: *)
(*\[Beta]\[Gamma]sol=Solve[{3+3 \[Alpha]+(4+3 k) \[Beta]+4 \[Gamma]==0,\[Alpha]+\[Beta]+2 (1+k) \[Gamma]==0},{\[Beta],\[Gamma]}][[1]]*)
(*\[Beta]-2\[Gamma]/.\[Beta]\[Gamma]sol//Simplify*)


(* ::Input:: *)
(*\[Chi]p=(OPE[J3,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*NO[Jplus,T]+1/6 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Expand*)
(*OPE[J3,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)
(*OPE[Jplus,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)
(*OPE[Jminus,NO[Jplus,T]+(1/2)/3 Derivative[1][Derivative[1][Jplus]]-NO[J3,Derivative[1][Jplus]]+NO[Derivative[1][J3],Jplus]//Simplify]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*\[Chi]p=(OPE[J3,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*(* this is actually the Subscript[v, sing] in Adamovic https://arxiv.org/abs/math/0401023 *)*)
(*\[Chi]p=(OPE[Jplus,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*k=-4/3;*)
(*\[Chi]p=(OPE[Jminus,\[Chi]]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.\[Beta]\[Gamma]sol//Simplify*)
(*OPE[J3,\[Chi]p]/.\[Beta]\[Gamma]sol/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.\[Beta]\[Gamma]sol/.{}/.\!\(TraditionalForm\`\(TraditionalForm\`{}\)\)//OPEToSeries//Simplify*)


(* ::Input:: *)
(*OPE[J3,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(*OPE[Jplus,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(*OPE[Jminus,(1+3 k) NO[J3,J3]+(-3+k) (NO[Jplus,Jminus]-Derivative[1][J3])]//OPEToSeries//Simplify*)
(**)
(**)
(*OPE[J3,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(*OPE[Jplus,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(*OPE[Jminus,2 NO[J3,Jplus][w]-Derivative[1][Jplus][w]]//OPEToSeries//Simplify*)
(**)
(**)
(*OPE[J3,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)
(*OPE[Jplus,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)
(*OPE[Jminus,2 NO[J3,Jminus][w]+Derivative[1][Jminus][w]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*Solve[{-2+k+k^2==0}]*)


(* ::Input:: *)
(*(* k = -4/3 *)*)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])(*+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])*)//Simplify*)
(**)
(*OPE[J3,\[Chi]]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*\[Chi]p=(%//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}*)
(**)
(*OPE[J3,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*OPE[Jplus,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(*OPE[Jminus,\[Chi]p]/.{k->-(4/3)}/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> \(-*)
(*\*FractionBox[\(1\), \(2\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(2\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(3\), \(4\)]\) - *)
(*\*FractionBox[\(3\ \[Alpha]\), \(4\)]}\)\)//OPEToSeries//Simplify*)
(**)


(* ::Subsection::Closed:: *)
(*level-4 zero vector*)


(* ::Input:: *)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 (2NO[Derivative[1][Derivative[1][J3]],J3]+NO[Derivative[1][Derivative[1][Jplus]],Jminus]+NO[Derivative[1][Derivative[1][Jminus]],Jplus])+\[Gamma](2NO[Derivative[1][J3],Derivative[1][J3]]+NO[Derivative[1][Jplus],Derivative[1][Jminus]]+NO[Derivative[1][Jminus],Derivative[1][Jplus]])+\[Delta](2NO[J3,NO[J3,T]]+NO[Jplus,NO[Jminus,T]]+NO[Jminus,NO[Jplus,T]])/.{\[Beta]->-(3/2)-(3 \[Alpha])/2,\[Gamma]->-((3 \[Alpha])/4),\[Delta]->-(3/4)}//Simplify*)


(* ::Title:: *)
(*D4*)


(* ::Chapter::Closed:: *)
(*data for Subscript[D, 4] Lie algebra*)


(* ::Input:: *)
(*hD4={DiagonalMatrix[{1/2,1/2,0,0,-1/2,-1/2,0,0}],DiagonalMatrix[{1/2,-1/2,0,0,-1/2,1/2,0,0}],DiagonalMatrix[{0,0,1/2,-1/2,0,0,-1/2,1/2}],DiagonalMatrix[{0,0,1/2,1/2,0,0,-1/2,-1/2}]}*)
(*ro1=Table[ArrayFlatten[{{Table[KroneckerDelta[k,i]KroneckerDelta[l,j],{k,1,4},{l,1,4}],0},{0,Table[-KroneckerDelta[k,j]KroneckerDelta[l,i],{k,1,4},{l,1,4}]}}],{i,1,3},{j,i+1,4}]//Flatten[#,1]&*)
(*ro2=Table[ArrayFlatten[{{0,Table[KroneckerDelta[k,i]KroneckerDelta[l,j]-KroneckerDelta[k,j]KroneckerDelta[l,i],{k,1,4},{l,1,4}]},{Table[0,{m,1,4},{n,1,4}],0}}],{i,1,3},{j,i+1,4}]//Flatten[#,1]&*)
(*lo1=Table[Transpose[ro1[[i]]],{i,1,6}]*)
(*lo2=Table[Transpose[ro2[[i]]],{i,1,6}]*)
(*g=Join[hD4,ro1,ro2,lo1,lo2]*)
(*gC=Join[hD4,1/2 lo1,1/2 lo2,1/2 ro1,1/2 ro2]*)


(* ::Input:: *)
(*transfer=ArrayFlatten[{{IdentityMatrix[4]/\[Sqrt]2,IdentityMatrix[4]/\[Sqrt]2},{I IdentityMatrix[4]/\[Sqrt]2,-I IdentityMatrix[4]/\[Sqrt]2}}]*)
(*MatrixForm[#]&/@Table[transfer . hD4[[i]] . Inverse[transfer],{i,1,4}]*)
(*MatrixForm[#]&/@Table[transfer . ro1[[i]] . Inverse[transfer],{i,1,6}]*)


(* ::Input:: *)
(**)
(*l=({*)
(* {1/\[Sqrt]2, b},*)
(* {c, d}*)
(*})*)
(*l . ({*)
(* {0, 1},*)
(* {1, 0}*)
(*}) . Transpose[l]-({*)
(* {1, 0},*)
(* {0, 1}*)
(*})//Flatten*)
(*Solve[Table[%[[i]]==0,{i,1,4}],{a,b,c,d}]*)


(* ::Input:: *)
(*Table[Tr[g[[i]] . gC[[j]]],{i,1,28},{j,1,28}]==IdentityMatrix[28]*)
(*K=Table[Tr[g[[i]] . g[[j]]],{i,1,28},{j,1,28}]*)
(*KD=Inverse[K]*)


(* ::Input:: *)
(*Com[x_,y_]:=x . y-y . x*)
(*fijUkD=Table[Tr[Com[g[[i]],g[[j]]] . gC[[k]]],{i,1,28},{j,1,28},{k,1,28}]*)


(* ::Input:: *)
(**)
(*Table[fijUkD[[1,j,k]],{j,5,28},{k,5,28}]//MatrixForm*)
(*\[Alpha][jj_][1]:=Table[fijUkD[[1,j,k]],{j,5,28},{k,5,28}][[jj,jj]]*)
(*Sum[\[Alpha][j][1]\[Alpha][j][1],{j,5-4,28-4}]*)


(* ::Input:: *)
(*weight=Table[Tr[Com[hD4[[i]],g[[j]]] . gC[[j]]],{j,1,28},{i,1,4}]*)


(* ::Input:: *)
(**)


(* ::Chapter::Closed:: *)
(*Null descendant*)


(* ::Input:: *)
(*Bosonic[J[i_]]*)
(*OPE[J[i_],J[j_]]:=MakeOPE[(k (1/2) K[[i,j]])/(z-w)^2+Sum[fijUkD[[i,j,k]]J[k],{k,1,28}]/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*T=1/(2(k+6)) Sum[1/(1/2) KD[[i,j]]NO[J[i],J[j]],{i,1,28},{j,1,28}]*)
(**)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-4}]&//Simplify*)
(*(OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-3}]&*)
(*(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-2}]&)-(2T//MakeVertexField[#][w]&)//Simplify*)
(*((OPE[T,T]//OPEToSeries//Simplify)//SeriesCoefficient[#,{z,w,-1}]&)-(T'//MakeVertexField[#][w]&)//Simplify*)*)


(* ::Input:: *)
(*(k*28)/(2(k+6))*)


(* ::Subsection:: *)
(*level-4 zero vector*)


(* ::Input:: *)
(*\[Chi]=(NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]*)
(*+\[Beta]/2 Sum[1/(1/2) KD[[i,j]]NO[Derivative[1][Derivative[1][J[i]]],J[j]],{i,28},{j,28}]*)
(*+\[Gamma] Sum[1/(1/2) KD[[i,j]]NO[Derivative[1][J[i]],Derivative[1][J[j]]],{i,28},{j,28}]*)
(*+\[Delta] Sum[1/(1/2) KD[[i,j]]NO[J[i],NO[J[j],T]],{i,28},{j,28}])/.{\[Beta]->-((1+\[Alpha])/(k+hcheck)),\[Gamma]->-(\[Alpha]/(2(k+hcheck))),\[Delta]->-(1/(2(k+hcheck)))}/.{hcheck->6}*)


(* ::Input:: *)
(*%//Simplify*)


(* ::Subsection::Closed:: *)
(*level-4 null descendant*)


(* ::Input:: *)
(*\[Chi]=(NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]*)
(*+\[Beta]/2 Sum[1/(1/2) KD[[i,j]]NO[Derivative[1][Derivative[1][J[i]]],J[j]],{i,28},{j,28}]*)
(*+\[Gamma] Sum[1/(1/2) KD[[i,j]]NO[Derivative[1][J[i]],Derivative[1][J[j]]],{i,28},{j,28}]*)
(*(*+\[Delta] Sum[1/(1/2)KD\[LeftDoubleBracket]i,j\[RightDoubleBracket]NO[J[i],NO[J[j],T]],{i,28},{j,28}]*))//Expand*)


(* ::Input:: *)
(**)
(*\[Chi]/.{\[Beta]->-((1+\[Alpha])/(k+hcheck)),\[Gamma]->-(\[Alpha]/(2(k+hcheck))),\[Delta]->-1/(2(k+hcheck))}/.{hcheck->6}//Simplify*)


(* ::Input:: *)
(*J1\[Chi]OPE=OPE[J[1],\[Chi]]/.{k->-2}//OPEToSeries//Simplify*)


(* ::Input:: *)
(*J1\[Chi]OPE//SeriesCoefficient[#,{z,w,-4}]&*)
(*J1\[Chi]OPE//SeriesCoefficient[#,{z,w,-3}]&*)
(*\[Chi]pField=J1\[Chi]OPE//SeriesCoefficient[#,{z,w,-2}]&*)


(* ::Input:: *)
(*Solve[{\!\(TraditionalForm\`\[Alpha] + 2\ \[Beta] + 4\ \[Gamma] + 1\)==0,\!\(TraditionalForm\`\(TraditionalForm\`\[Alpha] + 3\ \[Beta] + 2\ \[Gamma]\)\)==0},{\[Beta],\[Gamma]}]*)


(* ::Input:: *)
(*\[Chi]p=(\[Chi]pField/.\!\(TraditionalForm\`\(TraditionalForm\`{\[Beta] -> *)
(*\*FractionBox[\(1 - \[Alpha]\), \(4\)], \[Gamma] -> *)
(*\*FractionBox[\(1\), \(8\)]\ \((\(-\[Alpha]\) - 3)\)}\)\)//Simplify)/.{X_[w]:>X}*)
(*J1\[Chi]pOPE=OPE[J[1],\[Chi]p]/.{k->-2}//OPEToSeries//Simplify*)
(*\[Chi]ppField=J1\[Chi]pOPE//SeriesCoefficient[#,{z,w,-2}]&*)
(*\[Chi]pp=\[Chi]ppField/.{X_[w]:>X}*)
(*Table[OPE[J[n],\[Chi]pp]/.{k->-2}//OPEToSeries//Simplify,{n,28}]*)
(*%/.{k->-2}//Simplify*)


(* ::Input:: *)
(*Out[62]*)
(*Out[62]/.{k->-2}/.\!\(TraditionalForm\`{\[Beta] -> *)
(*\*FractionBox[\(1\), \(4\)] - *)
(*\*FractionBox[\(\[Alpha]\), \(4\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(\[Alpha]\), \(8\)]\) - *)
(*\*FractionBox[\(3\), \(8\)]}\)//SeriesCoefficient[#,{z,w,-4}]&//Simplify*)
(*Out[62]/.{k->-2}/.\!\(TraditionalForm\`{\[Beta] -> *)
(*\*FractionBox[\(1\), \(4\)] - *)
(*\*FractionBox[\(\[Alpha]\), \(4\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(\[Alpha]\), \(8\)]\) - *)
(*\*FractionBox[\(3\), \(8\)]}\)//SeriesCoefficient[#,{z,w,-3}]&//Simplify*)
(*Out[62]/.{k->-2}/.\!\(TraditionalForm\`{\[Beta] -> *)
(*\*FractionBox[\(1\), \(4\)] - *)
(*\*FractionBox[\(\[Alpha]\), \(4\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(\[Alpha]\), \(8\)]\) - *)
(*\*FractionBox[\(3\), \(8\)]}\)//SeriesCoefficient[#,{z,w,-2}]&//Simplify*)
(*Out[62]/.{k->-2}/.\!\(TraditionalForm\`{\[Beta] -> *)
(*\*FractionBox[\(1\), \(4\)] - *)
(*\*FractionBox[\(\[Alpha]\), \(4\)], \[Gamma] -> \(-*)
(*\*FractionBox[\(\[Alpha]\), \(8\)]\) - *)
(*\*FractionBox[\(3\), \(8\)]}\)//SeriesCoefficient[#,{z,w,-1}]&//Simplify*)


(* ::Input:: *)
(*Out[70]*)
(*(Out[70]/.{X_[w]:>X}//Expand)*)
(*OPE[J[1],%]/.{k->-2}//OPEToSeries//Simplify*)
(**)
(*(%/.{X_[w]:>X}//Expand)//SeriesCoefficient[#,{z,w,-2}]&*)
(*Table[OPE[J[n],%]/.{k->-2}//OPEToSeries//Simplify,{n,28}]*)


(* ::Input:: *)
(*For[m=1,m<29,m++,*)
(*(Out[70]/.{X_[w]:>X}//Expand);*)
(*OPE[J[m],%]/.{k->-2}//OPEToSeries//Simplify;*)
(*(%/.{X_[w]:>X}//Expand);*)
(*Print[Table[OPE[J[n],%]/.{k->-2}//OPEToSeries//Simplify,{n,28}]]*)
(*]*)


(* ::Title::Closed:: *)
(*Subscript[\[ScriptCapitalW], 3]*)


(* ::Input:: *)
(*(* \:6ce8\:518c\:4e24\:4e2a\:73bb\:8272\:7b97\:7b26 *)*)
(*Bosonic[T,W]*)


(* ::Input:: *)
(*(* \:68c0\:6d4b\:8fd8\:6ca1\:5b9a\:4e49 OPE \:7684\:72b6\:6001 *)*)
(*{OPE[T,W],NO[T,W],NO[W,T],NO[Derivative[1][W],Derivative[1][T]]}*)


(* ::Input:: *)
(*(* \:786e\:5b9a Subscript[\[ScriptCapitalW], 3] \:4ee3\:6570 OPE *)*)
(*\[CapitalLambda]=NO[T,T]-3/10 Derivative[1][Derivative[1][T]]*)
(*OPE[T,T]=MakeOPE[(c One[w])/(2(z-w)^4)+(2T[w])/(z-w)^2+Derivative[1][T][w]/(z-w)+Ord[z,w,0]]*)
(*OPE[T,W]=MakeOPE[(3W[w])/(z-w)^2+Derivative[1][W][w]/(z-w)+Ord[z,w,0]]*)
(*OPE[W,W]=MakeOPE[(c One[w])/(z-w)^6+(2T[w])/(z-w)^4+Derivative[1][T][w]/(z-w)^3+1/(z-w)^2 (2 \[Beta] \[CapitalLambda][w]+3/10 Derivative[1][Derivative[1][T]][w])*)
(*+1/(z-w) (\[Beta] Derivative[1][\[CapitalLambda]][w]+1/15 Derivative[1][Derivative[1][Derivative[1][T]]][w])+Ord[z,w,0]]*)
(**)


(* ::Input:: *)
(*(* \:8bbe\:5b9a OPE \:4e4b\:540e\:ff0c\:8ba1\:7b97 NO *) *)
(*{NO[T,W],NO[Derivative[1][T],W],NO[W,Derivative[1][T]],NO[Derivative[1][W],Derivative[1][T]],NO[NO[T,W],NO[T,W]]}//Simplify*)
(*{OPE[T,NO[T,T]-3/10 Derivative[1][Derivative[1][T]]],OPE[W,NO[T,T]-3/10 Derivative[1][Derivative[1][T]]]}//Simplify*)
(**)


(* ::Title::Closed:: *)
(*\[ScriptCapitalN]=4 small SCFA from reflection group*)


(* ::Input:: *)
(*(* \:6ce8\:518c bc \[Beta]\[Gamma]  *)*)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)
(**)


(* ::Input:: *)
(*(* \:6ce8\:518c OPE *)*)
(*OPE[b,c]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]];*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[-(One[w]/(z-w))+Ord[z,w,0]];*)


(* ::Input:: *)
(*(* \[ScriptCapitalN]=4 small SCFA \:7684 bc\[Beta]\[Gamma] \:5b9e\:73b0 *)*)
(*Jplus=\[Beta];*)
(*J0=NO[b,c]+2NO[\[Beta],\[Gamma]];*)
(*Jminus=-NO[\[Beta],NO[\[Gamma],\[Gamma]]]-NO[\[Gamma],NO[b,c]]+3/2 Derivative[1][\[Gamma]];*)
(**)
(*Gplus=b;*)
(*Gminus=NO[b,\[Gamma]];*)
(*Gtplus=NO[c,Derivative[1][\[Beta]]]+2NO[Derivative[1][c],\[Beta]];*)
(*Gtminus=-NO[b,Derivative[1][c],c]+2NO[\[Beta],\[Gamma],Derivative[1][c]]+NO[Derivative[1][\[Beta]],\[Gamma],c]-3/2 Derivative[1][Derivative[1][c]];*)
(**)
(*T=-3/2 NO[b,Derivative[1][c]]-NO[\[Beta],Derivative[1][\[Gamma]]]-1/2 NO[Derivative[1][b],c];*)


(* ::Input:: *)
(*OPE[Gminus,\[Beta]]//OPEToSeries*)
(*OPE[Gtplus,\[Beta]]//OPEToSeries*)


(* ::Input:: *)
(*OPE[Gy,Jplus]*)
(*OPE[Gy,J0 y]*)


(* ::Input:: *)
(*Jy=Jplus+J0 y+Jminus y^2;*)
(*Gy=Gplus+Gminus y;*)
(*Gty=Gtplus+Gtminus y;*)
(*OPE[Gy,Jy]//OPEToSeries//Simplify*)
(*OPE[Gty,Jy]//OPEToSeries//Simplify*)
(*OPE[Gty,OPE[Gy,Jy]]*)


(* ::Input:: *)
(*(OPE[J0,Jplus]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*(2Jplus)/(z-w)//Expand;*)
(*%-%%//Simplify*)
(**)
(*(OPE[J0,Jminus]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*-((2Jminus)/(z-w))//Expand;*)
(*%-%%//Simplify*)
(**)
(*(OPE[Jplus,Jminus]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*((-3/2) One)/(z-w)^2+J0/(z-w);*)
(*%-%%//Simplify*)


(* ::Input:: *)
(*T=-3/2 NO[b,Derivative[1][c]]-NO[\[Beta],Derivative[1][\[Gamma]]]-1/2 NO[Derivative[1][b],c];*)
(*(* check TT OPE *)*)
(*(OPE[T,T]//OPEToSeries)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*((-9/2)One)/(z-w)^4+(2T)/(z-w)^2+Derivative[1][T]/(z-w);*)
(*%-%%//Simplify*)
(*(* check Sugawara condition *)*)
(*1/(2(-3/2+2)) ( NO[Jplus,Jminus]+ NO[Jminus,Jplus]+1/2 NO[J0,J0])-T//Simplify//Collect[#,NO[x__]]&*)


(* ::Input:: *)
(*(* Jminus \:662f super-Virasoro primary\:ff0c\:5373 Subscript[{GJ^-}, n\[GreaterEqual]2]=Subscript[{Gbar J^-}, n\[GreaterEqual]2]=0 *)*)
(*OPE[Gplus,Jminus]//OPEToSeries//Simplify*)
(*OPE[Gtplus,Jminus]//OPEToSeries//Simplify*)
(*OPE[Gminus,Jminus]//OPEToSeries//Simplify*)
(*OPE[Gtminus,Jminus]//OPEToSeries//Simplify*)


(* ::Title::Closed:: *)
(*rank-1 \[ScriptCapitalN]=3 theory*)


(* ::Input:: *)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)


(* ::Input:: *)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)
(*OPE[b,c]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]];*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[-(One[w]/(z-w))+Ord[z,w,0]];*)


(* ::Input:: *)
(*Clear[h]*)
(*{h[b],h[c],h[\[Beta]],h[\[Gamma]]}={2,-1,3/2,-(1/2)};*)
(*{h[T],h[\[ScriptCapitalJ]],h[\[ScriptCapitalG]],h[\[ScriptCapitalG]bar],h[W],h[Wbar]}={2,1,3/2,3/2,3/2,3/2};*)
(*h[op_Plus]:=h/@op*)
(*h[\[Lambda]_?NumericQ op_]:=\[Lambda] h[op]*)
(*h[Derivative[n_][op_]]:=h[op]+n;*)
(*h[NO[op___]]:=h/@List[op]//Total*)


(* ::Input:: *)
(*(* h[\[Beta]]=3/2, h[\[Gamma]]=-(1/2), h[b]=2, h[c]=-1 *)*)
(*\[ScriptCapitalJ]=(p NO[\[Beta],\[Gamma]]+(p-1)NO[b,c])/2;*)
(*(* h[\[ScriptCapitalG]]=h[\[ScriptCapitalG]bar]=3/2 *)*)
(*\[ScriptCapitalG]=NO[b,\[Gamma]];*)
(*\[ScriptCapitalG]bar=p NO[\[Beta],Derivative[1][c]]+(p-1)NO[Derivative[1][\[Beta]],c];*)
(*T=-(1/2)p NO[\[Beta],Derivative[1][\[Gamma]]]+(1-1/2 p)NO[Derivative[1][\[Beta]],\[Gamma]]-1/2 (p+1)NO[b,Derivative[1][c]]+1/2 (1-p)NO[Derivative[1][b],c];*)
(*(* h[W]=3/2 *)*)
(*W=\[Beta];*)
(*(* h[Wbar] = 3/2 *)*)
(*Wbar=Array[aaa,8] . {Derivative[1][Derivative[1][\[Gamma]]],*)
(*NO[\[Gamma],b,Derivative[1][c]],NO[\[Gamma],Derivative[1][b],c],NO[Derivative[1][\[Gamma]],b,c],*)
(*NO[\[Beta],\[Gamma],Derivative[1][\[Gamma]]],NO[Derivative[1][\[Beta]],\[Gamma],\[Gamma]],*)
(*NO[\[Beta],\[Gamma],\[Gamma],b,c],*)
(*NO[\[Beta],\[Beta],\[Gamma],\[Gamma],\[Gamma]]*)
(*}/.aaa[8]->1/.{aaa[1]->10/9,aaa[2]->-(4/3),aaa[3]->2/3,aaa[4]->-(8/3),aaa[5]->-4,aaa[6]->2/3,aaa[7]->2};*)
(**)
(*(* \:8ba1\:7b97\:989d\:5916\:7684\:751f\:6210\:5143 \[ScriptCapitalG]W, h[\[ScriptCapitalG]W] = 2 *)*)
(*OPE[\[ScriptCapitalG],W]//OPEToSeries*)
(*\[ScriptCapitalG]W=b;*)
(*\[ScriptCapitalG]barWbar=-(8/3) NO[b,NO[c^\[Prime]\[Prime],c]]-3 NO[\[Beta],NO[\[Beta],NO[\[Gamma],NO[\[Gamma],Derivative[1][c]]]]]+4 NO[\[Beta],NO[\[Gamma],NO[b,NO[Derivative[1][c],c]]]]+4 NO[\[Beta],NO[\[Gamma],c^\[Prime]\[Prime]]]+4 NO[\[Beta],NO[Derivative[1][\[Gamma]],Derivative[1][c]]]+2/3 NO[Derivative[1][b],NO[Derivative[1][c],c]]-2 NO[Derivative[1][\[Beta]],NO[\[Beta],NO[\[Gamma],NO[\[Gamma],c]]]]+8/3 NO[Derivative[1][\[Beta]],NO[Derivative[1][\[Gamma]],c]]-2/3 NO[\[Beta]^\[Prime]\[Prime],NO[\[Gamma],c]]-(10 \!\(\*SuperscriptBox[\(c\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\))/9;*)


(* ::Input:: *)
(*NO[\[ScriptCapitalG],W]//Expand*)


(* ::Input:: *)
(*OPE[T,Derivative[1][T]]/.{p->3}//Simplify//OPEToSeries*)
(*-2 (4 NO[b,Derivative[1][c]][w]+3 NO[\[Beta],Derivative[1][\[Gamma]]][w]+2 NO[Derivative[1][b],c][w]+NO[Derivative[1][\[Beta]],\[Gamma]][w])//Expand*)


(* ::Input:: *)
(*OPE[NO[b,c],NO[b,c]]//OPEToSeries*)
(*OPE[NO[\[Beta],\[Gamma]],NO[\[Beta],\[Gamma]]]//OPEToSeries*)
(*OPE[NO[\[Beta],\[Gamma]],NO[b,c]]//OPEToSeries*)


(* ::Input:: *)
(*OPE[2 NO[b,c]+3 NO[\[Beta],\[Gamma]],2 NO[b,c]+3 NO[\[Beta],\[Gamma]]]*)
(*OPE[2\[ScriptCapitalJ],2\[ScriptCapitalJ]]/.p->1*)


(* ::Input:: *)
(*2\[ScriptCapitalJ]/.p->3//Expand*)
(*OPE[2NO[b,c],2NO[b,c]]+OPE[3 NO[\[Beta],\[Gamma]],3 NO[\[Beta],\[Gamma]]]*)


(* ::Input:: *)
(*(* W \:662f osp(2|2) super-Virasoro primary \:4e14\:662f chiral (h[W]=m[W]) *)*)
(*OPE[\[ScriptCapitalG],W]//OPEToSeries//Simplify*)
(*OPE[\[ScriptCapitalG]bar,W]/.p->3//OPEToSeries//Simplify*)
(**)
(*(* Wbar \:662f osp(2|2) super-Virasoro primary \:4e14\:662f anti-chiral (h[Wbar]=-m[Wbar]) *)*)
(*OPE[\[ScriptCapitalG],Wbar]//OPEToSeries//Simplify*)
(*OPE[\[ScriptCapitalG]bar,Wbar]/.p->3//OPEToSeries//Simplify*)


(* ::Input:: *)
(*(* Wbar Wbar regular *)*)
(*OPE[Wbar,Wbar]//OPEToSeries//Simplify*)
(**)
(*(* W Wbar \:7684 singular \:90e8\:5206\:53ea\:6709 One \:7684 super-Virasoro descendants *)*)
(*-9/20 OPE[W,Wbar]//OPEToSeries//Simplify*)
(*2{One,(3p)/cc \[ScriptCapitalJ],p/cc T+(3p(3p-1))/(2cc(cc-1)) (NO[\[ScriptCapitalJ],\[ScriptCapitalJ]]-2/3 T)}/.cc->-3(2p-1)/.{p->3}//Simplify*)


(* ::Input:: *)
(*(-4 NO[b,Derivative[1][c]][w]+9 NO[\[Beta],NO[\[Beta],NO[\[Gamma],\[Gamma]]]][w]+2 (6 NO[\[Beta],NO[\[Gamma],NO[b,c]]][w]-6 NO[\[Beta],Derivative[1][\[Gamma]]][w]+NO[Derivative[1][b],c][w]+2 NO[Derivative[1][\[Beta]],\[Gamma]][w])//Expand)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*12 NO[b,Derivative[1][c]]+9 NO[\[Beta],NO[\[Beta],NO[\[Gamma],\[Gamma]]]]+12 NO[\[Beta],NO[\[Gamma],NO[b,c]]]+3 NO[\[Beta],Derivative[1][\[Gamma]]]+12 NO[Derivative[1][b],c]+13 NO[Derivative[1][\[Beta]],\[Gamma]]//Expand*)


(* ::Input:: *)
(*primary=NO[\[Beta],\[Gamma],c,\[Gamma]];*)
(*(* Subscript[L, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[T,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG], n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG],primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG]bar, n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG]bar,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalJ], n\[GreaterEqual]1]c=0 *)*)
(*(OPE[\[ScriptCapitalJ],primary]//OPEToSeries)/.p->3*)
(*(* Subscript[W, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[W,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[Wbar, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[Wbar,primary]//OPEToSeries)/.p->3*)
(**)
(*(OPE[\[ScriptCapitalG]W,primary]//OPEToSeries)/.p->3//Simplify*)
(*(OPE[\[ScriptCapitalG]Wbar,primary]//OPEToSeries)/.p->3//Simplify*)


(* ::Input:: *)
(*primary=\[Gamma];*)
(*(* Subscript[L, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[T,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG], n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG],primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG]bar, n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG]bar,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalJ], n\[GreaterEqual]1]c=0 *)*)
(*(OPE[\[ScriptCapitalJ],primary]//OPEToSeries)/.p->3*)
(*(* Subscript[W, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[W,primary]//OPEToSeries)/.p->3*)
(*(* Subscript[Wbar, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[Wbar,primary]//OPEToSeries)/.p->3*)
(**)
(*(OPE[\[ScriptCapitalG]W,primary]//OPEToSeries)/.p->3*)
(*(OPE[\[ScriptCapitalG]Wbar,primary]//OPEToSeries)/.p->3*)


(* ::Input:: *)
(*(* Subscript[L, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[T,c]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG], n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG],c]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalG]bar, n\[GreaterEqual]1/2]c=0 *)*)
(*(OPE[\[ScriptCapitalG]bar,c]//OPEToSeries)/.p->3*)
(*(* Subscript[\[ScriptCapitalJ], n\[GreaterEqual]1]c=0 *)*)
(*(OPE[\[ScriptCapitalJ],c]//OPEToSeries)/.p->3*)
(*(* Subscript[W, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[W,c]//OPEToSeries)/.p->3*)
(*(* Subscript[Wbar, n\[GreaterEqual]1]c=0 *)*)
(*(OPE[Wbar,c]//OPEToSeries)/.p->3*)


(* ::Input:: *)
(*OPE[\[ScriptCapitalJ],\[Gamma]]//OPEToSeries*)
(*(OPE[T,\[Gamma]]//OPEToSeries)/.p->3*)
(*(OPE[\[ScriptCapitalG],\[Gamma]]//OPEToSeries)/.p->3*)
(*(OPE[\[ScriptCapitalG]bar,\[Gamma]]//OPEToSeries)/.p->3*)
(*(OPE[\[ScriptCapitalG]W,\[Gamma]]//OPEToSeries)/.p->3*)
(*(OPE[W,\[Gamma]]//OPEToSeries)/.p->3*)
(*(OPE[Wbar,\[Gamma]]//OPEToSeries)/.p->3*)


(* ::Input:: *)
(*h[\[Beta]]=3/2;*)
(*h[\[Gamma]]=-(1/2);*)
(*h[b]=2;*)
(*h[c]=-1;*)
(*m[\[Beta]]=3/2;*)
(*m[\[Gamma]]=-3/2;*)
(*m[b]=1;*)
(*m[c]=-1;*)
(*h[pt]=1;*)
(*m[pt]=0;*)
(**)
(*Table[{{\[Beta]^n[\[Beta]],\[Gamma]^n[\[Gamma]],b^n[b],c^n[c],pt^n[pt]},Sum[n[\[ScriptCapitalO]]h[\[ScriptCapitalO]],{\[ScriptCapitalO],{\[Beta],\[Gamma],b,c,pt}}]==3/2,Sum[n[\[ScriptCapitalO]]m[\[ScriptCapitalO]],{\[ScriptCapitalO],{\[Beta],\[Gamma],b,c,pt}}]==-(3/2)},*)
(*{n[\[Beta]],0,10},{n[\[Gamma]],0,10},{n[b],0,1},{n[c],0,1},{n[pt],0,10}]/.{*)
(*{{__},False,_}->Null,*)
(*{{__},_,False}->Null*)
(*}//Flatten[#,4]&//DeleteCases[#,Null]&*)


(* ::Input:: *)
(*(* TT OPE *)*)
(*(OPE[T,T]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]/.p->3*)
(*((-3 (2p-1)/2)One)/(z-w)^4+(2T)/(z-w)^2+Derivative[1][T]/(z-w)/.p->3*)
(*%-%%//Simplify*)


(* ::Input:: *)
(*(* T\[ScriptCapitalJ] OPE *)*)
(*(OPE[T,\[ScriptCapitalJ]]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*(\[ScriptCapitalJ]/(z-w)^2+Derivative[1][\[ScriptCapitalJ]]/(z-w))/.{p->3};*)
(*%-%%//Simplify*)


(* ::Input:: *)
(*(* T\[Beta] OPE *)*)
(*(OPE[T,\[Beta]]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*((3/2 \[Beta])/(z-w)^2+Derivative[1][\[Beta]]/(z-w))/.{p->3};*)
(*%-%%//Simplify*)
(**)
(*(* T\[Gamma] OPE *)*)
(*(OPE[T,\[Gamma]]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*((-(1/2)\[Gamma])/(z-w)^2+Derivative[1][\[Gamma]]/(z-w))/.{p->3};*)
(*%-%%//Simplify*)
(**)
(*(* Tb OPE; h[b]=2 *)*)
(*(OPE[T,b]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*((2 b)/(z-w)^2+Derivative[1][b]/(z-w))/.{p->3};*)
(*%-%%//Simplify*)
(**)
(*(* Tc OPE; h[c]=-1 *)*)
(*(OPE[T,c]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO];*)
(*(-c/(z-w)^2+Derivative[1][c]/(z-w))/.{p->3};*)
(*%-%%//Simplify*)


(* ::Input:: *)
(*(* T\[ScriptCapitalG]=(h \[ScriptCapitalG])/(z-w)^2+\[ScriptCapitalG]^\[Prime]/(z-w) *)*)
(*(OPE[T,\[ScriptCapitalG]]/.p->3//OPEToSeries)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*(h \[ScriptCapitalG])/(z-w)^2+Derivative[1][\[ScriptCapitalG]]/(z-w)/.p->3*)
(*%-%%/.h->3/2//Simplify*)


(* ::Input:: *)
(*(* h[W]=3/2 *)*)
(*(OPE[T,W]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]/.p->3;*)
(*(h W)/(z-w)^2+Derivative[1][W]/(z-w)/.{h->3/2};*)
(*%-%%*)
(*(* h[Wbar]=3/2 *)*)
(*(OPE[T,Wbar]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]/.p->3;*)
(*(3/2 Wbar)/(z-w)^2+Derivative[1][Wbar]/(z-w);*)
(*%-%%//Simplify//Series[#,{z,w,0}]&//Simplify*)


(* ::Input:: *)
(*(OPE[\[ScriptCapitalJ],T]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*(OPE[\[ScriptCapitalJ],\[Beta]]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*(OPE[\[ScriptCapitalJ],\[Gamma]]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*(OPE[\[ScriptCapitalJ],b]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(*(OPE[\[ScriptCapitalJ],c]//OPEToSeries//Normal)/.{p->3}/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)


(* ::Input:: *)
(*(* m[W]=3/2; h[W] = m[W] \[Equal]> chiral op *)*)
(*(OPE[\[ScriptCapitalJ],W]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]/.p->3;*)
(*(3/2 W)/(z-w);*)
(*%-%%//Simplify//Series[#,{z,w,0}]&//Simplify*)
(*(* m[Wbar]=-(3/2); h[W] = -m[W] \[Equal]> anti-chiral op *)*)
(*(OPE[\[ScriptCapitalJ],Wbar]//OPEToSeries//Normal)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]/.p->3;*)
(*(-3/2 Wbar)/(z-w);*)
(*%-%%//Simplify//Series[#,{z,w,0}]&//Simplify*)


(* ::Input:: *)
(*(* chiral op \:6ee1\:8db3 shortening condition Subscript[{\[ScriptCapitalG]barO}, 1]=0 *)*)
(*(* Subscript[\[ScriptCapitalG]bar, -(1/2)]W=b, Subscript[\[ScriptCapitalG]bar, 1/2]W=0 *)*)
(*OPE[\[ScriptCapitalG]bar,W]/.p->3//OPEToSeries*)
(*(* anti-chiral op  *)*)
(*OPE[\[ScriptCapitalG],Wbar]/.p->3//OPEToSeries*)
(**)
(*OPE[\[ScriptCapitalG],W]/.p->3//OPEToSeries*)
(*OPE[\[ScriptCapitalG]bar,Wbar]/.p->3//OPEToSeries*)


(* ::Input:: *)
(*Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]/.p->3//Expand*)


(* ::Input:: *)
(*Completion[\[ScriptCapitalJ]/.p->3,Completion[\[ScriptCapitalJ]/.p->3,\[ScriptCapitalJ]/.p->3][0]][0]//Expand*)


(* ::Input:: *)
(*OPE[T,\[ScriptCapitalJ]]//FullSimplify*)


(* ::Input:: *)
(*OPE[T,Completion[Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0],\[ScriptCapitalJ]][0]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*1/2 (-4 NO[b,Derivative[1][c]][w]+9 NO[\[Beta],NO[\[Beta],NO[\[Gamma],\[Gamma]]]][w]+12 NO[\[Beta],NO[\[Gamma],NO[b,c]]][w]-9 NO[\[Beta],Derivative[1][\[Gamma]]][w]+4 NO[Derivative[1][b],c][w]+9 NO[Derivative[1][\[Beta]],\[Gamma]][w])/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]*)
(**)
(*-NO[b,Derivative[1][c]]+9/4 NO[\[Beta],NO[\[Beta],NO[\[Gamma],\[Gamma]]]]+3 NO[\[Beta],NO[\[Gamma],NO[b,c]]]-9/4 NO[\[Beta],Derivative[1][\[Gamma]]]+NO[Derivative[1][b],c]+9/4 NO[Derivative[1][\[Beta]],\[Gamma]]*)
(*%/%%//Simplify*)


(* ::Input:: *)
(*Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]//Expand*)


(* ::Input:: *)
(*OPE[T,\[ScriptCapitalJ]]//OPEToSeries*)
(*OPE[T,Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]]//OPEToSeries//Simplify*)
(*OPE[T,Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]//Expand*)
(*OPE[T,Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]]//OPEToSeries//FullSimplify*)


(* ::Input:: *)
(*Completion[W,Wbar][0]-1/27 Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]+2/9 Completion[\[ScriptCapitalG],\[ScriptCapitalG]bar][0]+4/9 Completion[\[ScriptCapitalJ],T][0]//Simplify*)


(* ::Input:: *)
(*Bracket[T,Completion[W,Wbar][0]][3]/.p->3*)
(*Bracket[T,Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][3]/.p->3*)


(* ::Input:: *)
(*Bracket[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]/.p->3//Simplify*)
(*Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]*)


(* ::Chapter:: *)
(*Null states*)


(* ::Title::Closed:: *)
(*General affine Lie algebra*)


(* ::Input:: *)
(*OPE[J[a_],J[b_]]:=Module[{c},MakeOPE[(k KK[a,b]One[w])/(z-w)^2+I f[a,b][c]*1/(z-w) J[c][w]+Ord[z,w,0]]]*)
(*OPE[T,J[b_]]:=MakeOPE[J[b][w]/(z-w)^2+1/(z-w) Derivative[1][J[b]][w]+Ord[z,w,0]]*)
(*OPE[T,T]:=MakeOPE[(c/2)/(z-w)^4+2/(z-w)^2 T[w]+1/(z-w) Derivative[1][T][w]+Ord[z,w,0]]*)


(* ::Input:: *)
(*(OPE[J[a],T]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[J[a],T]//OPEToSeries//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[J[a],T]//OPEToSeries//SeriesCoefficient[#,{z,w,-1}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)


(* ::Subsection:: *)
(*Sugawara stress tensor*)


(* ::Input:: *)
(*T1=1/(2(k+hcheck)) KKD[a$1,b$1]NO[J[a$1],J[b$1]];*)
(*T2=1/(2(k+hcheck)) KKD[a$2,b$2]NO[J[a$2],J[b$2]];*)
(**)
(*TTOPE=OPE[T1,T2];*)
(**)
(*TTOPE//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-4}]&//Expand//simplifyLieAlgebra//Simplify*)
(*TTOPE//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-3}]&//Expand//simplifyLieAlgebra//Simplify*)
(**)
(*(((TTOPE//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-2}]&//Expand//simplifyLieAlgebra//Simplify)/.{X_[w]:>X}//Expand)/.{KKD[a_,b_]NO[J[a_],J[b_]]:>KKD[aa,bb]NO[J[aa],J[bb]],KKD[b_,a_]NO[J[a_],J[b_]]:>KKD[aa,bb]NO[J[aa],J[bb]]}//simplifyLieAlgebra//Simplify//.{KKD[a_,b_]NO[J[a_],J[b_]][w]:>KKD[aa,bb]NO[J[aa],J[bb]][w],KKD[b_,a_]NO[J[a_],J[b_]][w]:>KKD[aa,bb]NO[J[aa],J[bb]][w]}//simplifyLieAlgebra)-2T1*)
(**)


(* ::Input:: *)
(*(OPE[J[a],T1]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[J[a],T1]//OPEToSeries//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[J[a],T1]//OPEToSeries//SeriesCoefficient[#,{z,w,-1}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(**)
(*(OPE[T1,J[a]]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[T1,J[a]]//OPEToSeries//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(*(OPE[T1,J[a]]//OPEToSeries//SeriesCoefficient[#,{z,w,-1}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)


(* ::Input:: *)
(*\[Chi]=NO[T,T]+\[Alpha]/2 Derivative[1][Derivative[1][T]]+\[Beta]/2 KKD[ap,bp]NO[Derivative[1][Derivative[1][J[ap]]],J[bp]]+\[Gamma] KKD[ap,bp]NO[Derivative[1][J[ap]],Derivative[1][J[bp]]];*)
(**)
(*((OPE[T,\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-4}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify)/.{KKD[aa_,bb_] NNO[J[aa_],J[bb_]]:>2(k+hcheck)T}//Simplify*)
(**)
(*((OPE[T,\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify//simplifyLieAlgebra//Simplify)/.{T->1/(2(k+hcheck)) KKD[ap,bp]NO[J[ap],J[bp]]}/.{NO->NNO}/.{NNO[J[ap],Derivative[1][J[bp]]]->NNO[Derivative[1][J[ap]],J[bp]]}//Collect[#,{NNO[Derivative[1][J[ap]],J[bp]],NNO[J[ap],Derivative[1][J[bp]]]},Simplify]&//Simplify*)
(**)


(* ::Input:: *)
(*\[Chi]=NO[T1,T2]+\[Alpha]/2 Derivative[1][Derivative[1][T1]]+\[Beta]/2 KKD[ap,bp]NO[Derivative[1][Derivative[1][J[ap]]],J[bp]]+\[Gamma] KKD[ap,bp]NO[Derivative[1][J[ap]],Derivative[1][J[bp]]];*)
(**)
(*((OPE[T1,\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify)//Simplify*)
(**)
(*((OPE[T1,\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-4}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify)//Simplify*)
(**)
(*((OPE[T1,\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-5}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify)//Simplify*)


(* ::Input:: *)
(*\[Chi]=NO[T1,T2]+\[Alpha]/2 Derivative[1][Derivative[1][T1]]+\[Beta]/2 KKD[ap,bp]NO[Derivative[1][Derivative[1][J[ap]]],J[bp]]+\[Gamma] KKD[ap,bp]NO[Derivative[1][J[ap]],Derivative[1][J[bp]]];*)
(**)
(*(OPE[J[a],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-4}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify*)
(**)
(*(OPE[J[a],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-3}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Simplify//simplifyLieAlgebra*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Solve[{(8+c+6 \[Alpha]+6 (hcheck+k) \[Beta])==0,*)
(*(3+5 \[Alpha]+3 hcheck \[Beta]+3 k \[Beta]+4 hcheck \[Gamma]+4 k \[Gamma]/.{(*c\[Rule](k dimg)/(k+hcheck)*)})==0},{\[Beta],\[Gamma]}]//Simplify*)
(*\[Beta]-2\[Gamma]/.%[[1]]//Simplify*)


(* ::Input:: *)
(*Solve[-((22 hcheck+(22+5 dimg) k)/(12 (hcheck+k)^2))==(-2k)/3,k]//Simplify;*)
(*%/.{hcheck->3,dimg->8}//FullSimplify*)


(* ::Input:: *)
(*((OPE[J[a],\[Chi]]//OPEToSeries//SeriesCoefficient[#,{z,w,-2}]&)/.{X_[w]:>X}/.{NO->NNO}//Expand//simplifyLieAlgebra//Expand)*)


(* ::Input:: *)
(*Out[135]/.{*)
(*KKD[ap_,bp_]NNO[J[a],NNO[J[ap_],J[bp_]]]:>KKD[aa,bb]NNO[J[a],NNO[J[aa],J[bb]]],*)
(*KKD[ap_,bp_]NNO[J[ap_],NNO[J[a],J[bp_]]]:>KKD[aa,bb]NNO[J[aa],NNO[J[a],J[bb]]],*)
(*KKD[ap_,bp_]NNO[J[ap_],NNO[J[bp_],J[a]]]:>KKD[aa,bb]NNO[J[aa],NNO[J[bb],J[a]]],*)
(**)
(*NNO[Derivative[1][J[ap_]],J[bp_]]f[a][ap_,bp_]:>NNO[Derivative[1][J[aa]],J[bb]]f[a][aa,bb],*)
(*NNO[Derivative[1][J[ap_]],J[bp_]]f[a][bp_,ap_]:>-NNO[Derivative[1][J[aa]],J[bb]]f[a][aa,bb],*)
(*NNO[J[ap_],Derivative[1][J[bp_]]]f[a][ap_,bp_]:>NNO[J[aa],Derivative[1][J[bb]]]f[a][aa,bb],*)
(*NNO[J[ap_],Derivative[1][J[bp_]]]f[a][bp_,ap_]:>-NNO[J[aa],Derivative[1][J[bb]]]f[a][aa,bb],*)
(**)
(*f[a][ap_,bp_]f[ap_][cp_,dp_]NNO[J[cp_],NNO[J[dp_],J[bp_]]]:>f[a][aa,bb]f[aa][cc,dd]NNO[J[cc],NNO[J[dd],J[bb]]],*)
(*f[a][bp_,ap_]f[ap_][cp_,dp_]NNO[J[cp_],NNO[J[bp_],J[dp_]]]:>-f[a][aa,bb]f[aa][cc,dd]NNO[J[cc],NNO[J[bb],J[dd]]],*)
(*f[a][bp_,ap_]f[ap_][cp_,dp_]NNO[J[bp_],NNO[J[cp_],J[dp_]]]:>-f[a][aa,bb]f[aa][cc,dd]NNO[J[bb],NNO[J[cc],J[dd]]],*)
(*f[][ap_,bp_,cp_]f[a,cp_][dp_]NNO[J[ap_],NNO[J[bp_],J[dp_]]]:>f[][aa,bb,cc]f[a,cc][dd]NNO[J[aa],NNO[J[bb],J[dd]]]*)
(**)
(*}//Collect[#,{NNO[J[_],NNO[J[_],J[_]]],J[a]^\[Prime]\[Prime]},Simplify]&*)


(* ::Title::Closed:: *)
(*other*)


(* ::Input:: *)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)
(*OPE[b,c]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]];*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]];*)


(* ::Input:: *)
(*J=NO[\[Gamma],\[Beta]]*)
(*NO[J,J]*)
(*T=1/2 NO[\[Beta],Derivative[1][\[Gamma]]]-1/2 NO[\[Gamma],Derivative[1][\[Beta]]]*)


(* ::Input:: *)
(*T=-NO[b,Derivative[1][c]];*)
(*OPE[T,NO[b,c]]//OPEToSeries*)
(*OPE[T,One]//OPEToSeries*)


(* ::Input:: *)
(*ee=1/2 NO[b,b]*)


(* ::Input:: *)
(*NO[Derivative[1][\[Psi]],\[Psi]]*)
(*NO[\[Psi],Derivative[1][\[Psi]]]*)


(* ::Input:: *)
(*T=-1/2 NO[\[Psi],\[Psi]']*)
(*OPE[T,T]//OPEToSeries*)


(* ::Input:: *)
(*D[1/(z-w),z]D[1/(z-w),w]*)
(*-(1/(z-w))D[1/(z-w),z,w]*)
(*D[1/(z-w),z]*)
(*-D[1/(z-w),z,w]*)


(* ::Input:: *)
(*1/4*D[1/(z-w),w]*)


(* ::Input:: *)
(*-1/4 (z-w)D[1/(z-w),z,w]*)


(* ::Input:: *)
(*Bosonic[J[i_]]*)


(* ::Input:: *)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[bb,cc]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*Jp=\[Beta];*)
(*J3=1/2 (NO[bb,cc]-2NO[\[Beta],\[Gamma]]);*)
(*h=2J3;*)
(*Jm=-NO[\[Beta],NO[\[Gamma],\[Gamma]]]+NO[\[Gamma],NO[bb,cc]]-3/2 Derivative[1][\[Gamma]];*)
(*OPE[J3,Jp]*)
(*OPE[J3,Jm]*)
(*OPE[Jp,Jm]*)
(*OPE[h,h]*)


(* ::Input:: *)
(*OPE[Jm,cc]*)
(*OPE[Jm,-NO[\[Gamma],cc]]*)
(*OPE[Jm,2 NO[\[Gamma],NO[\[Gamma],cc]]]*)
(*OPE[Jm,-6 NO[\[Gamma],NO[\[Gamma],NO[\[Gamma],cc]]]]*)


(* ::Title::Closed:: *)
(*compare pyope*)


(* ::Input:: *)
(*Bosonic[\[CurlyPhi]A,\[CurlyPhi]B]*)
(*OPE[\[CurlyPhi]A,\[CurlyPhi]B]=MakeOPE[A[w]/(z-w)^2+B[w]/(z-w)+Ord[z,w,0]];*)
(*result=OPE[D[\[CurlyPhi]A],\[CurlyPhi]B];*)
(*{OPEPole[3][result],OPEPole[2][result]}*)
(*Derivative[1][\[CurlyPhi]A]*)


(* ::Title:: *)
(*Global definitions (Lie algebras)*)


(* ::Input:: *)
(*c=.;*)
(*d=.;*)


(* ::Input::Initialization:: *)
Clear[simplifyLieAlgebra]
simplifyLieAlgebra[exp_]:=(exp//Expand)//.{
KK[aa_,bb_]KKD[cc_,dd_]:>\[Delta]@@Complement[Union[{aa,bb},{cc,dd}],Intersection[{aa,bb},{cc,dd}]]/;Length[Intersection[{aa,bb},{cc,dd}]]==1,

KK[a_,b_]KKD[c_,d_]:>dimg/;Length[Intersection[{a,b},{c,d}]]==2,

\[Delta][a_,b_]x_:>(x/.{a->b})/;(!TrueQ[D[x,a]==0]),
\[Delta][a_,b_]x_:>(x/.{b->a})/;(!TrueQ[D[x,b]==0]),

\[Delta][a_,a_]:>dimg,
\[Delta][a_,b_]\[Delta][c_,d_]:>dimg/;(Intersection[{a,b},{c,d}]//Length)==2,
\[Delta][a_,b_]^2:>dimg,


(* f[a,b,c][] = f^abc, f[a,b][c] = Subscript[f^ab, c], f[a][b,c] = Subscript[f^a, bc], f[][a,b,c] = Subscript[f, abc] *)

KK[aa_,bb_]f[ii___][jj___]:>KfContraction[KK[aa,bb],f[ii][jj]]/;Length[Intersection[{aa,bb},{jj}]]>=1,
KKD[aa_,bb_]f[ii___][jj___]:>KfContraction[KKD[aa,bb],f[ii][jj]]/;Length[Intersection[{aa,bb},{ii}]]>=1,

(* Subscript[(f^a) , ab] = Subscript[(f^ab) , b] = 0*)
f[a_,b_][c_]:>0/;Length[Intersection[{a,b},{c}]]>0,
f[a_][b_,c_]:>0/;Length[Intersection[{a},{b,c}]]>0,

f[dd_,ee_,ff_][]:>Signature[{dd,ee,ff}](f@@Sort[{dd,ee,ff}])[]/;!TrueQ[{dd,ee,ff}==Sort[{dd,ee,ff}]],
f[][dd_,ee_,ff_]:>Signature[{dd,ee,ff}](f[]@@Sort[{dd,ee,ff}])/;!TrueQ[{dd,ee,ff}==Sort[{dd,ee,ff}]],

(* Subscript[f^ab, c] Subscript[f^d, ab] = 2h^vSubscript[\[Delta], c]^d *)
f[i___][j___]f[k___][l___]:>ffContraction[f[i][j],f[k][l]]/;Length[Intersection[{i}~Join~{j},{k}~Join~{l}]]==2,

(* f^abc Subscript[f, abc] = 2h^v dimg *)
f[i___][j___]f[k___][l___]:>ffContraction[f[i][j],f[k][l]]/;Length[Intersection[{i}~Join~{j},{k}~Join~{l}]]==3,


(* fff relations, Subscript[f^a, bc] Subscript[f^ba^\[Prime], d]Subscript[f^cd, e] = h^vSubscript[f^aa^\[Prime], c] *)
f[i1___][j1___]f[i2___][j2___]f[i3___][j3___]:>fffContraction[f[i1][j1],f[i2][j2],f[i3][j3]]/; Length[Intersection[{i1}~Join~{j1},{i2}~Join~{j2}]]==1&&Length[Intersection[{i3}~Join~{j3},{i2}~Join~{j2}]]==1&&Length[Intersection[{i1}~Join~{j1},{i3}~Join~{j3}]]==1,


f[a_][b_,c_]NNO[J[b_],J[c_]]:>I hcheck Derivative[1][J[a]],
f[a_][b_,c_]NNO[J[c_],J[b_]]:>-I hcheck Derivative[1][J[a]],

f[][d_,e_,f_]NNO[J[a_],NNO[J[b_],J[c_]]]:>Signature[{d,e,f}]Signature[{a,b,c}](I hcheck KKD[cc,dd]NNO[Derivative[1][J[dd]],J[cc]]/.{cc->c[RandomInteger[10000]],dd->d[RandomInteger[10000]]})/;Length[Intersection[{d,e,f},{a,b,c}]]==3

};

Clear[ffContraction]


(* full contraction of ff *)
ffContraction[f[i___][j___],f[k___][l___]]:=Module[{},
Signature[{i}~Join~{j}]Signature[{k}~Join~{l}]2hcheck  dimg
]/;Length[Intersection[{i}~Join~{j},{k}~Join~{l}]]==3;


KfContraction[KK[a_,b_],f[i___][j___]]:=Module[{raised,raisedTo},
raised=Intersection[{a,b},{j}][[1]];
raisedTo=Complement[{a,b},{raised}][[1]];
Signature[{j}]Signature[{raised}~Join~Complement[{j},{raised}]]f[i,raisedTo]@@Complement[{j},{raised}]
];

KfContraction[KKD[a_,b_],f[i___][j___]]:=Module[{lowered,loweredTo},
lowered=Intersection[{a,b},{i}][[1]];
loweredTo=Complement[{a,b},{lowered}][[1]];
Signature[{i}]Signature[Complement[{i},{lowered}]~Join~{lowered}](f@@Complement[{i},{lowered}])@@({loweredTo}~Join~{j})
];


fffContraction[f[i1___][j1___],f[i2___][j2___],f[i3___][j3___]]:=Module[{up1,down1,up2,down2,up3,down3,remain1,remain2,remain3,
inter12,inter13,inter23,up,down
},
up1={i1};
down1={j1};
up2={i2};
down2={j2};
up3={i3};
down3={j3};
inter12=Intersection[{i1}~Join~{j1},{i2}~Join~{j2}];
inter13=Intersection[{i1}~Join~{j1},{i3}~Join~{j3}];
inter23=Intersection[{i2}~Join~{j2},{i3}~Join~{j3}];

remain1=Complement[up1~Join~down1,inter12~Join~inter13];
remain2=Complement[up2~Join~down2,inter12~Join~inter23];
remain3=Complement[up3~Join~down3,inter13~Join~inter23];

up={};
down={};
If[MemberQ[up1,remain1[[1]]],up=up~Join~remain1,down=down~Join~remain1];
If[MemberQ[up2,remain2[[1]]],up=up~Join~remain2,down=down~Join~remain2];
If[MemberQ[up3,remain3[[1]]],up=up~Join~remain3,down=down~Join~remain3];

Signature[up1~Join~down1]Signature[up2~Join~down2]Signature[up3~Join~down3]Signature[inter13~Join~inter12~Join~remain1]Signature[inter23~Join~remain2~Join~inter12]Signature[remain3~Join~inter13~Join~inter23]Signature[up~Join~down]Signature[remain1~Join~remain2~Join~remain3]*(-hcheck)(f@@up)@@down

]/; Length[Intersection[{i1}~Join~{j1},{i2}~Join~{j2}]]==1&&Length[Intersection[{i3}~Join~{j3},{i2}~Join~{j2}]]==1&&Length[Intersection[{i1}~Join~{j1},{i3}~Join~{j3}]]==1


(* partial contraction of ff *)
ffContraction[f[i___][j___],f[k___][l___]]:=Module[{up1,down1,up2,down2,intersection,remain1,remain2,$m},
up1={i};
down1={j};
up2={k};
down2={l};

intersection=Intersection[{i}~Join~{j},{k}~Join~{l}];
remain1=Complement[{i}~Join~{j},intersection][[1]];
remain2=Complement[{k}~Join~{l},intersection][[1]];
If[MemberQ[up1,remain1]&&MemberQ[up2,remain2],$m=KK[remain1,remain2]];
If[MemberQ[down1,remain1]&&MemberQ[down2,remain2],$m=KKD[remain1,remain2]];
If[MemberQ[up1,remain1]&&MemberQ[down2,remain2]||MemberQ[down1,remain1]&&MemberQ[up2,remain2],$m=\[Delta][remain1,remain2]];
2hcheck Signature[{i}~Join~{j}]Signature[{k}~Join~{l}]Signature[intersection~Join~{remain1}]Signature[intersection~Join~{remain2}]$m

]/;Length[Intersection[{i}~Join~{j},{k}~Join~{l}]]==2;


(* ::Input:: *)
(*Clear[simplifyLieAlgebra2]*)
(*simplifyLieAlgebra2[exp_]:=exp//.{*)
(*KK[aa_,bb_]KKD[cc_,dd_]:>\[Delta]@@Complement[Union[{aa,bb},{cc,dd}],Intersection[{aa,bb},{cc,dd}]]/;Length[Intersection[{aa,bb},{cc,dd}]]==1,*)
(**)
(*KK[a_,b_]KKD[c_,d_]:>dimg/;Length[Intersection[{a,b},{c,d}]]==2,*)
(**)
(*\[Delta][a_,b_]x_:>(x/.{a->b})/;(!TrueQ[D[x,a]==0]),*)
(*\[Delta][a_,b_]x_:>(x/.{b->a})/;(!TrueQ[D[x,b]==0]),*)
(**)
(*\[Delta][a_,a_]:>dimg,*)
(*\[Delta][a_,b_]\[Delta][c_,d_]:>dimg/;(Intersection[{a,b},{c,d}]//Length)==2,*)
(*\[Delta][a_,b_]^2:>dimg,*)
(**)
(**)
(*(* f[a,b,c][] = f^abc, f[a,b][c] = Subscript[f^ab, c], f[a][b,c] = Subscript[f^a, bc], f[][a,b,c] = Subscript[f, abc] *)*)
(*(* f^abc Subscript[f, abc] = 2h^vdimg*)*)
(*f[a_,b_,c_][]f[][d_,e_,f_]:>Signature[{a,b,c}]Signature[{d,e,f}]2hcheck  dimg/;Length[Intersection[{a,b,c},{d,e,f}]]==3,*)
(**)
(**)
(*(* K^abSubscript[f^cd, b] = f^cda *)*)
(*KK[a_,b_]f[c_,d_][e_]:>(f[c,d,Complement[{a,b},Intersection[{a,b},{e}]][[1]]][])/;Length[Intersection[{a,b},{e}]]==1,*)
(*(* Subscript[K, ab]Subscript[f^cb, e] = Subscript[f^c, ae] *)*)
(*KKD[a_,b_]f[c_,d_][e_]:>(f[c][Complement[{a,b},Intersection[{a,b},{d}]][[1]],e])/;Length[Intersection[{a,b},{d}]]==1,*)
(*(* Subscript[K, ab]Subscript[f^bc, e] = - Subscript[K, ab]Subscript[f^cb, e] = -Subscript[f^c, ae]*)*)
(*KKD[a_,b_]f[d_,c_][e_]:>-KKD[a,b]f[c,d][e]/;Length[Intersection[{a,b},{d}]]==1,*)
(*(* Subscript[K, ab]Subscript[f^b, de] = Subscript[f, ade] *)*)
(*KKD[a_,b_]f[c_][d_,e_]:>(f[][Complement[{a,b},Intersection[{a,b},{c}]][[1]],d,e])/;Length[Intersection[{a,b},{c}]]==1,*)
(**)
(*KKD[a_,b_]f[c_,d_,e_][]:>f[c,d][Complement[{a,b,e},Intersection[{a,b},{e}]][[1]]]/;Length[Intersection[{a,b},{e}]]==1,*)
(*KKD[a_,b_]f[c_,d_,e_][]:>-f[c,e][Complement[{a,b,d},Intersection[{a,b},{d}]][[1]]]/;Length[Intersection[{a,b},{d}]]==1,*)
(*KKD[a_,b_]f[c_,d_,e_][]:>f[d,e][Complement[{a,b,c},Intersection[{a,b},{c}]][[1]]]/;Length[Intersection[{a,b},{c}]]==1,*)
(**)
(**)
(**)
(*f[a_,b_][c_]:>0/;Length[Intersection[{a,b},{c}]]>0,*)
(*f[a_][b_,c_]:>0/;Length[Intersection[{a},{b,c}]]>0,*)
(**)
(*f[dd_,ee_,ff_][]:>Signature[{dd,ee,ff}]Signature[Sort[{dd,ee,ff}]](f@@Sort[{dd,ee,ff}])[]/;!TrueQ[{dd,ee,ff}==Sort[{dd,ee,ff}]],*)
(*f[][dd_,ee_,ff_]:>Signature[{dd,ee,ff}]Signature[Sort[{dd,ee,ff}]](f[]@@Sort[{dd,ee,ff}])/;!TrueQ[{dd,ee,ff}==Sort[{dd,ee,ff}]],*)
(**)
(*(* Subscript[f^ab, c] Subscript[f^d, ab] = 2h^vSubscript[\[Delta], c]^d *)*)
(*f[aa_,bb_][cc_]f[dd_][ee_,ff_]:>2hcheck Signature[{aa,bb}]Signature[{ee,ff}]\[Delta][cc,dd]/;Length[Intersection[{aa,bb},{ee,ff}]]==2&&!TrueQ[cc==dd],*)
(**)
(*f[a_,b_][c_]f[c_][a_,b_]:>f[a,b,c][]f[][c,a,b],*)
(*f[a_,b_][c_]f[c_][b_,a_]:>-f[a,b,c][]f[][c,a,b],*)
(**)
(*f[a_][b_,c_]f[b_][a_,d_]:>-2hcheck KKD[c,d],*)
(*f[a_][c_,b_]f[b_][a_,d_]:>2hcheck KKD[c,d],*)
(*f[a_][b_,c_]f[b_][d_,a_]:>2hcheck KKD[c,d],*)
(*f[a_][c_,b_]f[b_][d_,a_]:>-2hcheck KKD[c,d],*)
(**)
(**)
(*f[a_][b_,c_]f[b_,d_][a_]:>2hcheck \[Delta][c,d],*)
(*f[a_][c_,b_]f[b_,d_][a_]:>-2hcheck \[Delta][c,d],*)
(*f[a_][b_,c_]f[d_,b_][a_]:>-2hcheck \[Delta][c,d],*)
(*f[a_][c_,b_]f[d_,b_][a_]:>2hcheck \[Delta][c,d],*)
(**)
(*f[aa_,bb_,cc_][]f[][dd_,bb_,cc_]:>2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[aa_,bb_,cc_][]f[][bb_,dd_,cc_]:>-2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[aa_,bb_,cc_][]f[][bb_,cc_,dd_]:>2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][aa_,bb_,cc_]f[dd_,bb_,cc_][]:>2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][aa_,bb_,cc_]f[bb_,dd_,cc_][]:>-2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][aa_,bb_,cc_]f[bb_,cc_,dd_][]:>2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][bb_,aa_,cc_]f[dd_,bb_,cc_][]:>-2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][bb_,aa_,cc_]f[bb_,dd_,cc_][]:>2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(*f[][bb_,aa_,cc_]f[bb_,cc_,dd_][]:>-2 hcheck \[Delta][aa,dd]/;!TrueQ[aa==dd],*)
(**)
(**)
(**)
(**)
(*f[][aa_,bb_,cc_] f[bb_,cc_][dd_]:>2hcheck KKD[aa,dd]/;!TrueQ[aa==dd],*)
(*f[aa_,bb_,cc_][] f[dd_][bb_,cc_]:>2hcheck KK[aa,dd]/;!TrueQ[aa==dd],*)
(**)
(**)
(**)
(*(* fff relations, Subscript[f^a, bc] Subscript[f^ba^\[Prime], d]Subscript[f^cd, e] = h^vSubscript[f^aa^\[Prime], c] *)*)
(*f[a_][b_,c_] f[aa_,b_][d_] f[d_,c_][e_]:>-hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][b_,c_] f[b_,aa_][d_] f[d_,c_][e_]:>+hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][b_,c_] f[b_,aa_][d_] f[c_,d_][e_]:>-hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][b_,c_] f[aa_,b_][d_] f[c_,d_][e_]:>+hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(**)
(*f[a_][c_,b_] f[aa_,b_][d_] f[d_,c_][e_]:>+hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_,aa_][d_] f[d_,c_][e_]:>-hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_,aa_][d_] f[c_,d_][e_]:>+hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[aa_,b_][d_] f[c_,d_][e_]:>-hcheck f[a,aa][e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(**)
(**)
(*f[a_][b_,c_] f[b_][aa_,d_] f[d_,c_][e_]:>+hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][b_,c_] f[b_][aa_,d_] f[c_,d_][e_]:>-hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_][aa_,d_] f[d_,c_][e_]:>-hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_][aa_,d_] f[c_,d_][e_]:>+hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(**)
(*f[a_][b_,c_] f[b_][d_,aa_] f[d_,c_][e_]:>-hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][b_,c_] f[b_][d_,aa_] f[c_,d_][e_]:>+hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_][d_,aa_] f[d_,c_][e_]:>+hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(*f[a_][c_,b_] f[b_][d_,aa_] f[c_,d_][e_]:>-hcheck f[a][aa,e]/;!TrueQ[a==e]&&!TrueQ[aa==e],*)
(**)
(**)
(*f[a_][b_,c_]NNO[J[b_],J[c_]]:>I hcheck Derivative[1][J[a]],*)
(*f[a_][b_,c_]NNO[J[c_],J[b_]]:>-I hcheck Derivative[1][J[a]],*)
(**)
(*f[][d_,e_,f_]NNO[J[a_],NNO[J[b_],J[c_]]]:>Signature[{d,e,f}]Signature[{a,b,c}](I hcheck KKD[cc,dd]NNO[Derivative[1][J[dd]],J[cc]]/.{cc->c[RandomInteger[10000]],dd->d[RandomInteger[10000]]})/;Length[Intersection[{d,e,f},{a,b,c}]]==3*)
(*};*)


(* ::Title:: *)
(*Global definitions*)


(* ::Input::Initialization:: *)
Clear[Declare]
(* Here bosons refer to bosonic vertex fields *)
(* fermions refer to fermonic vertex fields *)
Declare[freefields_:{}][bosons_:{}][fermions_:{}]:=Module[{},
freeFundFields=freefields;
bosonicFields=bosons;
Bosonic@@bosons;
fermionicFields=fermions;
Fermionic@@fermions;
]


(* ::Input::Initialization:: *)
IsBosonicField[f_]:=Module[{pattern},
pattern=Table[{
bosonicFields[[i]][z_]:>True,
\!\(\*SuperscriptBox[\(bosonicFields\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[z_]:>True,
\[Lambda]_ bosonicFields[[i]][z_]:>True,
\[Lambda]_ \!\(\*SuperscriptBox[\(bosonicFields\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[z_]:>True,
a_. e1_[z_]+e2_:>(bosonicFields[e1[z]]||bosonicFields[a])&&bosonicFields[e2],
One[z_]:>True
},
{i,Length[bosonicFields]}]//Flatten;
(f//Expand)/.pattern//TrueQ
];
IsFermionicField[f_]:=Module[{pattern},
pattern=Table[{
fermionicFields[[i]][z_]:>True,
\!\(\*SuperscriptBox[\(fermionicFields\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[z_]:>True,
\[Lambda]_ fermionicFields[[i]][z_]:>True,
\[Lambda]_ \!\(\*SuperscriptBox[\(fermionicFields\[LeftDoubleBracket]i\[RightDoubleBracket]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[z_]:>True,
a_. e1_[z_]+e2_:>(fermionicFields[e1[z]]||fermionicFields[a])&&fermionicFields[e2]
},
{i,Length[fermionicFields]}]//Flatten;
(f//Expand)/.pattern//TrueQ
];

(*
Bosonic[\[CurlyPhi]];
Fermionic[\[Psi]]

{IsFermionicField[\[Psi][z]]\[Equal]True,
IsFermionicField[(\[Psi]^\[Prime])[z]]\[Equal]True,
IsFermionicField[3(\[Psi]^\[Prime])[z]]\[Equal]True,
IsFermionicField[3\[CurlyPhi][z]]\[Equal]False,
IsFermionicField[3(\[CurlyPhi]^\[Prime])[z]]\[Equal]False,
IsFermionicField[3(\[ExponentialE]^\[Delta])[z]]\[Equal]False
}*)



(* ::Input::Initialization:: *)
Clear[MakeVertexField]
MakeVertexField[exp_][z_]:=exp/.{NO[o___]:>NO[o][z],o_:>o[z]/;IsBosonicField[o[z]]||IsFermionicField[o[z]]}


(* ::Input::Initialization:: *)
Clear[OPEFull]
OPEFull[o1_+o2_,o3_][z_,w_][n_]:=OPEFull[o1,o3][z,w][n]+OPEFull[o2,o3][z,w][n];
OPEFull[o1_,o2_+o3_][z_,w_][n_]:=OPEFull[o1,o2][z,w][n]+OPEFull[o1,o3][z,w][n];
OPEFull[o1_,o2_][z_,w_][n_]:=Module[{ExtractVertexOps,V1,V2,LeadingTerm,LeadingPower,
freeOPE,leadingPowerOfFreeOPE,leadingPowerOfVertexOPE,vertexOPE,vertexFullOPE,ztemp},
ExtractVertexOps[o_]:=MakeField[o][ztemp]/.{f_[ztemp]:>1};
V1=ExtractVertexOps[o1];
V2=ExtractVertexOps[o2];
V1=V1/.ExtractVertexOps->Identity;
V2=V2/.ExtractVertexOps->Identity;
LeadingTerm[ope_,z]:=Normal[ope/.z->(z+O[z]^2)];
LeadingPower[ope_]:=z D[Log[LeadingTerm[ope,z]],z];
freeOPE=OPEFree[o1/V1,o2/V2][z,w][n];
vertexOPE=OPE[V1,V2]//OPEToSeries;
leadingPowerOfFreeOPE=LeadingPower[freeOPE/.{w->0}]//Simplify;
leadingPowerOfVertexOPE=LeadingPower[vertexOPE/.{w->0}]//Simplify;

vertexFullOPE=(OPE[V1,V2]//OPEToSeries//Normal)+Sum[(z-w)^nn/Factorial[nn] (NO[Derivative[nn][V1],V2]//MakeVertexField[#][w]&),{nn,0,-leadingPowerOfFreeOPE}];

OPEFree[o1/V1,o2/V2][z,w][-leadingPowerOfVertexOPE]*vertexFullOPE//Series[#,{z,w,n}]&
];


(* ::Input::Initialization:: *)
Clear[Bracket,Completion]
(* Subscript[{O1O2}, -n]=1/n!(\[PartialD]^nO1 O2)*)
Bracket[O1_,O2_][n_/;n<0]:=1/Abs[n]! NO[Derivative[Abs[n]][O1],O2]
Bracket[O1_,O2_][0]:=NO[O1,O2]
Bracket[O1_,O2_][n_/;n>0]:=(OPE[O1,O2]//OPEToSeries//Simplify//SeriesCoefficient[#,{z,w,-n}]&)/.\[ScriptCapitalO]_[w]:>\[ScriptCapitalO]
Completion[O1_,O2_][m_]:=Sum[
If[Pochhammer[2h[O1]+2h[O2]-2m-n-1,n]!=0,
((-1)^n Pochhammer[2h[O1]-n-m,n])/(n! Pochhammer[2h[O1]+2h[O2]-2m-n-1,n]),0]Derivative[n][Bracket[O1,O2][m+n]],
{n,0,2h[O1]-m}]


(* ::Input:: *)
(*OPE[T,Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]]//OPEToSeries//Simplify*)


(* ::Input:: *)
(*p=3;*)
(*Completion[W,Wbar][0]-1/27 Completion[\[ScriptCapitalJ],Completion[\[ScriptCapitalJ],\[ScriptCapitalJ]][0]][0]+2/9 Completion[\[ScriptCapitalG],\[ScriptCapitalG]bar][0]+4/9 Completion[\[ScriptCapitalJ],T][0]//Simplify*)


(* ::Input:: *)
(*Bosonic[\[Beta],\[Gamma]];*)
(*Fermionic[b,c];*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[-One[w]/(z-w)+Ord[z,w,0]]*)
(*OPE[b,c]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*Jplus=\[Beta];*)
(*J0=NO[b,c]+2NO[\[Beta],\[Gamma]];*)
(*Jminus=NO[\[Beta],NO[\[Gamma],\[Gamma]]]+NO[\[Gamma],NO[b,c]]-3/2 Derivative[1][\[Gamma]];*)
(*OPEFull[J0,Jplus][z,w][-1]*)
(*OPEFull[J0,Jminus][z,w][-1]-(-2Jminus//MakeVertexField[#][w]&)/(z-w)//FullSimplify*)
(*OPEFull[Jplus,Jminus][z,w][-1]-((3/2 One[w])/(z-w)^2+(-J0//MakeVertexField[#][w]&)/(z-w))//FullSimplify*)
(*T=-3/2 NO[b,Derivative[1][c]]-NO[\[Beta],Derivative[1][\[Gamma]]]-1/2 NO[Derivative[1][b],c];*)
(*OPEFull[T,T][z,w][-1]-(-((9 One[w])/(2 (z-w)^4))+(2T//MakeVertexField[#][w]&)/(z-w)^2+(Derivative[1][T]//MakeVertexField[#][w]&)/(z-w))//Simplify*)


(* ::Input:: *)
(*Bosonic[\[Beta],\[Gamma]]*)
(*OPE[\[Beta],\[Gamma]]=MakeOPE[One[w]/(z-w)+Ord[z,w,0]]*)


(* ::Input:: *)
(*OPE[NO[\[Beta],\[Gamma]],NO[\[Beta],Derivative[1][\[Gamma]]]]//OPEToSeries*)
(*OPE[NO[\[Gamma],\[Beta]],NO[\[Beta],Derivative[1][\[Gamma]]]]//OPEToSeries*)


(* ::Title::Closed:: *)
(*Abstract VOA*)


(* ::Input:: *)
(*Clear[Y,YM]*)
(*Y[a_,z_]:=Sum[a[n]z^(-n-h[a]),{m,-Infinity,Infinity}]//Quiet*)
(*YM[a_,z_]:=Sum[a[n,Text[M]]z^(-n-1),{m,-Infinity,Infinity}]//Quiet*)
(*ToPhysics[f_]:=f/.{a[n_,Text[M]]:>a[n-h[a]+1]}//Quiet*)
(*ToMath[f_]:=f/.{a[n_]:>a[n+h[a]-1,"M"]}//Quiet*)
(**)
(*YM[a,z]//ToPhysics*)
(*Y[a,z]//ToMath*)


(* ::Input:: *)
(*Clear[Curly,Borcherds]*)
(*CC[{*)
(* {-1},*)
(* {nn_}*)
(*}]:=(-1)^nn;*)
(*CC[{*)
(* {m_/;IntegerQ[m]&&m>=0},*)
(* {n_/;IntegerQ[n]&&n>=0}*)
(*}]:=Factorial[m]/(Factorial[m-n]Factorial[n]);*)
(*Borcherds[f_]:=f/.{*)
(*Curly[a_,b_,k_][l_]:>*)
(*Sum[CC[{*)
(* {k-1},*)
(* {nn}*)
(*}](-1)^nn (a[-h[a]+k-nn] . b[l+h[a]-k+nn]*)
(*+\[Epsilon][a,b](-1)^k b[l+h[a]-nn-1] . a[-h[a]+nn+1]),{nn,0,+Infinity}]*)
(*};*)
(*\[Epsilon][\[Psi],\[Psi]]=-1;*)
(*h[\[Psi]]=1/2;*)
(**)
(*(Curly[\[Psi],\[Psi],1][2]//Borcherds)/.{Infinity->100}//Expand*)
(*(Curly[\[Psi],\[Psi],1][1]//Borcherds)/.{Infinity->100}*)
(*(Curly[\[Psi],\[Psi],1][0]//Borcherds)/.{Infinity->100}//Expand*)
(*(Curly[\[Psi],\[Psi],1][-1]//Borcherds)/.{Infinity->100}//Expand*)


(* ::Input:: *)
(*h[a_[n_][b_]]:=h[b]-n*)
(*h[\[Psi][1/2][\[Psi]]]*)
(**)


(* ::Input:: *)
(*order=101;*)
(*n=0;*)
(*Sum[\[Psi][k] . \[Psi][n-1-k],{k,-(1/2),-(order/2),-1}]-Sum[\[Psi][n-1-k] . \[Psi][k],{k,1/2,order/2,1}]*)
(**)
(*n=1;*)
(*Sum[\[Psi][k] . \[Psi][n-1-k],{k,-(1/2),-(order/2),-1}]-Sum[\[Psi][n-1-k] . \[Psi][k],{k,1/2,order/2,1}]*)
(**)
(*n=-1;*)
(*Sum[\[Psi][k] . \[Psi][n-1-k],{k,-(1/2),-(order/2),-1}]-Sum[\[Psi][n-1-k] . \[Psi][k],{k,1/2,order/2,1}]*)


(* ::Input:: *)
(*Clear[J,c]*)
(*c[h_,j_,m_]:=SeriesCoefficient[Log[1+z]^m (1+z)^(h-1),{z,0,j}]*)
(*J[{m_}]:=(2\[Pi] I)^(-m-1) Sum[c[1,j+1-1,m+1-1]J[j],{j,m,m+5}]*)
(*J[{-1}]*)
(*J[{-2}]*)
(*J[{-3}]*)



