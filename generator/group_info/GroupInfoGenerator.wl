(* ::Package:: *)

BeginPackage["Lebedev`"]

GroupInfoGenerator::usage =
        "GroupInfoGenerator[templateFilename_, spaceGroupName_, spaceGroupGenerators_, baseAsymUnitExpr_, asymUnitID_, sectionTemplateElementStr_, sectionTemplateIndexXYZStr_, sectionTemplateTotalElements_] 
        generates the C++ code of transform data between ModeO and ModeS.
        ModeO indicates that the data is compressed in the orientational parameter, 
        while ModeS indecates that the data is compressed in the spatial parameter"

Begin["`Private`"]


MoveOps[op_]:=Composition[TranslationTransform[-Floor[op[{0,0,0}]]],op];


ExpandOpsSingle[generator_]:=
DeleteDuplicates[
Join[
{AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}]},
Table[
MoveOps[Composition@@Table[generator,i]]
,{i,1,5}]
]
];


OpsProduct[a_,b_]:=Flatten[Table[Composition[i,j],{i,a},{j,b}],1];


ExpandOps[generators_]:=Module[
{exops},
exops=ExpandOpsSingle/@generators;
MoveOps/@Fold[Composition[OpsProduct[#2,#1]]&,{AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}]},exops]
];


RemoveTranslationalPart[t_]:=Composition[TranslationTransform[-t[{0,0,0}]],t];


MakeUPointList[u_,uops_]:=DeleteDuplicates[(#@u)&/@uops];


GatherUPoint[u_,sgops_,uops_]:=Module[
{sgrotops,allupos,isequiv,ans},
sgrotops=DeleteDuplicates[RemoveTranslationalPart/@sgops];
allupos=MakeUPointList[u,uops];
allupos=Table[{i,allupos[[i]]},{i,1,Length[allupos]}];
isequiv[a_,b_]:=AnyTrue[sgrotops,(SameQ[#[a],b])&];
ans=Gather[allupos,isequiv[#1[[2]],#2[[2]]]&];
SortBy[#,First]&/@ans
];


TagList[list_]:=Transpose[{Range[Length[list]],list}];


TagListToRule[list_]:=(#[[2]]->#[[1]])&/@list;


MakeMappingRelations[lp_,sgops_]:=Module[
{allrelatedpoints,gathered},
allrelatedpoints=RemoveTranslationalPart[#][lp]&/@sgops;
gathered=SortBy[#,First]&/@GatherBy[TagList[allrelatedpoints],#[[2]]&];
{#[[1,2]],#[[All,1]]}&/@gathered
];


MakeSimplifiedMatchingIndex[utemp_,sgops_,uops_]:=Module[
{gatheredu,leadingu,rule},
gatheredu=GatherUPoint[utemp,sgops,uops];
leadingu=gatheredu[[All,1]];
rule=TagListToRule[TagList[MakeUPointList[utemp,uops]]];
SortBy[#,First]&/@((MakeMappingRelations[#[[2]],sgops]&/@leadingu)/.rule)
];
(*
\:51fd\:6570MakeSimplifiedMatchingIndex\:5b9e\:9645\:4e0a\:7ed9\:51fa\:4e86\:ff1a
Subscript[\:5728\:4efb\:610f\:7b26\:5408\:70b9\:7fa4G, u]\:ff08\:5143\:7d20\:5bf9\:5e94\:7684\:53d8\:6362\:4e3a\:5217\:8868uops\:ff09\:7684\:7403\:9762\:91c7\:6837\:70b9\:4e0a\:7684\:4f20\:64ad\:5b50\:ff0c
Subscript[\:5728\:7b26\:5408\:7a7a\:95f4\:7fa4G, s]\:ff08\:5143\:7d20\:5bf9\:5e94\:7684\:53d8\:6362\:4e3a\:5217\:8868sgops\:ff09\:7684\:5bf9\:79f0\:6027\:ff08\:6216\:8054\:5408\:5bf9\:79f0\:6027\:ff09\:7684\:5916\:573a\:4e0b\:ff0c
1.\:4ece \:5728\:53d6\:5411\:53d8\:91cf\:538b\:7f29\:7684\:5b58\:50a8\:65b9\:5f0f \:5230 \:5728\:7a7a\:95f4\:53d8\:91cf\:538b\:7f29\:7684\:5b58\:50a8\:65b9\:5f0f \:7684\:5bf9\:5e94\:5173\:7cfb
\:5bf9\:7ed3\:679c\:4e2d\:6bcf\:4e00\:4e2a\:5217\:8868\:ff0c\:5217\:8868\:5143\:7d20\:7684\:7b2c\:4e00\:4e2a\:6570\:503c\:4e3a\:7403\:9762\:91c7\:6837\:70b9\:7684\:7f16\:53f7\:ff0c\:5176\:5bf9\:5e94\:7684\:4f20\:64ad\:5b50\:90fd\:901a\:8fc7\:5217\:8868\:5143\:7d20\:7684\:7b2c\:4e8c\:4e2a\:503c\:ff08\:4e5f\:662f\:4e2a\:5217\:8868\:ff09\:4e2d\:7684\:4efb\:610f\:4e00\:4e2a\:5bf9\:79f0\:6027\:64cd\:4f5c\:5bf9\:5e94\:5230\:8be5\:5217\:8868\:7684\:7b2c\:4e00\:4e2a\:91c7\:6837\:70b9\:4e0a\:7684\:6570\:636e\:ff0c\:5373\:5b9e\:9645\:4f7f\:7528\:65f6\:4f7f\:7528\:4efb\:610f\:4e00\:4e2a\:5373\:53ef\:83b7\:5f97\:76ee\:6807\:4e0d\:5bf9\:79f0\:5355\:5143\:4e2d\:7684\:6570\:636e\:3002
2.\:4ece \:5728\:7a7a\:95f4\:53d8\:91cf\:538b\:7f29\:7684\:5b58\:50a8\:65b9\:5f0f \:5230 \:5728\:53d6\:5411\:53d8\:91cf\:538b\:7f29\:7684\:5b58\:50a8\:65b9\:5f0f \:7684\:5bf9\:5e94\:5173\:7cfb
\:5bf9\:5e94\:5173\:7cfb\:4e0e1.\:4e2d\:7c7b\:4f3c\:ff0c\:4f46\:6ce8\:610f\:ff0c\:7531\:4e8e\:9700\:8981\:6062\:590d\:5b8c\:6574\:7684\:7a7a\:95f4\:53d8\:91cf\:4e0a\:7684\:6570\:636e\:ff0c\:5217\:8868\:5143\:7d20\:4e2d\:7684\:7b2c\:4e8c\:4e2a\:503c\:6240\:5bf9\:5e94\:7684\:6bcf\:4e2a\:5bf9\:79f0\:6027\:64cd\:4f5c\:90fd\:9700\:8981\:8c03\:7528\:ff0c\:4ee5\:83b7\:5f97\:5b8c\:6574\:7684\:5b9e\:7a7a\:95f4\:6570\:636e\:3002
*)


GetBoundingBox[xyzexp_]:={Minimize[{#,xyzexp},{x,y,z}][[1]],Maximize[{#,xyzexp},{x,y,z}][[1]]}&/@{x,y,z};


ApplyTransformToExpr[xyzexp_,opt_]:=Simplify[xyzexp/.((#[[1]]->#[[2]])&/@Transpose[{{x,y,z},opt[{x,y,z}]}])];


GetSamplePointsMiddle[range_]:=Module[
{lower,upper},
lower=Floor[range[[1]]];
upper=Ceiling[range[[2]]];
Select[Range[lower,upper,1/2],Mod[#,1]!=0&&#>=range[[1]]&&#<= range[[2]]&]
];


GetSampleShiftAndLength[range_,sizelist_]:=Module[
{posrange,samplesxyz},
posrange=range*sizelist;
samplesxyz=GetSamplePointsMiddle/@posrange;
{#[[1]]-1/2,Length[#]}&/@samplesxyz
];


MakeSizeSymbol[list_]:=StringRiffle[list/.{x->"X",y->"Y",z->"Z"},{"N","",""}];


DetectSizeType[ops_]:=Module[
{rots,exps},
rots=DeleteDuplicates[RemoveTranslationalPart/@ops];
exps=#[{x,y,z}]&/@rots;
exps=DeleteDuplicates[exps/.{-x->x,-y->y,-z->z}];
MakeSizeSymbol/@Sort/@DeleteDuplicates/@Transpose[exps]
];


FormatMatrix[mat_]:=StringRiffle[StringRiffle[#,", "]&/@Partition[(If[#<100,If[#<10,"  "," "],""]<>ToString[#])&/@mat,UpTo[12]],{"{{\n    ",", \n    ","\n}};"}];


FormatMatrixSingleCB[mat_]:=StringRiffle[StringRiffle[#,", "]&/@Partition[(If[#<100,If[#<10,"  "," "],""]<>ToString[#])&/@mat,UpTo[12]],{"{\n    ",", \n    ","\n};"}];


FormatMatrixInline[mat_]:=StringRiffle[StringRiffle[#,", "]&/@{(ToString[#])&/@mat},{"{{",", ","}}"}];


MakeExtractPosList[match_]:=Module[
{temp,ans},
temp=match[[All,All,1]];
ans=SortBy[Flatten[Table[{match[[i,j]],i},{i,1,Length[temp]},{j,1,Length[temp[[i]]]}],1],First][[All,2]];
FormatMatrixSingleCB[ans-1]
];


GetGroupRankEx[match_]:=Length[match]*Length[Flatten[match[[1,All,2]],1]];


MakeReconstructLebedevPointPosList[match_]:=Module[
{temp,ans},
temp=Sort@Flatten[Table[
{i,k,match[[i,j,1]]}
,{i,1,Length[match]},{j,1,Length[match[[i]]]},{k,match[[i,j,2]]}],2];
ans=temp[[All,3]];
FormatMatrixSingleCB[ans-1]
];


MakeReconstructFieldPosList[match_]:=Module[
{temp,ans},
temp=Sort@Flatten[Table[
{i,k,match[[i,j,1]]}
,{i,1,Length[match]},{j,1,Length[match[[i]]]},{k,match[[i,j,2]]}],2];
ans=temp[[All,1]];
FormatMatrixSingleCB[ans-1]
];


MoveOpsAccordingToBBox[ops_,bbox_]:=Composition[TranslationTransform[-
(
Min/@Transpose[{Floor[#[bbox[[All,1]]]],Floor[#[bbox[[All,2]]]]}]
)
],#]&/@ops;


MakeOpExp[op_,sizelist_]:=Expand@Simplify[op[({x,y,z}+{1/2,1/2,1/2})/sizelist]*sizelist-{1/2,1/2,1/2}];


GenerateExtractStr[opexps_]:=StringRiffle[#,"\n"]&@
Table[
"case "<>ToString[i-1]<>":GROUP__GROUPID___INDEX_CASES("<>ToString[i-1]<>", ("<>ToString[InputForm[opexps[[i,1]]]]<>"), ("<>ToString[InputForm[opexps[[i,2]]]]<>"), ("<>ToString[InputForm[opexps[[i,3]]]]<>"));",
{i,1,Length[opexps]}
];


DetectNewVars[exps_,sizelist_,bbox_]:=Module[
{allexps,sizelistlite,maxden,newvar,ans},
allexps=Expand/@Flatten[Join[exps,{sizelist*bbox[[All,1]],sizelist*bbox[[All,2]]}]];
sizelistlite=sizelist;
maxden=Table[
Max@Abs@Denominator@(Coefficient[#,size]&/@allexps)
,{size,sizelistlite}];
newvar=ToExpression[ToString[#[[1]]]<>ToString[#[[2]]]]&/@Transpose[{sizelistlite,maxden}];
ans=sizelist/.(#[[1]]->#[[2]]&/@Transpose[{sizelistlite,maxden*newvar}]);
{ans,newvar,maxden}
];


MakeSectionTemplateDerivedSizes[oldsizes_,newsizesinfo_]:=StringRiffle[
Join[
("constexpr SignedT "<>ToString[#1]<>" = "<>ToString[#2]<>" / "<>ToString[#3]<>";")&@@@DeleteDuplicates[Transpose[{newsizesinfo[[2]],oldsizes,newsizesinfo[[3]]}]],
("static_assert("<>ToString[#2]<>" == "<>ToString[#1]<>" * "<>ToString[#3]<>", \"template parameter "<>ToString[#2]<>" must be divided by "<>ToString[#3]<>"\");")&@@@DeleteDuplicates[Transpose[{newsizesinfo[[2]],oldsizes,newsizesinfo[[3]]}]]
],
"\n"
];


MakeSectionTemplateExtractSwitches[matches_,inverseops_,bbox_,sizelist_]:=Module[
{movedinversesgops,tempmatch},
movedinversesgops=MoveOpsAccordingToBBox[inverseops,bbox];
tempmatch=SortBy[Flatten[matches,1],First];
GenerateExtractStr[MakeOpExp[#,sizelist]&/@movedinversesgops[[Flatten[tempmatch[[All,2,1]]]]]]
];


MakeSectionTemplateReconstructSwitches[matches_,ops_,bbox_,sizelist_]:=Module[
{movedops,nfield},
movedops=MoveOpsAccordingToBBox[ops,bbox];
nfield=Length[matches];
GenerateExtractStr[MakeOpExp[#,sizelist]&/@Flatten[Table[movedops,{i,1,nfield}]]]
];


ModifyTarget[line_,rule_]:=Module[
{s},
s=StringReplace[line,StartOfString ~~s:Whitespace...~~rule[[1]]:>s];
StringReplace[StringTrim[rule[[2]]],StartOfLine~~a_:>s<>a]
];


SectionTemplateReplateFunc[str_,rules_]:=Module[
{lines},
lines=StringSplit[str,"\n"];
Table[
lines=Table[
If[StringContainsQ[line,rule[[1]]],
ModifyTarget[line,rule]
,
line]
,{line,lines}];
,{rule,rules}];
StringRiffle[lines,"\n"]
];


TransformCRLF2LF[str_]:=StringReplace[str,"\r\n"->"\n"];


checkManullyDerivedFuncs[testSize_,elementStr_,indexStr_,totalStr_,sizesymbols_,newvarinfo_,bbox_,asymunitexpr_]:=Module[
{usedTestSize,testSizeListStr,helpersizes,rangestr,testFuncStr,testapp,calcans,fullTestsizes,refans,compare},
Needs["CCompilerDriver`"];
usedTestSize=testSize[[1;;Length[DeleteDuplicates[sizesymbols]]]];
testSizeListStr=StringRiffle[("int "<>ToString[#1]<>" = "<>ToString[#2])&@@@Transpose[{DeleteDuplicates[sizesymbols],usedTestSize}],{"",";\n ",";\n"}];
helpersizes=StringRiffle[
Join[
("constexpr SignedT "<>ToString[#1]<>" = "<>ToString[#2]<>" / "<>ToString[#3]<>";")&@@@DeleteDuplicates[Transpose[{newvarinfo[[2]],sizesymbols,newvarinfo[[3]]}]]
],
"\n"
];
rangestr=ToString[InputForm[#]]&/@Flatten[bbox*newvarinfo[[1]]];
testFuncStr=StringReplace[
"#include <stdio.h>\n"
<>"\nint element(int x, int y, int z){\n"
<>testSizeListStr
<>elementStr
<>"\n}\n"
<>"\nint index(int x, int y, int z){\n"
<>testSizeListStr
<>indexStr
<>"\n}\n"
<>"\nint total(){\n"
<>testSizeListStr
<>totalStr
<>"\n}\n"
<>"\nint main(){
    "<>testSizeListStr<>"
    "<>helpersizes<>"
	int x, y, z;
	for(x = "<>rangestr[[1]]<>"; x < ("<>rangestr[[2]]<>"); ++x){
	for(y = "<>rangestr[[3]]<>"; y < ("<>rangestr[[4]]<>"); ++y){
	for(z = "<>rangestr[[5]]<>"; z < ("<>rangestr[[6]]<>"); ++z){
		printf(\"%d \", (int)element(x,y,z)*(int)index(x,y,z));
	}
	}
	}
	printf(\"\\n%d\", (int)total());
}
",
{
"index_xyz<TEMPLATE_reduced_size_list, SignedT>"->"index",
"size_t"->"int",
"SignedT"->"int",
"constexpr"->"",
"true"->"1",
"false"->"0"
}
];
testapp=Check[CCompilerDriver`CreateExecutable[testFuncStr, "testapp"],False];
If[testapp==False,
Print["Error occurs whem compile c source:\n", testFuncStr];Abort[]
,Nothing[]];
calcans=ToExpression@StringSplit[Import["!"<>CCompilerDriver`QuoteFile[testapp],"Text"]];
fullTestsizes=sizesymbols/.Flatten[{((#1->#2)&@@@Transpose[{DeleteDuplicates[sizesymbols],usedTestSize}])}];
refans=Table[
asymunitexpr/.{x->((i+1/2)/fullTestsizes[[1]]),y->((j+1/2)/fullTestsizes[[2]]),z->((k+1/2)/fullTestsizes[[3]])},
{i,bbox[[1,1]]*fullTestsizes[[1]],bbox[[1,2]]*fullTestsizes[[1]]-1},
{j,bbox[[2,1]]*fullTestsizes[[2]],bbox[[2,2]]*fullTestsizes[[2]]-1},
{k,bbox[[3,1]]*fullTestsizes[[3]],bbox[[3,2]]*fullTestsizes[[3]]-1}];
refans=Flatten[refans]/.{True->1,False->0};
refans=Join[refans*(Accumulate[refans]-1),{Total[refans]}];
compare=SameQ[calcans,refans];
{compare,fullTestsizes,calcans,refans}
]


GroupInfoGenerator[templateFilename_, targetDirname_, spaceGroupName_, spaceGroupGenerators_, baseAsymUnitExpr_, asymUnitID_, sectionTemplateElementStr_, sectionTemplateIndexXYZStr_, sectionTemplateTotalElementsStr_]:=Module[
{
esc,
privateBaseAsymUnitExpr=baseAsymUnitExpr,
ohop0,ohop1,ohop2,ohop3,ohop4,ohop5,ohops48,
sgops, inversesgops,asymunitexpr,bbox,sizesymbolnames,sizesymbols,newvarinfo,newsizesymbols,
compareans,
matches,sectionTemplateExtractFieldPosMats,groupRandExs,sectionTemplateDerivedSizes,sectionTemplateReconstructLebedevPointPosMats,sectionTemplateReconstructFieldPosMats,sectionTemplateExtractSwitches,sectionTemplateReconstructSwitches,
template,targetname,ansstr,
deviceansstr,hostansstr
},

privateBaseAsymUnitExpr=baseAsymUnitExpr/.{Global`x->x,Global`y->y,Global`z->z};

(*operations for oh point group*)
ohop0=AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}];
ohop1=AffineTransform[{DiagonalMatrix[{-1,-1,1}],{0,0,0}}];
ohop2=AffineTransform[{DiagonalMatrix[{-1,1,-1}],{0,0,0}}];
ohop3=AffineTransform[{{{0,0,1},{1,0,0},{0,1,0}},{0,0,0}}];
ohop4=AffineTransform[{{{0,1,0},{1,0,0},{0,0,-1}},{0,0,0}}];
ohop5=AffineTransform[{DiagonalMatrix[{-1,-1,-1}],{0,0,0}}];
ohops48={ohop0};
ohops48=Flatten[{#,Function[op,Composition[ohop1,op]]/@#},1]&@ohops48;
ohops48=Flatten[{#,Function[op,Composition[ohop2,op]]/@#},1]&@ohops48;
ohops48=Flatten[{#,Function[op,Composition[ohop3,op]]/@#,Function[op,Composition[ohop3,ohop3,op]]/@#},1]&@ohops48;
ohops48=Flatten[{#,Function[op,Composition[ohop4,op]]/@#},1]&@ohops48;
ohops48=Flatten[{#,Function[op,Composition[ohop5,op]]/@#},1]&@ohops48;

(*operations for target space group*)
sgops=ExpandOps[spaceGroupGenerators];
inversesgops=MoveOps/@InverseFunction/@sgops;

(*inspect the asymmetric unit*)
asymunitexpr=ApplyTransformToExpr[privateBaseAsymUnitExpr,sgops[[asymUnitID]]];
bbox=GetBoundingBox[asymunitexpr];


(*deduce suitiable unit cell vars*)
sizesymbolnames=DetectSizeType[sgops];
sizesymbols=ToExpression/@sizesymbolnames;
newvarinfo=DetectNewVars[MakeOpExp[#,sizesymbols]&/@sgops,sizesymbols,bbox];
newsizesymbols=newvarinfo[[1]];

(*print informations to let the user validate the inputs*)
esc = Association[
"reset" -> "\033[1;0m",
"black" -> "\033[1;30m", "red" -> "\033[1;31m",
"green" -> "\033[1;32m", "yellow" -> "\033[1;33m",
"blue" -> "\033[1;34m", "magenta" -> "\033[1;35m"
];

Print["============================================================================================"];
Print["generating the C++ file for space group ",esc["red"],spaceGroupName,esc["reset"]];
Print["the selected asymmetric unit is ",esc["red"],StringReplace[ToString[InputForm[asymunitexpr]],"Lebedev`Private`"->""],esc["reset"]];
Print["data in the asymmetric unit fall in the range of ",esc["red"],bbox,esc["reset"]];
Print["the sizes of the computational box is labeled as ",esc["red"],sizesymbolnames,esc["reset"]];
Print["****************************************************************************************"];
Print[esc["magenta"], "please VALIDATE the following MANUALLY DERIVED parts according to these settings", esc["reset"]];
Print["****************************************************************************************"];
Print["The sectionTemplateElementStr is set to be:\n",esc["red"],sectionTemplateElementStr,esc["reset"]];
Print["The sectionTemplateIndexXYZStr is set to be:\n",esc["red"],sectionTemplateIndexXYZStr,esc["reset"]];
Print["The sectionTemplateTotalElementsStr is set to be:\n",esc["red"],sectionTemplateTotalElementsStr,esc["reset"]];

compareans=checkManullyDerivedFuncs[{32,48,64},sectionTemplateElementStr,sectionTemplateIndexXYZStr,sectionTemplateTotalElementsStr,sizesymbols,newvarinfo,bbox,asymunitexpr];

If[compareans[[1]],
Print[esc["green"], "Message: These MANUALLY DERIVED parts have passed the validation with the size of "<>ToString[compareans[[2]]]<>"\n", esc["reset"]],
Print[
esc["magenta"], "error: These MANUALLY DERIVED parts have failed the validation with the size of "<>ToString[compareans[[2]]]<>"\n       Although the file will be generated anyway, but it is probably incorrect!\n", esc["reset"],
"The calculated results are ",ToString[compareans[[3]]]<>"\n",
"while the reference results are",ToString[compareans[[4]]]<>"\n"
]
];

matches=MakeSimplifiedMatchingIndex[#,sgops,ohops48]&/@
{
{0,0,c},{0,b,b},{a,a,a},{a,a,c},{a,b,0},{a,b,c}
};
sectionTemplateExtractFieldPosMats=MakeExtractPosList/@matches;
groupRandExs=GetGroupRankEx/@matches;
sectionTemplateDerivedSizes=MakeSectionTemplateDerivedSizes[sizesymbols,newvarinfo];
sectionTemplateReconstructLebedevPointPosMats=MakeReconstructLebedevPointPosList/@matches;
sectionTemplateReconstructFieldPosMats=MakeReconstructFieldPosList/@matches;
sectionTemplateExtractSwitches=MakeSectionTemplateExtractSwitches[#,inversesgops,bbox,newsizesymbols]&/@matches;
sectionTemplateReconstructSwitches=MakeSectionTemplateReconstructSwitches[#,inversesgops,bbox,newsizesymbols]&/@matches;

template=Import[templateFilename,"Text"];
targetname=FileNameJoin[{targetDirname,StringReplace[FileNameTake[templateFilename],{"XXX"->spaceGroupName,".template"->""}]}];
ansstr=SectionTemplateReplateFunc[template,
Join[
{
"SECTION_TEMPLATE_element"->TransformCRLF2LF@sectionTemplateElementStr,
"SECTION_TEMPLATE_index_xyz"->TransformCRLF2LF@sectionTemplateIndexXYZStr,
"SECTION_TEMPLATE_total_elements"->TransformCRLF2LF@sectionTemplateTotalElementsStr,
"SECTION_TEMPLATE_derived_sizes"->sectionTemplateDerivedSizes
},
Table[("SECTION_TEMPLATE_extract_field_pos_mat"<>ToString[i])->sectionTemplateExtractFieldPosMats[[i+1]],{i,0,5}],
Table[("SECTION_TEMPLATE_reconstruct_lebedev_point_pos_mat"<>ToString[i])->sectionTemplateReconstructLebedevPointPosMats[[i+1]],{i,0,5}],
Table[("SECTION_TEMPLATE_reconstruct_field_pos_mat"<>ToString[i])->sectionTemplateReconstructFieldPosMats[[i+1]],{i,0,5}],
Table[("SECTION_TEMPLATE_extract_switch"<>ToString[i])->sectionTemplateExtractSwitches[[i+1]],{i,0,5}],
Table[("SECTION_TEMPLATE_reconstruct_switch"<>ToString[i])->sectionTemplateReconstructSwitches[[i+1]],{i,0,5}]
]
];
ansstr=StringReplace[ansstr,
Join[
{
"TEMPLATE_operation_number"->ToString[Length[sgops]],
"TEMPLATE_fields_required_table"->FormatMatrixInline[Length/@matches],
"TEMPLATE_lengthx"->ToString[InputForm[(bbox[[1,2]]-bbox[[1,1]])*NX]],
"TEMPLATE_lengthy"->ToString[InputForm[(bbox[[2,2]]-bbox[[2,1]])*NY]],
"TEMPLATE_lengthz"->ToString[InputForm[(bbox[[3,2]]-bbox[[3,1]])*NZ]],
"TEMPLATE_asymunit_expr"->ToString[InputForm[asymunitexpr]] ,
"TEMPLATE_size_list"->StringRiffle[("size_t "<>ToString[#])&/@DeleteDuplicates[sizesymbols],", "],
"TEMPLATE_reduced_size_list"->StringRiffle[ToString/@DeleteDuplicates[sizesymbols],", "],
"TEMPLATE_asymunit_indexrangex"->"["<>ToString[InputForm[bbox[[1,1]]]]<>", "<>ToString[InputForm[bbox[[1,2]]]]<>"]",
"TEMPLATE_asymunit_indexrangey"->"["<>ToString[InputForm[bbox[[2,1]]]]<>", "<>ToString[InputForm[bbox[[2,2]]]]<>"]",
"TEMPLATE_asymunit_indexrangez"->"["<>ToString[InputForm[bbox[[3,1]]]]<>", "<>ToString[InputForm[bbox[[3,2]]]]<>"]",
"TEMPLATE_Nx"->ToString[sizesymbols[[1]]],
"TEMPLATE_Ny"->ToString[sizesymbols[[2]]],
"TEMPLATE_Nz"->ToString[sizesymbols[[3]]],
"TEMPLATE_Sx"->ToString[InputForm[Expand[bbox[[1,1]]*newsizesymbols[[1]]]]],
"TEMPLATE_Sy"->ToString[InputForm[Expand[bbox[[2,1]]*newsizesymbols[[2]]]]],
"TEMPLATE_Sz"->ToString[InputForm[Expand[bbox[[3,1]]*newsizesymbols[[3]]]]]
},
Table[("TEMPLATE_field_indexes"<>ToString[i])->FormatMatrixInline[matches[[i+1,All,1,1]]-1],{i,0,5}],
Table[("TEMPLATE_field_submultiplicities"<>ToString[i])->FormatMatrixInline[Length/@matches[[i+1]]],{i,0,5}]
]

];
ansstr=StringReplace[ansstr,
Join[
{
"__GROUPID__"->spaceGroupName
},
Table[("__GROUPRANKEX_"<>ToString[i]<>"__")->ToString[groupRandExs[[i+1]]],{i,0,5}]
]

];
ansstr=StringReplace[ansstr,"Lebedev`Private`"->""];

deviceansstr=SectionTemplateReplateFunc[ansstr,
{
"SECTION_TEMPLATE_EXTRACT_ITERATER_OUTER_LOOP_START"->"{",
"SECTION_TEMPLATE_EXTRACT_ITERATER_OUTER_LOOP_END"->"}",
"SECTION_TEMPLATE_RECONSTRUCT_ITERATER_OUTER_LOOP_START"->"{",
"SECTION_TEMPLATE_RECONSTRUCT_ITERATER_OUTER_LOOP_END"->"}",
"SECTION_TEMPLATE_ADDITIONAL_INCLUDE"->
"#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include \"device_helper.h\""
}
];

deviceansstr=StringReplace[deviceansstr,
{
"TEMPLATE_ITERATER_OUTER"->"blockIdx.x",
"TEMPLATE_ITERATER_X_START"->"(SignedT)(blockIdx.y)",
"TEMPLATE_ITERATER_X_STEP"->"(SignedT)(gridDim.y)",
"TEMPLATE_ITERATER_Y_START"->"(SignedT)(threadIdx.x)",
"TEMPLATE_ITERATER_Y_STEP"->"(SignedT)(blockDim.x)",
"TEMPLATE_ITERATER_Z_START"->"(SignedT)(threadIdx.y)",
"TEMPLATE_ITERATER_Z_STEP"->"(SignedT)(blockDim.y)",
"TEMPLATE_GLOBAL"->"LEBEDEV_GROUP_INFO_GLOBAL",
"TEMPLATE_BOTH_CALLABLE"->"LEBEDEV_GROUP_INFO_BOTH_CALLABLE",
"TEMPLATE_BOTH_INLINE"->"LEBEDEV_GROUP_INFO_BOTH_INLINE",
"TEMPLATE_TARGET"->"device",
"TEMPLATE_DIM3_ELE_TYPE"->"decltype(dim3::x)"
}
];

hostansstr=SectionTemplateReplateFunc[ansstr,
{
"SECTION_TEMPLATE_EXTRACT_ITERATER_OUTER_LOOP_START"->
"for(SignedT outer = 0; outer < extract_num<PType, SignedT>; ++outer)
{",
"SECTION_TEMPLATE_EXTRACT_ITERATER_OUTER_LOOP_END"->"}",
"SECTION_TEMPLATE_RECONSTRUCT_ITERATER_OUTER_LOOP_START"->
"for(SignedT outer = 0; outer < reconstruct_num<PType, SignedT>; ++outer)
{",
"SECTION_TEMPLATE_RECONSTRUCT_ITERATER_OUTER_LOOP_END"->"}",
"SECTION_TEMPLATE_ADDITIONAL_INCLUDE"->""
}
];

hostansstr=StringReplace[hostansstr,
{
"TEMPLATE_ITERATER_OUTER"->"outer",
"TEMPLATE_ITERATER_X_START"->"0",
"TEMPLATE_ITERATER_X_STEP"->"1",
"TEMPLATE_ITERATER_Y_START"->"0",
"TEMPLATE_ITERATER_Y_STEP"->"1",
"TEMPLATE_ITERATER_Z_START"->"0",
"TEMPLATE_ITERATER_Z_STEP"->"1",
"TEMPLATE_GLOBAL "->"",(*Note that there's a additional space here*)
"TEMPLATE_BOTH_CALLABLE "->"",(*Note that there's a additional space here*)
"TEMPLATE_BOTH_INLINE"->"inline",
"TEMPLATE_TARGET"->"host",
"TEMPLATE_DIM3_ELE_TYPE"->"int"
}
];

Export[targetname,deviceansstr,"Text",CharacterEncoding->"UTF8"];
Export[StringReplace[targetname,".cuh"~~EndOfString->".h"],hostansstr,"Text",CharacterEncoding->"UTF8"];
Print["the output C++ file is ",esc["red"],targetname,esc["reset"], " and ",esc["red"],StringReplace[targetname,".cuh"~~EndOfString->".h"],esc["reset"]];
Print["============================================================================================"];
]
End[ ]

EndPackage[ ]
