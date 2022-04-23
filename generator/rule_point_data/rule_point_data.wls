#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Script assumptions :
    1.  Script has not changed directories, as with SetDirectory[]
    2.  Mathematica 11.0, or later
    3.  This file is only launched from console, not from interactive input mode
*)
lebedevSrcFilename=($ScriptCommandLine)[[2]];
templateFilename=($ScriptCommandLine)[[3]];
template=Import[templateFilename];
outputFilename=($ScriptCommandLine)[[4]];

scriptPath = ExpandFileName[First[$ScriptCommandLine]]; 
scriptName = Last[FileNameSplit[scriptPath]]; 

Print["Script ", scriptPath, " is using source file ", lebedevSrcFilename," and template file ", templateFilename, " to generate ", outputFilename];

ModifyTarget[line_,rule_]:=Module[
{s},
s=StringReplace[line,StartOfString ~~s:Whitespace..~~rule[[1]]:>s];
StringReplace[
    StringReplace[StringTrim[rule[[2]]],StartOfLine~~a_:>s<>a],
    StartOfLine~~Whitespace..~~EndOfLine->""
]
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

GetLebedevFileData[filename_]:=Module[
{filedata},
filedata=Import[filename,"Text"];
filedata=StringTrim/@StringSplit[filedata,"\n"];
filedata=Select[filedata,StringLength[#]>0&&StringTake[#,1]!="/"&];
filedata
];
SliceFileData[filedata_]:=Module[
{poslist,names,lastpos,fileslices},
poslist=Position[StringMatchQ[filedata,"void ld"~~_~~_~~_~~_~~"(double* x, double* y, double* z, double* w)"],True][[All,1]];
names=ToExpression@StringReplace[filedata[[poslist]],"void ld"~~a:_~~b:_~~c:_~~d:_~~"(double* x, double* y, double* z, double* w)":>a<>b<>c<>d];
lastpos=Position[filedata,"int order_table(int rule)"][[1,1]];
AppendTo[poslist,lastpos];
fileslices=Table[
filedata[[poslist[[i]]+7;;poslist[[i+1]]-4]],{i,1,Length[poslist]-1}];
{names,fileslices}
];
GetPointData[type_,str_]:=Module[
{interp,a,b,c,x,y,z,v,values,symbols},
(*Note this part cannot be wirtten as Interpreter["Number"][...], since this will results in a delayed estimation*)
(*Probably a bug in current WolframEngine?*)
interp=Interpreter["Number"];
values=interp[StringReplace[#,val:_~~" = "~~num:__~~";":>num]]&/@str;
symbols=StringTake[#,1]&/@str;
Table[
Switch[symbols[[i]],
"a",a=values[[i]],
"b",b=values[[i]],
"v",v=values[[i]]

]
,{i,1,Length[symbols]}];
Switch[type,
1, 
x=0.;y=0.;z=1.;,
2,
x=0.;y=Sqrt[2.]/2.;z=Sqrt[2.]/2.;,
3,
x=Sqrt[3.]/3.;y=Sqrt[3.]/3.;z=Sqrt[3.]/3.;,
4,
b=Sqrt[1-2a^2];x=a;y=a;z=b;,
5,
b=Sqrt[1-a^2];x=a;y=b;z=0.;,
6,
c=Sqrt[1-a^2-b^2];x=a;y=b;z=c;
];
{type,x,y,z,v}

];
AnalyzeSlices[src_]:=Module[
{subposlist,typelist,pointdata},
subposlist=Position[StringMatchQ[src,"n = n + gen_oh("~~_~~", a, b, v, x + n, y + n, z + n, w + n);"],True][[All,1]];
typelist=ToExpression/@StringReplace[src[[subposlist]],"n = n + gen_oh("~~a:_~~", a, b, v, x + n, y + n, z + n, w + n);":>a];
PrependTo[subposlist,0];
pointdata=Table[src[[subposlist[[i]]+1;;subposlist[[i+1]]-1]],{i,1,Length[subposlist]-1}];
Table[GetPointData[typelist[[i]],pointdata[[i]]],{i,1,Length[typelist]}]
];
MyNumberForm[num_]:=ToString[PaddedForm[num,24,ScientificNotationThreshold->{-10,9}]];
FormatLine[list_]:=StringRiffle[MyNumberForm/@list,{"{{",",","}}"}];
FormatAns[ans_]:=StringRiffle[FormatLine/@ans,{"{{\n    ",",\n    ","\n}}"}];
PrettyFormatAns[ruleindex_,ans_]:=Module[
{header,footer,type,pos,weight},
header=
"
template<>
struct RulePointData<"<>ToString[ruleindex]<>">{
    static constexpr size_t Rule_="<>ToString[ruleindex]<>";
    static constexpr size_t point_type_total_=point_type_accum_count(Rule_)[5];";
footer="
};";
type=StringRiffle[ans[[All,1]]/.{1->"LEBEDEV_POINT_TYPE::OPTRN6_00C",2->"LEBEDEV_POINT_TYPE::OPTRN12_0BB",3->"LEBEDEV_POINT_TYPE::OPTRN8_AAA",4->"LEBEDEV_POINT_TYPE::OPTRN24_AAC",5->"LEBEDEV_POINT_TYPE::OPTRN24_AB0",6->"LEBEDEV_POINT_TYPE::OPTRN48_ABC"},{
"
    static constexpr c_array<LEBEDEV_POINT_TYPE, point_type_total_> type_list_={{
        ",
",
        ",
"
    }};"

}];
pos=StringRiffle[FormatLine/@ans[[All,{2,3,4}]],{
"
    static constexpr c_array<c_array<double, 3>, point_type_total_> pos_list_={{
      ",
",
      ",
"
    }};"
}];
weight=StringRiffle[(StringTrim@*MyNumberForm)/@ans[[All,5]],{
"
    static constexpr c_array<double, point_type_total_> weight_list_={{
        ",
",
        ",
"
    }};"
}];
header<>type<>pos<>weight<>footer
];
filedata=GetLebedevFileData[lebedevSrcFilename];
{names,fileslices}=SliceFileData[filedata];
anss=AnalyzeSlices[#]&/@fileslices;
pointCountTable={6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,386,434,482,530,590,650,698,770,830,890,974,1046,1118,1202,1274,1358,1454,1538,1622,1730,1814,1910,2030,2126,2222,2354,2450,2558,2702,2810,2930,3074,3182,3314,3470,3590,3722,3890,4010,4154,4334,4466,4610,4802,4934,5090,5294,5438,5606,5810};
poss=Position[pointCountTable,#][[1,1]]&/@names;
ans=StringRiffle[Table[PrettyFormatAns[poss[[i]],anss[[i]]],{i,1,Length[anss]}],{"\n","\n","\n"}];

ans=SectionTemplateReplateFunc[template, {"SECTION_TEMPLATE_POINT_DATA"->ans}];

Export[
outputFilename,
ans,
"Text"
];
