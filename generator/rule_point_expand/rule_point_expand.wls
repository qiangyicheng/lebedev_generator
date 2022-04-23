#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Script assumptions :
    1.  Script has not changed directories, as with SetDirectory[]
    2.  Mathematica 11.0, or later
    3.  This file is only launched from console, not from interactive input mode
*)
templateFilename=($ScriptCommandLine)[[2]];
template=Import[templateFilename];
outputFilename=($ScriptCommandLine)[[3]];

scriptPath = ExpandFileName[First[$ScriptCommandLine]]; 
scriptName = Last[FileNameSplit[scriptPath]]; 

Print["Script ", scriptPath, " is using template file ", templateFilename, " to generate ", outputFilename];

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

genohrule[type_]:=Module[
{typepos,typename,typecount,uniquepos,prefix,surfix,
ohop0,ohop1,ohop2,ohop3,ohop4,ohop5,ohops48},
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
typepos=type/.{
0->{0,0,c},
1->{0,b,b},
2->{a,a,a},
3->{a,a,c},
4->{a,b,0},
5->{a,b,c}
};
typename=type/.{0->"OPTRN6_00C",1->"OPTRN12_0BB",2->"OPTRN8_AAA",3->"OPTRN24_AAC",4->"OPTRN24_AB0",5->"OPTRN48_ABC"};
typecount=type/.{0->6,1->12,2->8,3->24,4->24,5->48};
typecount=ToString[typecount];
uniquepos=DeleteDuplicates[(#@typepos)&/@ohops48];
prefix="
template <>
inline constexpr c_array<c_array<double, 3>, "<>typecount<>"> point_expand<LEBEDEV_POINT_TYPE::"<>typename<>">(c_array<double, 3> const &p)
{
    c_array<c_array<double, 3>, "<>typecount<>"> ans{};";
surfix="
    return ans;
}";
prefix<>
StringRiffle[
Table[
"
    ans["<>ToString[i-1]<>"][0] = "<>StringReplace[ToString[uniquepos[[i,1]]],{"a"->"p[0]","b"->"p[1]","c"->"p[2]"}]<>";
    ans["<>ToString[i-1]<>"][1] = "<>StringReplace[ToString[uniquepos[[i,2]]],{"a"->"p[0]","b"->"p[1]","c"->"p[2]"}]<>";
    ans["<>ToString[i-1]<>"][2] = "<>StringReplace[ToString[uniquepos[[i,3]]],{"a"->"p[0]","b"->"p[1]","c"->"p[2]"}]<>";"
,{i,1,Length[uniquepos]}]
,{"\n","\n","\n"}]
<>
surfix
];
ans=StringRiffle[genohrule/@{0,1,2,3,4,5},"\n\n"];

ans=SectionTemplateReplateFunc[template, {"SECTION_TEMPLATE_POINT_EXPAND"->ans}];

Export[
outputFilename,
ans,
"Text"
];
