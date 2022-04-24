#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Script assumptions :
    1.  Script has not changed directories, as with SetDirectory[]
    2.  Mathematica 11.0, or later
    3.  This file is only launched from console, not from interactive input mode
*)

GroupInfo::para="Not sufficient parameters. At least `1` required while `2` is received. Received parameters are `3`";
If[Length[$ScriptCommandLine]<4,
Message[GroupInfo::para,4,Length[$ScriptCommandLine],$ScriptCommandLine];Exit[10],
Nothing[]
]

templateFilename=($ScriptCommandLine)[[2]];
outputFoldername=($ScriptCommandLine)[[3]];
requiredGroupList=Flatten[StringSplit[#,{";",",",Whitespace}]&/@(($ScriptCommandLine)[[4;;]])];

scriptPath = ExpandFileName[First[$ScriptCommandLine]]; 
scriptName = Last[FileNameSplit[scriptPath]]; 

Print["Script ", scriptPath, " is using template file ", templateFilename, " to generate headers for group ", requiredGroupList, " in folder ", outputFoldername];

packagePath=(DirectoryName[scriptPath] <> "GroupInfoGenerator.wl");
Get[OpenRead[packagePath]];


settings=Association[];
AddSettings[sgname_,sgops_,baseasymunitexpr_,selectasymunitid_,sectionTemplateElement_, sectionTemplateIndexXYZ_, sectionTemplateTotalElements_]:=
AppendTo[settings,sgname->{sgname,sgops,baseasymunitexpr,selectasymunitid,sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements}];


(* ::Text:: *)
(* *)


(*space group 230, namely Gyroid*)
Clear[sgname,sgops,baseasymunitexpr,selectasymunitid,sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];
sgname="3D230";
sgops={
AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,1}],{1/2,0,1/2}}],
AffineTransform[{DiagonalMatrix[{-1,1,-1}],{0,1/2,1/2}}],
AffineTransform[{{{0,0,1},{1,0,0},{0,1,0}},{0,0,0}}],
AffineTransform[{{{0,1,0},{1,0,0},{0,0,-1}},{3/4,1/4,1/4}}],
AffineTransform[{DiagonalMatrix[{-1,-1,-1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{1,1,1}],{1/2,1/2,1/2}}]
}
baseasymunitexpr=(-1/8<=x<=1/8)&&(-1/8<=y<=1/8)&&(0<=z<=1/4)&&(Max[x,-x,y,-y]<=z);
(*select suitable asymmetric unit*)
selectasymunitid=9;
sectionTemplateElement=
"
constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
return x >= NXYZ8 ? true :
    (y >= -x - 1 && y <= x && z >= -x - 1 && z <= x);
";
sectionTemplateIndexXYZ=
"
constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
constexpr SignedT half_count = (2 * NXYZ8 * (1 + NXYZ8) * (2 * NXYZ8 + 1)) / 3;
return (x < NXYZ8) ?
    (2 * x * (1 + x) * (2 * x + 1)) / 3 + 2 * (y + x + 1) * (x + 1) + (z + x + 1) :
    half_count + (x - NXYZ8) * NXYZ4 * NXYZ4 + (y + NXYZ8) * NXYZ4 + (z + NXYZ8);
";
sectionTemplateTotalElements=
"
constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return index_xyz<TEMPLATE_reduced_size_list, SignedT>(NXYZ4 - 1, NXYZ8 - 1, NXYZ8 - 1 + 1/*one-past-last*/);
";
AddSettings[sgname, sgops, baseasymunitexpr, selectasymunitid, sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];


(*space group 229, namely BCC*)
Clear[sgname,sgops,baseasymunitexpr,selectasymunitid,sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];
sgname="3D229";
sgops={
AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,1,-1}],{0,0,0}}],
AffineTransform[{{{0,0,1},{1,0,0},{0,1,0}},{0,0,0}}],
AffineTransform[{{{0,1,0},{1,0,0},{0,0,-1}},{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,-1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{1,1,1}],{1/2,1/2,1/2}}]
}
baseasymunitexpr=(0<=x<=1/2)&&(0<=y<=1/2)&&(0<=z<=1/4)&&(y<=x)&&(z<=Min[1/2-x,y]);
(*select suitable asymmetric unit*)
selectasymunitid=9;
sectionTemplateElement=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return z <= y && x <= (2 * NXYZ4 - y - 1) && x <= z;
";
sectionTemplateIndexXYZ=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return (x * (5 + 18 * NXYZ4 + 12 * NXYZ4 * NXYZ4 - 9 * x - 12 * NXYZ4 * x + 4 * x * x)) / 6
    +((-x + y) * (1 - x + y)) /2
    +(z - x);
";
sectionTemplateTotalElements=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return index_xyz<TEMPLATE_reduced_size_list, SignedT>(NXYZ4 - 1, NXYZ4, NXYZ4 + 1/*one-past-last*/);
";
AddSettings[sgname, sgops, baseasymunitexpr, selectasymunitid, sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];


(*space group 225, namely FCC*)
Clear[sgname,sgops,baseasymunitexpr,selectasymunitid,sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];
sgname="3D225";
sgops={
AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,1,-1}],{0,0,0}}],
AffineTransform[{{{0,0,1},{1,0,0},{0,1,0}},{0,0,0}}],
AffineTransform[{{{0,1,0},{1,0,0},{0,0,-1}},{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,-1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{1,1,1}],{0,1/2,1/2}}],
AffineTransform[{DiagonalMatrix[{1,1,1}],{1/2,0,1/2}}]
}
baseasymunitexpr=(0<=x<=1/2)&&(0<=y<=1/4)&&(0<=z<=1/4)&&(y<=Min[x,1/2-x])&&(z<=y);
(*select suitable asymmetric unit*)
selectasymunitid=48;
sectionTemplateElement=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return y <= z && y <= (2 * NXYZ4 - z - 1) && x <= y;
";
sectionTemplateIndexXYZ=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return (x * (2 + 6 * NXYZ4 + 3 * NXYZ4 * NXYZ4 - 3 * x - 3 * NXYZ4 * x + x * x))/3
    +(x - y) * (x + y - 1 - 2 * NXYZ4)
    +(z - y);
";
sectionTemplateTotalElements=
"
constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
return index_xyz<TEMPLATE_reduced_size_list, SignedT>(NXYZ4 - 1, NXYZ4 - 1, NXYZ4 + 1/*one-past-last*/);
";
AddSettings[sgname, sgops, baseasymunitexpr, selectasymunitid, sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];


(*plane group 9, namely C6 without triple axis.*)
(*we follow the rule that the extra dimension *)
Clear[sgname,sgops,baseasymunitexpr,selectasymunitid,sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];
sgname="2D009";
sgops={
AffineTransform[{DiagonalMatrix[{1,1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,-1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{-1,1,1}],{0,0,0}}],
AffineTransform[{DiagonalMatrix[{1,1,1}],{1/2,1/2,0}}]
}
baseasymunitexpr=(0<=x<=1/4)&&(0<=y<=1/2)&&(0<=z<=1);
(*select suitable asymmetric unit*)
selectasymunitid=1;
sectionTemplateElement=
"
return true;
";
sectionTemplateIndexXYZ=
"
constexpr SignedT NX4 = (SignedT)(NX / 4);
constexpr SignedT NY2 = (SignedT)(NY / 2);
return (x * NY2 + y) * NZ + z;
";
sectionTemplateTotalElements=
"
constexpr SignedT NX4 = (SignedT)(NX / 4);
constexpr SignedT NY2 = (SignedT)(NY / 2);
return index_xyz<TEMPLATE_reduced_size_list, SignedT>(NX4 - 1, NY2 - 1, NZ - 1 + 1/*one-past-last*/);
";
AddSettings[sgname, sgops, baseasymunitexpr, selectasymunitid, sectionTemplateElement, sectionTemplateIndexXYZ, sectionTemplateTotalElements];


(* ::Text:: *)
(* *)


settingsForGroupInfoGenerator::undefined="The settings for required group `1` is not defined in file "<>scriptPath;
Table[
If[
KeyExistsQ[settings,i],
GroupInfoGenerator[templateFilename, outputFoldername, ##]&@@settings[i],
Message[settingsForGroupInfoGenerator::undefined,i];
Exit[10]
]
,{i,requiredGroupList}]
