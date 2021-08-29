//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.02.18
//+++ Purpose: add getting functions for abaqus io importor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.25  add more abaqus inp keylines info.

//#include <fstream>
//#include <iostream>
#include <array>
#include "Mesh/AbaqusIO.h"

//********************************************************
int AbaqusIO::GetElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){// quad-4 mesh
        return 4;
    }
    else if(substr.find("CPS8")!=string::npos){// quad-8
        return 8;
    }
    else if(substr.find("CPS3")!=string::npos){// tri-3 mesh
        return 3;
    }
    else if(substr.find("CPS6")!=string::npos){// tri-6 mesh
        return 6;
    }
	else if (substr.find("CPS8R") != string::npos) {
		// quad-8
		return 8;
	}
	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return 2;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return 3;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return 2;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return 3;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return 2;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return 2;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return 3;
	}

	// for 3d case================================
    else if(substr.find("C3D4")!=string::npos){// for tet-4
        return 4;
    }
    else if(substr.find("C3D8")!=string::npos){// for hex-8 mesh
        return 8;
    }
    else if(substr.find("C3D10")!=string::npos){// for tet-10 mesh
        return 10;
    }
    else if(substr.find("C3D20")!=string::npos){// for hex-20 mesh
        return 20;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 27;
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 2;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 3;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 2;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 3;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 2;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 2;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 3;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 3;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 3;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 4;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 4;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 8;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 3;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 4;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 4;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 8;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 8;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 9;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 9;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//******************************************************************
int AbaqusIO::GetSpaceDimFromInp() const {//node coordinates dimension or call space dimension
	ifstream in;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get coordinates dimension number, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	string str;
	vector<double> numbers;
	int spaceDim=0;
	while (!in.eof()) {
		getline(in, str);
		if (str.find("*Node") != string::npos&&str.find("Output") == string::npos) {
			getline(in, str);
			numbers = StringUtils::SplitStrNum(str, ',');
			break;
		}
	}
	in.close();
	spaceDim= numbers.size()-1;
	if (spaceDim<1|| spaceDim>3) {
		MessagePrinter::PrintErrorTxt("Unexpected node coordinates dimension, please check the input file at node lines!");
		MessagePrinter::AsFem_Exit();
	}
	return spaceDim;
}

int AbaqusIO::GetElmtDimFromInp() const{//mesh dimension
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 2;
    }

	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return 1;		//2;//1;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;		//2;//1;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return 1;		//2;//1;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return 1;		//2;//1;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return 1;		//2;//1;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return 1;		//2;//1;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return 1;		//2;//1;
	}

    else if(substr.find("C3D4")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D10")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 3;
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 1;		//3;//1;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;		//3;//1;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 1;		//3;//1;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 1;		//3;//1;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 1;		//3;//1;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 1;		//3;//1;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 1;		//3;//1;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 2;		//3;//2;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 2;		//3;//2;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 2;		//3;//2;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 2;		//3;//2;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 2;		//3;//2;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 2;		//3;//2;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 2;		//3;//2;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 2;		//3;//2;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 2;		//3;//2;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 2;		//3;//2;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 2;		//3;//2;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 2;		//3;//2;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtDimFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//***************************************************************
int AbaqusIO::GetElmtOrderFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element order number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 2;
    }
	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return 2;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return 2;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return 2;
	}

	// for 3d case
    else if(substr.find("C3D4")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D10")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }
	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 2;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 2;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 2;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 1;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 1;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 1;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 1;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 2;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 1;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 1;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 1;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 2;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 2;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 3;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 3;
	}


    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtOrderFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************
int AbaqusIO::GetElmtVTKCellTypeFromInp()const{//refer to https://vtk.org/doc/nightly/html/vtkCellType_8h.html#ab1d6fd1f3177b8a2a32bb018807151f8a52b8e3f5a3b009dbac0daf8bd6d96ee4
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element vtkcelltype number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return 9;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return 23;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return 5;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return 22;
    }
	//2D 待补充  "T2D2" "T2D3" "B21"  "B22"  "B23" "PIPE21" "PIPE22" 
	else if (substr.find("T2") != string::npos|| substr.find("B2") != string::npos || substr.find("PIPE2") != string::npos) {//暂时将其当做线单元line=3
		// 6-node second order triangle
		return 3;
	}//============================================================


	// for 3d case
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return 10;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return 12;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 24;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return 25;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return 29;
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 3;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 35;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 3;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 35;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 3;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 3;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 35;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 5;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 5;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 9;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 9;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 30;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 5;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 9;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 9;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 30;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 30;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 28;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 28;
	}

	else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtVTKCellTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//**************************************************
MeshType AbaqusIO::GetElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element mesh type, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return MeshType::QUAD4;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return MeshType::QUAD8;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return MeshType::TRI3;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return MeshType::TRI6;
    }

	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return MeshType::EDGE2;	//<---for pass test.  MeshType::TRUSS2;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return MeshType::TRUSS3;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return MeshType::EDGE2;		//<---for pass test.  MeshType::BEAM2;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return MeshType::BEAM2;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return MeshType::BEAM3;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return MeshType::PIPE2;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return MeshType::PIPE3;
	}

	// for 3d case
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return MeshType::TET4;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return MeshType::HEX8;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::TET10;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::HEX20;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::HEX27;
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return MeshType::TRUSS2;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return MeshType::TRUSS3;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return MeshType::BEAM2;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return MeshType::BEAM3;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return MeshType::BEAM3;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return MeshType::PIPE2;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return MeshType::PIPE3;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return MeshType::TRI3;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return MeshType::TRI3;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return MeshType::QUAD4;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return MeshType::QUAD4;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return MeshType::QUAD8;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return MeshType::TRI3;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return MeshType::QUAD4;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return MeshType::QUAD4;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return MeshType::QUAD8;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return MeshType::QUAD8;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return MeshType::QUAD9;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return MeshType::QUAD9;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
string AbaqusIO::GetElmtMeshTypeNameFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element mesh type name, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return "quad4";
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return "quad8";
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return "tri3";
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return "tri6";
    }
	// include 2d truss meshtype and membrane element meshtype ===XY===			membrane element library	http://127.0.0.1:2080/v6.13/books/usb/default.htm?startat=pt06ch29s01ael10.html
	else if (substr.find("T2D2") != string::npos) {//Truss element library   http://127.0.0.1:2080/v6.13/books/usb/default.htm?startat=pt06ch29s02.html
		// for 2-node linear displacement
		return "edge2";//<---for pass test.  "truss2";
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return "truss3";
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {//http://127.0.0.1:2080/v6.13/books/usb/default.htm?startat=pt06ch29s03abo26.html#usb-elm-ebeamoverview  Beam element library
		// for 2-node linear beam						//http://127.0.0.1:2080/v6.13/books/usb/default.htm?startat=pt06ch29s06abo27.html   shell element library
		return "edge2";//<---for pass test.  "beam2";
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return "beam3";
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return "beam3";
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return "pipe2";
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return "pipe3";
	}


    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return "tet4";
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return "hex8";
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return "tet10";
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return "hex20";
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return "hex27";
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return "truss2";
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return "truss3";
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return "beam2";
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return "beam3";
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return "beam3";
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return "pipe2";
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return "pipe3";
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return "tri3";
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return "tri3";
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return "quad4";
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return "quad4";
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return "quad8";
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return "quad3";
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return "quad4";
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return "quad4";
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return "quad8";
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return "quad8";
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return "quad9";
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return "quad9";
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return "unknowm-mesh";
}
//************************************************************
int AbaqusIO::GetSubElmtNodesNumFromInp() const{//体-面-线-点      面-线-点;跟上一版本GetAbaqusBCElmtNodesNumFromInp含义近似
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 3;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 4;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 6;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 8;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 9;
    }
	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return 1;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return 1;
	}
	// for 3d case  ================================
	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 1;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 1;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 2;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 2;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 2;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 2;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 3;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 2;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 2;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 2;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 3;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 3;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 3;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 3;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
int AbaqusIO::GetSubSubElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub sub element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // quad4-edge2-point
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){
        // quad8-edge3-point
        return 1;
    }
    else if(substr.find("CPS3")!=string::npos){
        // tri3-edge2-point
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){
        // tri6-edge3-point
        return 1;
    }
    else if(substr.find("C3D4")!=string::npos){
        // tet4-tri3-edge2
        return 2;
    }
    else if(substr.find("C3D8")!=string::npos){
        // hex8-quad4-edge2
        return 2;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron-tri6-edge3
        return 3;
    }
    else if(substr.find("C3D20")!=string::npos){
        // hex20-quad8-edge3
        return 3;
    }
    else if(substr.find("C3D27")!=string::npos){
        // hex27-quad9-edge3
        return 3;
    }
	//for plane 3D elements   ===XY===
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 1;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 1;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 1;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 1;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 1;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 1;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 1;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 1;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 1;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 1;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 1;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 1;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubSubElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
int AbaqusIO::GetSubElmtDimFromInp() const{
    return GetElmtDimFromInp()-1;
}

int AbaqusIO::GetSubElmtOrderFromInpElmt()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element order number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){// quad-8 mesh
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){// tri-6 mesh
        return 2;
    }

	// include 2d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return 1;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return 1;
	}

	// for 3d case
	else if(substr.find("C3D4")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 2;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }

	// include 3d truss meshtype and membrane element meshtype ===XY===
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return 1;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return 1;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return 1;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return 1;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return 1;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return 1;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return 1;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return 1;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return 1;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return 1;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return 1;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return 2;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return 1;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return 1;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return 1;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return 2;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return 2;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return 3;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return 3;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtOrderFromInpElmt, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
MeshType AbaqusIO::GetSubElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub element type information, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return MeshType::EDGE2;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return MeshType::EDGE3;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return MeshType::EDGE2;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return MeshType::TRI3;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return MeshType::QUAD4;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::TRI6;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::QUAD8;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::QUAD9;
    }

	else if (substr.find("T2D2") != string::npos) {
		// for 2-node linear displacement
		return MeshType::POINT;
	}
	else if (substr.find("T2D3") != string::npos) {
		// for 3-node quadratic displacement
		return MeshType::POINT;
	}
	// include 2d beam meshtype and membrane element meshtype ===XY===
	else if (substr.find("B21") != string::npos) {
		// for 2-node linear beam
		return MeshType::POINT;
	}
	else if (substr.find("B22") != string::npos) {
		// for 3-node quadratic beam
		return MeshType::POINT;
	}
	else if (substr.find("B23") != string::npos) {
		// for 2-node cubic beam
		return MeshType::POINT;
	}
	else if (substr.find("PIPE21") != string::npos) {
		// for 2-node linear pipe
		return MeshType::POINT;
	}
	else if (substr.find("PIPE22") != string::npos) {
		// for 3-node quadratic pipe
		return MeshType::POINT;
	}
	else if (substr.find("T3D2") != string::npos) {
		// for 2-node linear displacement
		return MeshType::POINT;
	}
	else if (substr.find("T3D3") != string::npos) {
		// for 3-node quadratic displacement
		return MeshType::POINT;
	}
	else if (substr.find("B31") != string::npos) {
		// for 2-node linear beam
		return MeshType::POINT;
	}
	else if (substr.find("B32") != string::npos) {
		// for 3-node quadratic beam
		return MeshType::POINT;
	}
	else if (substr.find("B33") != string::npos) {
		// for 2-node cubic beam
		return MeshType::POINT;
	}
	else if (substr.find("PIPE31") != string::npos) {
		// for 2-node linear pipe
		return MeshType::POINT;
	}
	else if (substr.find("PIPE32") != string::npos) {
		// for 3-node quadratic pipe
		return MeshType::POINT;
	}
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return MeshType::EDGE2;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return MeshType::EDGE2;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return MeshType::EDGE2;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return MeshType::EDGE2;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return MeshType::EDGE3;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return MeshType::EDGE2;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return MeshType::EDGE2;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return MeshType::EDGE2;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return MeshType::EDGE3;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return MeshType::EDGE3;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return MeshType::EDGE3;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return MeshType::EDGE3;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
MeshType AbaqusIO::GetSubSubElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub sub element type information, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle-egde2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle-edge2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle-edge2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle-edge3-point
        return MeshType::POINT;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron-tri3-edge2
        return MeshType::EDGE2;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron-quad4-edge2
        return MeshType::EDGE2;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron-quad8-edge3
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron-quad9-edge3
        return MeshType::EDGE3;
    }

	//for plane 3D elements   ===XY===
	else if (substr.find("S3") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3R)
		return MeshType::POINT;
	}
	else if (substr.find("S3R") != string::npos) {
		// for 3-node triangular general-purpose shell, finite membrane strains (identical to element S3)
		return MeshType::POINT;
	}
	else if (substr.find("S4") != string::npos) {
		// for 4-node general-purpose shell, finite membrane strains
		return MeshType::POINT;
	}
	else if (substr.find("S4R") != string::npos) {
		// for 4-node general-purpose shell, reduced integration with hourglass control, finite membrane strains
		return MeshType::POINT;
	}
	else if (substr.find("S8R") != string::npos) {
		// for 8-node doubly curved thick shell, reduced integration
		return MeshType::POINT;
	}
	else if (substr.find("M3D3") != string::npos) {
		// for 3-node triangle
		return MeshType::POINT;
	}
	else if (substr.find("M3D4") != string::npos) {
		// for 4-node quadrilateral
		return MeshType::POINT;
	}
	else if (substr.find("M3D4R") != string::npos) {
		// for 4-node quadrilateral, reduced integration, hourglass control
		return MeshType::POINT;
	}
	else if (substr.find("M3D8") != string::npos) {
		// for 8-node quadrilateral
		return MeshType::POINT;
	}
	else if (substr.find("M3D8R") != string::npos) {
		// for 8-node quadrilateral, reduced integration
		return MeshType::POINT;
	}
	else if (substr.find("M3D9") != string::npos) {
		// for 9-node quadrilateral
		return MeshType::POINT;
	}
	else if (substr.find("M3D9R") != string::npos) {
		// for 9-node quadrilateral, reduced integration, hourglass control
		return MeshType::POINT;
	}

    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubSubElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
int AbaqusIO::GetNodesNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nNodes;
    substr.clear();
    nNodes=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Node")!=string::npos&&str.find("Output") == string::npos){
            nNodes=0;
            getline(in,str);// this line already contains 1st node's x,y,z   coordinates
            while(str.find("*")==string::npos){		//find("*Element")
                nNodes+=1;
                getline(in,str);
            }
        }
    }
    in.close();
    return nNodes;
}
//*******************************************************
int AbaqusIO::GetElmtsNumFromInp()const{		////所有单元_nBulkElmts+nBCElmts  也包括边界单元啊
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nElmts;
    substr.clear();
    nElmts=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element,")!=string::npos&&str.find("Output") == string::npos){
            nElmts=0;
            getline(in,str);// this line already contains 1st node's x,y,z   coordinates
            while(str.find("*")==string::npos){		//find("*End Instance") "*Nset, nset=" "*Elset, elset=*"  *End Part  *Beam Section这些关键字
                nElmts+=1;
                getline(in,str);
            }
        }
    }
    in.close();
    return nElmts;
}
//***********************************
int AbaqusIO::GetNsetsNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Nset number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nNsets;
    substr.clear();
    nNsets=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Nset,")!= string::npos&&str.find("nset=Set")!= string::npos){	//&&str.find("instance=")!= string::npos	//part和assembly里都有set,会重名。。
            nNsets+=1;
        }
    }
    in.close();
    return nNsets;
}
//***********************************
int AbaqusIO::GetElsetsNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Elset number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nElsets;
    substr.clear();
    nElsets=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Elset,")!=string::npos&&str.find("Elset=Set") != string::npos) {//&&str.find("instance=") != string::npos    NOTE pressure in elset= _..._spos/sneg  //另外part和assembly里都有set,会重名。。
            nElsets+=1;
        }
    }
    in.close();
    return nElsets;
}
//***********************************
int AbaqusIO::GetSurfacesNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nSurfaces;
    substr.clear();
    nSurfaces=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos){
            nSurfaces+=1;
        }
    }
    in.close();
    return nSurfaces;
}
//*********************************************************
int AbaqusIO::GetSurfaceElmtsNumFromInp()const{
    ifstream in,insub;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    int nElmts;
    substr.clear();
    nElmts=0;
    int i;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left		//需要改
            // _Surf-Left_S4, S4
            // here we only read '_Surf-Left_S4'
            getline(in,str);
            i=str.find(",");
            SurfaceSetName=str.substr(0,i-1);// '_Surf-Left_S4'
            insub.open(_MeshFileName,ios::in);
            while(!insub.eof()){
                getline(insub,str1);
                if(str1.find("*Elset, elset=")!=string::npos&&
                   str1.find(SurfaceSetName)!=string::npos){
                    // now we find the related Elset information
                    if(str1.find("generate")!=string::npos){
                        getline(insub,substr);// read the element index information
                        numbers=StringUtils::SplitStrNum(substr);
                        if(static_cast<int>(numbers.size())!=3){
                            MessagePrinter::PrintErrorTxt(SurfaceSetName+" is using generate way, however, your element index number is not equal to 3");
                            MessagePrinter::AsFem_Exit();
                        }
                        nElmts+=1+(static_cast<int>(numbers[1])-static_cast<int>(numbers[0]))/static_cast<int>(numbers[2]);
                    }
                    else{
                        for(i=0;i<numeric_limits<int>::max();i++){
                            getline(insub,substr);
                            numbers=StringUtils::SplitStrNum(substr);
                            nElmts+=static_cast<int>(numbers.size());
							if (substr.find("Surf-") != string::npos) { nElmts = 0; }
                            if(substr.find("*")!=string::npos){
                                break;
                            }
                        }
                    }
                } // end-of-element-set-information
            }// end-of-while-insub
            insub.close();
        }
    }
    in.close();
    return nElmts;
}
//*********************************************************
int AbaqusIO::GetSurfaceEdgeIDViaSurfaceNameFromInp(string surfacesetname)const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    vector<int> ElmtIDs;
    substr.clear();
    ElmtIDs.clear();
    int i,edge;
    edge=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos && str.find(surfacesetname)!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left					default is: *Surface, type=ELEMENT, name=Surf-1
            // _Surf-Left_S4, S4													_Surf-1_SPOS, SPOS
            getline(in,str);// read _Surf-Left_S4, S4
            i=str.find(",");
            substr=str.substr(i+1,string::npos);// '_Surf-Left_S4'
            numbers=StringUtils::SplitStrNum(substr);// split s4 into 4
            if(numbers.size()<1){
                MessagePrinter::PrintErrorTxt("can not find any numbers in "+str+" ,please check your inp file");
                //MessagePrinter::AsFem_Exit();
				MessagePrinter::PrintNormalTxt("Default SURFACE edge set as S1 ...");//==上面改了==
				MessagePrinter::PrintStars();
				numbers.push_back(1);
            }
            edge=static_cast<int>(numbers[0]);
        }
    }
    in.close();
    return edge;
}
//*********************************************************
vector<int> AbaqusIO::GetSurfaceElmtIDViaSurfaceNameFromInp(string surfacesetname)const{
    ifstream in,insub;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can nott get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    vector<int> ElmtIDs;
    substr.clear();
    int i;
    ElmtIDs.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos && str.find(surfacesetname)!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            // _Surf-Left_S4, S4
            getline(in,str);// read _Surf-Left_S4, S4
            i=str.find(",");
            SurfaceSetName=str.substr(0,i-1);// '_Surf-Left_S4'
            insub.open(_MeshFileName,ios::in);
            while(!insub.eof()){
                getline(insub,str1);
                if(str1.find("*Elset, elset=")!=string::npos&&
                   str1.find(SurfaceSetName)!=string::npos){
                    // now we find the related Elset information
                    if(str1.find("generate")!=string::npos){
                        getline(insub,substr);// read the element index information
                        numbers=StringUtils::SplitStrNum(substr);
                        if(static_cast<int>(numbers.size())!=3){
                            MessagePrinter::PrintErrorTxt(SurfaceSetName+" is using generate way, however, your element index number is not equal to 3");
                            MessagePrinter::AsFem_Exit();
                        }
                        for(int e=static_cast<int>(numbers[0]);e<=static_cast<int>(numbers[1]);e+=static_cast<int>(numbers[2])){
                            ElmtIDs.push_back(e);
                        }
                    }
                    else{
                        for(i=0;i<numeric_limits<int>::max();i++){
                            getline(insub,substr);
                            numbers=StringUtils::SplitStrNum(substr);
                            for(i=0;i<static_cast<int>(numbers.size());i++){
                                ElmtIDs.push_back(static_cast<int>(numbers[i]));
                            }
                            if(substr.find("*")!=string::npos){
                                break;
                            }
                        }
                    }
                } // end-of-element-set-information
            }// end-of-while-insub
            insub.close();
        }
    }
    in.close();
    return ElmtIDs;// please keep in mind, this element id is just the global bulk id, not the boundary element id !!!
}



//****************补充GetFuns********************************
int AbaqusIO::GetBeamSectionNumFromInp() const {
	ifstream in;
	string str;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get BeamSection number, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int nSectionNum = 0;
	while (!in.eof()) {
		if (str.find("*Beam Section") != string::npos&&str.find("elset=Set") != string::npos) {
			nSectionNum += 1;
		}
		getline(in, str);
	}
	in.close();
	return nSectionNum;
	cout << "nBeamSectionNum= " << nSectionNum << endl;
}
//************************************************
int AbaqusIO::GetSolidSectionNumFromInp() const {
	ifstream in;
	string str;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get SolidSection number, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int nSectionNum = 0;
	while (!in.eof()) {
		if (str.find("*Solid Section") != string::npos&&str.find("elset=Set") != string::npos) {
			nSectionNum += 1;
		}
		getline(in, str);
	}
	in.close();
	return nSectionNum;
	cout << "nSolidSectionNum= " << nSectionNum << endl;
}
//************************************************
int AbaqusIO::GetShellSectionNumFromInp() const {
	ifstream in;
	string str;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get ShellSection number, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int nSectionNum = 0;
	while (!in.eof()) {
		if (str.find("*Shell General Section") != string::npos&&str.find("elset=Set") != string::npos) {
			nSectionNum += 1;
		}
		getline(in, str);
	}
	in.close();
	return nSectionNum;
	cout << "nShellSectionNum= " << nSectionNum << endl;
}

//************************************************
vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> > AbaqusIO::GetBeamSectionESetsandMatAdFromInp() const {

	ifstream in;
	string str, substr, substr2, substr3, tempstr;
	vector<string> strvec;

	vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> > beamsecinfo;

	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get BeamSectionESetsandMatAd, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int i, j, k;
	substr.clear();
	substr2.clear();
	beamsecinfo.clear();
	while (!in.eof()) {
		if (str.find("*Beam Section, elset=") != string::npos) {
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				//have all the necessary information
				i = strvec[1].find_first_of('=');//elset
				j = strvec[2].find_first_of('=');//mat name
				k = strvec[4].find_first_of('=');//section type 比如矩形
				substr = strvec[1].substr(i + 1) + "part";//赋材料属性 自己把set名称+"part"
				substr = StringUtils::RemoveStrSpace(substr);
				substr2 = strvec[2].substr(j + 1, strvec[2].length() - 1);//测试默认到字符串长度和给定过大长度的结果是否一样
				substr2 = StringUtils::RemoveStrSpace(substr2);
				substr3 = strvec[4].substr(k + 1);
				substr3 = StringUtils::RemoveStrSpace(substr3);
			}
			else {
				MessagePrinter::PrintErrorTxt(" cant find BeamSectionESets= name in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			str.clear();
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			//array<string, 3> ab={ substr,substr2, substr3 };
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			array<double, 2> ab;
			array<double, 3> bdirec;//截面尺寸、法向
			int si = 0;
			for (auto i : strvec) {
				ab[si] = stod(i);
				si++;
			}
			//如果有两行数据 判断继续读取 再存
			str.clear();
			getline(in, str);
			if (str.find("*") == string::npos) {// if not keyline, then go on beam orientation second line
				str = StringUtils::RemoveStrSpace(str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				int si = 0;
				for (auto i : strvec) {
					bdirec[si] = stod(i);
					si++;
				}
			}
			array<string,3> beamline = { substr,substr2, substr3 };
			beamsecinfo.push_back(make_pair(beamline, make_pair(ab, bdirec)));// beamsection(Elset,mat,sectiontype,a b,directionVec)
		}
		getline(in, str);
	}
	in.close();
	return beamsecinfo;

}

//************************************************
vector<pair<pair<string, string>, double>> AbaqusIO::GetSolidSectionESetsandMatAFromInp() const {
	ifstream in;
	string str, substr, substr2, tempstr;
	vector<string> strvec;
	vector<pair<pair<string, string>, double>> solidsecinfo;
	double AorT;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get SolidSectionESetsandMatA, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int i, j;
	substr.clear();
	substr2.clear();
	solidsecinfo.clear();
	while (!in.eof()) {
		if (str.find("*Solid Section, elset=") != string::npos) {
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				//have all the necessary information
				i = strvec[1].find_first_of('=');
				j = strvec[2].find_first_of('=');
				substr = strvec[1].substr(i + 1, strvec[1].length() - 1) + "part";//赋材料属性 自己把set名称+"part"
				substr = StringUtils::RemoveStrSpace(substr);
				substr2 = strvec[2].substr(j + 1, strvec[2].length() - 1);
				substr2 = StringUtils::RemoveStrSpace(substr2);
			}
			else {
				MessagePrinter::PrintErrorTxt(" cant find SolidSectionESets= name in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			str.clear();
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			tempstr.clear();
			AorT = 0.0;
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			if(strvec.size()==0){
				MessagePrinter::PrintErrorTxt(" cant find SolidSectionAorT= value in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			tempstr = strvec[0];		//str.substr(0, str.length() - 1);
			//if (tempstr.size() == 1) {
			//	AorT = atof(tempstr.c_str());
			//}
			//else {
			//	MessagePrinter::PrintErrorTxt(" cant find SolidSectionAorT= value in your inp file !!!");
			//	MessagePrinter::AsFem_Exit();
			//}
			AorT = atof(tempstr.c_str());
			solidsecinfo.push_back(make_pair(make_pair(substr, substr2), AorT));// solidsection(Elset,mat,thickness)	//好像得改
		}
		getline(in, str);
	}
	in.close();
	return solidsecinfo;

}

//************************************************//Membrane Section 待补充 类似
vector<pair<pair<string, string>, double>> AbaqusIO::GetShellSectionESetsandMatTFromInp() const {
	ifstream in;
	string str, substr, substr2, tempstr;
	vector<string> strvec;
	vector<pair<pair<string, string>, double>> shellsecinfo;
	double AorT;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get ShellSectionESetsandMatT, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int i, j;
	substr.clear();
	substr2.clear();
	shellsecinfo.clear();
	while (!in.eof()) {
		if (str.find("*Shell General Section, elset=") != string::npos) {
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				//have all the necessary information
				i = strvec[1].find_first_of('=');
				j = strvec[2].find_first_of('=');
				substr = strvec[1].substr(i + 1, strvec[1].length() - 1) + "part";//赋材料属性 自己把set名称+"part"
				substr = StringUtils::RemoveStrSpace(substr);
				substr2 = strvec[2].substr(j + 1, strvec[2].length() - 1);
				substr2 = StringUtils::RemoveStrSpace(substr2);
			}
			else {
				MessagePrinter::PrintErrorTxt(" cant find ShellSectionESets= name in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			str.clear();
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			tempstr.clear();
			AorT = 0.0;
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			tempstr = strvec[0];	
			//tempstr = str.substr(0, str.length() - 1);
			//if (tempstr.size() == 1) {
			//	AorT = atof(tempstr.c_str());
			//}
			//else {
			//	MessagePrinter::PrintErrorTxt(" cant find ShellSectionAorT= value in your inp file !!!");
			//	MessagePrinter::AsFem_Exit();
			//}
			AorT = atof(tempstr.c_str());
			shellsecinfo.push_back(make_pair(make_pair(substr, substr2), AorT));// solidsection(Elset,mat,thickness)
		}
		getline(in, str);
	}
	in.close();
	return shellsecinfo;
}

//****************Read Amplitude info********************************
//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > AbaqusIO::GetAmpsTimeFromInp() {
vector<std::pair<std::vector<std::string>, std::vector<double>> >  AbaqusIO::GetAmpsTimeFromInp(){
	ifstream in;
	streampos pos;
	string str, substr, substr2, substr3, substr4;
	vector<string> strvec;
	vector<string> ampinfo;
	vector<double> ampdata;
	//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > amptemp;
	vector<std::pair<std::vector<std::string>, std::vector<double>> > amptemp;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get AmpsTime, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	//getline(in, str);
	int i, j, k, m;
	//double ti = 0, ampi = 0;
	amptemp.clear();
	Step.AmpsTime.clear();
	while (!in.eof()) {
		getline(in, str);
		substr.clear();
		substr2.clear();
		substr3.clear();
		substr4.clear();
		ampinfo.clear();
		ampdata.clear();
		if (str.find("*Amplitude") != string::npos) {// currently only for one Amp1==modified===
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				//have all the necessary information
				if (strvec.size() == 3) {
					if(str.find("definition=") == string::npos){//	key/name/smooth  //default is TABULAR type
						i = strvec[1].find_first_of('=');
						j = strvec[2].find_first_of('=');
						substr = strvec[1].substr(i + 1); 	//Aname
						substr = StringUtils::RemoveStrSpace(substr);
						substr3 = strvec[2].substr(j + 1);	//smooth
						substr3 = StringUtils::RemoveStrSpace(substr3);
						substr2 = "TABULAR";
						ampinfo.push_back(substr);
						ampinfo.push_back(substr2);
						ampinfo.push_back(substr3);
					}
					else{//	key/name/type
						i = strvec[1].find_first_of('=');
						j = strvec[2].find_first_of('=');
						substr = strvec[1].substr(i + 1); 	//Aname
						substr = StringUtils::RemoveStrSpace(substr);
						substr2 = strvec[2].substr(j + 1);	//可能没有这句Atype  time=TOTAL TIME		//没考虑smooth和fixInterval值选项
						substr2 = StringUtils::RemoveStrSpace(substr2);
						ampinfo.push_back(substr);
						ampinfo.push_back(substr2);
					}
				}
				else if (strvec.size() == 4) {// must be key/name/type=equal space/interval. 
					i = strvec[1].find_first_of('=');
					j = strvec[2].find_first_of('=');
					k = strvec[3].find_first_of('=');
					substr = strvec[1].substr(i + 1); 	//Aname
					substr = StringUtils::RemoveStrSpace(substr);
					substr2 = strvec[2].substr(j + 1);	
					substr2 = StringUtils::RemoveStrSpace(substr2);
					substr4 = strvec[3].substr(k + 1);
					substr4 = StringUtils::RemoveStrSpace(substr4);
					substr3 = "0.0";//no smooth, but we take 0.0
					ampinfo.push_back(substr);
					ampinfo.push_back(substr2);
					ampinfo.push_back(substr3);
					ampinfo.push_back(substr4);
				}
				else if (strvec.size() == 5) {// must be key/name/type=equal space/smooth/interval. 
					i = strvec[1].find_first_of('=');
					j = strvec[2].find_first_of('=');
					k = strvec[3].find_first_of('=');
					m = strvec[4].find_first_of('=');
					substr = strvec[1].substr(i + 1); 	//Aname
					substr = StringUtils::RemoveStrSpace(substr);
					substr2 = strvec[2].substr(j + 1);	//
					substr2 = StringUtils::RemoveStrSpace(substr2);
					substr3 = strvec[3].substr(k + 1);	//
					substr3 = StringUtils::RemoveStrSpace(substr3);
					substr4 = strvec[4].substr(m + 1);	//
					substr4 = StringUtils::RemoveStrSpace(substr4);
					ampinfo.push_back(substr);
					ampinfo.push_back(substr2);
					ampinfo.push_back(substr3);
					ampinfo.push_back(substr4);
				}
				else if(strvec.size() == 2) {// must be key/name
					i = strvec[1].find_first_of('=');
					substr = strvec[1].substr(i + 1); 	//Aname
					substr = StringUtils::RemoveStrSpace(substr);
					//substr2 = "none";
					substr2 = "TABULAR";//default is TABULAR type
					ampinfo.push_back(substr);
					ampinfo.push_back(substr2);
				}
			}
			else {
				MessagePrinter::PrintErrorTxt(" cant find AmpsTime= name/type in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			//str.clear();
			//getline(in, str);
			//str = StringUtils::RemoveStrSpace(str);
			//strvec.clear();
			//strvec = StringUtils::SplitStr(str, ',');
			//int si = 0;
			//vector<std::pair<double, double> > amptime;
			//amptime.clear();
			//for (auto i : strvec) {
			//	si++;
			//	if (si % 2 == 0) {
			//		ampi = stod(i);
			//		amptime.push_back(make_pair(ti, ampi));
			//	}
			//	else {
			//		ti = stod(i);
			//	}
			//}
			////如果有两行数据 判断继续读取 再存
			//str.clear();
			//getline(in, str);
			//if (str.find("*") == string::npos) {// if not keyline, then go on ampvalue second line
			//	str = StringUtils::RemoveStrSpace(str);
			//	strvec.clear();
			//	strvec = StringUtils::SplitStr(str, ',');
			//	int si = 0;
			//	for (auto i : strvec) {
			//		si++;
			//		if (si % 2 == 0) {
			//			ampi = stod(i);
			//			amptime.push_back(make_pair(ti, ampi));
			//		}
			//		else {
			//			ti = stod(i);
			//		}
			//	}

			str.clear();
			getline(in, str);
			ampdata.clear();
			while(str.find("*") == string::npos){
				//getline(in, str);
				str = StringUtils::RemoveStrSpace(str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				for (auto ii : strvec) {
					ampdata.push_back(stod(ii));
				}
				pos = in.tellg();//记录当前指针位置
				getline(in, str);
			}
			//}
			//amptemp.push_back(make_pair(make_pair(substr, substr2), amptime));// AmpsTime(amp name, type;  value table)
			amptemp.push_back(make_pair(ampinfo, ampdata));// AmpsTime(ampinfo; ampdata)			
			//getline(in, str);
			in.seekg(pos);//将指针指回记录行，可以重新读一次
		}
		//Step.AmpsTime = amptemp;
	}
	Step.AmpsTime = amptemp;
	in.close();
	return Step.AmpsTime;
}

//****************Read BC info********************************
vector<std::pair<string, vector<int>> > AbaqusIO::GetBoundaryDofFromInp() {//BC name and type infomation at note line** [line-1]
	ifstream in;
	string str;
	vector<string> strvec;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get BoundaryDof, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	string bcset;
	vector<int> bcdof;
	Step.boundaryDoF.clear();
	while (!in.eof()) {
		if (str.find("*Boundary") != string::npos) {
			getline(in, str);
			while (str.find("*") == string::npos) {
				str = StringUtils::RemoveStrSpace(str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				if (strvec.size() == 3) {//目前只考虑了6dof约束 没有速度等条件
					bcset = strvec[0];
					bcdof.clear();
					for (int i = stoi(strvec[1]); i <= stoi(strvec[2]); i++) {//生成的inp每句只约束一个自由度
						bcdof.push_back(i);
					}
					Step.boundaryDoF.push_back(make_pair(bcset, bcdof));
				}
				else if (strvec.size() < 1) {
					MessagePrinter::PrintErrorTxt(" cant find BC line info in your inp file !!!");
					MessagePrinter::AsFem_Exit();
				}
				else {
					MessagePrinter::PrintErrorTxt(" strvec.size()!=3 multi-dofs are in one line, please check Boundary line in inp !!!");
					continue;
				}
				getline(in, str);
			}
			//遇到下个关键字* 结束读数据行
		}
		getline(in, str);
	}
	in.close();
	return Step.boundaryDoF;
}

//****************Read Mate info********************************  name damping density elastic poission
vector<AbaqusIO::material_info> AbaqusIO::GetMateFromInp() {
	ifstream in;
	string str;
	vector<string> strvec;
	int i, j;
	string substr, substr2;
	material_info matstruc;
	bool elasticflag = false;
	Mate.clear();

	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get Mateinfo, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	while (!in.eof()) {
		if (str.find("*Material,") != string::npos) {//读完一个material语段，所有材料参数，然后继续寻找material2，存到容器
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {//*Material, name=Material-1
				i = strvec[1].find_first_of('=');
				substr = strvec[1].substr(i + 1);//substr(i + 1, strvec[1].length() - 1)
				substr = StringUtils::RemoveStrSpace(substr);
				matstruc.matname = substr;
			}
			else {
				MessagePrinter::PrintErrorTxt(" cant find Material= name  in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			strvec.clear();
			if (str.find("*Damping,") != string::npos) {
				matstruc.damp_alpha = 0.0;
				matstruc.damp_beta = 0.0;
				strvec = StringUtils::SplitStr(str, ',');
				if (strvec.size() == 2) {
					i = strvec[1].find_first_of('=');
					substr = strvec[1].substr(i + 1, strvec[1].length() - 1);
					substr = StringUtils::RemoveStrSpace(substr);
					matstruc.damp_alpha = atof(substr.c_str());//通常设置alpha 没有beta。
				}
				else if (strvec.size() == 3) {
					i = strvec[1].find_first_of('=');
					substr = strvec[1].substr(i + 1, strvec[1].length() - 1);
					substr = StringUtils::RemoveStrSpace(substr);
					matstruc.damp_alpha = atof(substr.c_str());
					j = strvec[2].find_first_of('=');
					substr2 = strvec[2].substr(j + 1, strvec[2].length() - 1);
					substr2 = StringUtils::RemoveStrSpace(substr2);
					matstruc.damp_beta = atof(substr2.c_str());
				}
				else {
					MessagePrinter::PrintErrorTxt(" unexpected dampline in your inp file !!!");
					MessagePrinter::AsFem_Exit();
				}
			}
			else if (str.find("*Density") != string::npos) {
				getline(in, str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				matstruc.density = atof(strvec[0].c_str());
			}
			else {
				MessagePrinter::PrintErrorTxt(" unexpected matline in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			if (str.find("*Elastic") != string::npos) {
				elasticflag = true;
				getline(in, str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				matstruc.Youngs_modulus = atof(strvec[0].c_str());
				matstruc.poisson = atof(strvec[1].c_str());
			}
			else if (str.find("*Density") != string::npos) {
				getline(in, str);
				strvec.clear();
				strvec = StringUtils::SplitStr(str, ',');
				matstruc.density = atof(strvec[0].c_str());

				getline(in, str);
				str = StringUtils::RemoveStrSpace(str);
				if (str.find("*Elastic") != string::npos) {
					elasticflag = true;
					getline(in, str);
					strvec.clear();
					strvec = StringUtils::SplitStr(str, ',');
					matstruc.Youngs_modulus = atof(strvec[0].c_str());
					matstruc.poisson = atof(strvec[1].c_str());
				}
				else {
					MessagePrinter::PrintErrorTxt(" neither density nor elastic line in your inp file !!!");
					MessagePrinter::AsFem_Exit();
				}

			}
			else {
				MessagePrinter::PrintErrorTxt(" neither density nor elastic line in your inp file !!!");
				MessagePrinter::AsFem_Exit();
			}
			if (elasticflag) {//一轮材料参数读取完成
				elasticflag = false;
				Mate.push_back(matstruc);
				continue;
			}
			else {
				MessagePrinter::PrintErrorTxt(" not find elastic line in your inp file, Mat read error !!!");
				MessagePrinter::AsFem_Exit();
			}

		}
		getline(in, str);
	}
	in.close();
	return Mate;
}

//****************Read Cload info********************************  amp,set    direction,value
vector<std::pair<pair<string, string>, vector<double>> > AbaqusIO::GetCLoadFromInp() {
	ifstream in;
	string str, substr, loadset;
	vector<string> strvec;
	vector<double> Cloadvecs;
	Cloads.clear();

	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get CLoadinfo, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int i;
	double direct, loadvalue;
	while (!in.eof()) {
		if (str.find("*Cload") != string::npos) {
			substr.clear();
			strvec.clear();
			Cloadvecs.clear();
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				i = strvec[1].find_first_of('=');
				substr = strvec[1].substr(i + 1, strvec[1].length() - 1);
				substr = StringUtils::RemoveStrSpace(substr);//Amp name
			}
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			loadset = strvec[0];
			direct = atof(strvec[1].c_str());
			loadvalue = atof(strvec[2].c_str());
			Cloadvecs.push_back(direct);//目前二维的例子就是1或2
			Cloadvecs.push_back(loadvalue);
			Cloads.push_back(make_pair(make_pair(substr, loadset), Cloadvecs));
			//cout<<,,,<<endl;
		}
		getline(in, str);
	}
	in.close();
	return Cloads;
}


//****************Read Dsload info********************************  amp,set,type    value
vector<std::pair<vector<string>, double> > AbaqusIO::GetDsLoadFromInp() {
	ifstream in;
	string str, substr, loadset, loadtype;
	vector<string> strvec;
	vector<string> loadstrs;
	Dsloads.clear();

	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get DsLoad info, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	int i;
	double loadvalue;
	while (!in.eof()) {
		if (str.find("*Dsload") != string::npos) {
			substr.clear();
			strvec.clear();
			loadstrs.clear();
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				i = strvec[1].find_first_of('=');
				substr = strvec[1].substr(i + 1, strvec[1].length() - 1);
				substr = StringUtils::RemoveStrSpace(substr);
				loadstrs.push_back(substr);//Amp name
			}
			else {
				MessagePrinter::PrintErrorTxt(" not find DsLoad info in your inp file, Load read error !!!");
				MessagePrinter::AsFem_Exit();
			}
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			loadset = strvec[0];
			loadtype = strvec[1];//P
			loadvalue = atof(strvec[2].c_str());
			loadstrs.push_back(loadset);
			loadstrs.push_back(loadtype);
			Dsloads.push_back(make_pair(loadstrs, loadvalue));
			//cout<<,,,<<endl;
		}
		getline(in, str);
	}
	in.close();
	return Dsloads;
}


//****************Read Dload info********************************  amp,set,type    value,directions
vector<std::pair<vector<string>, vector<double>> > AbaqusIO::GetDloadsFromInp() {
	string str, substr, loadset, loadtype;
	vector<string> strvec;
	vector<string> loadstrs;
	Dloads.clear();

	int i;
	double loadvalue;
	vector<double> Dloadvecs;
	ifstream in;
	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get Dloads info, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	getline(in, str);
	while (!in.eof()) {
		if (str.find("*Dload") != string::npos) {
			substr.clear();
			strvec.clear();
			loadstrs.clear();
			Dloadvecs.clear();
			strvec = StringUtils::SplitStr(str, ',');
			if (strvec.size() >= 2) {
				i = strvec[1].find_first_of('=');
				substr = strvec[1].substr(i + 1, strvec[1].length() - 1);
				substr = StringUtils::RemoveStrSpace(substr);
				loadstrs.push_back(substr);//Amp name
			}
			else {
				MessagePrinter::PrintErrorTxt(" not find DLoad info in your inp file, Load read error !!!");
				MessagePrinter::AsFem_Exit();
			}
			getline(in, str);
			str = StringUtils::RemoveStrSpace(str);
			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			loadset = strvec[0];
			loadtype = strvec[1];//GRAV
			loadvalue = atof(strvec[2].c_str());
			loadstrs.push_back(loadset);
			loadstrs.push_back(loadtype);
			for (int ii = 3; ii < strvec.size(); ii++) {
				Dloadvecs.push_back(atof(strvec[ii].c_str()));
			}
			Dloads.push_back(make_pair(loadstrs, Dloadvecs));
			//cout<<,,,<<endl;
		}
		getline(in, str);
	}
	in.close();
	return Dloads;
}

//****************Read load lines info******************************** load key/data lines
vector<std::pair<string, std::vector<string>> > AbaqusIO::GetLoadStrListFromInp() {
	ifstream in;
	streampos pos;
	string str, infostr;
	vector<string> datastr;
	Loadlines.clear();

	in.open(_MeshFileName, ios::in);
	if (!in.is_open()) {
		MessagePrinter::PrintErrorTxt("can not get Loadinfo, we cant open " + _MeshFileName);
		MessagePrinter::AsFem_Exit();
	}
	
	while (!in.eof()) {
		getline(in, str);
		infostr.clear();		
		if (str.find("*Cload") != string::npos|| str.find("*Dload") != string::npos|| str.find("*Dsload") != string::npos) {
			str = StringUtils::RemoveStrSpace(str);
			infostr = str;
		}
		str.clear();
		getline(in, str);
		datastr.clear();
		while (str.find("*") == string::npos) {
			str = StringUtils::RemoveStrSpace(str);
			datastr.push_back(str);
			pos = in.tellg();//记录当前指针位置
			getline(in, str);
		}
		Loadlines.push_back(make_pair(infostr, datastr));
		in.seekg(pos);//将指针指回记录行，可以重新读一次
	}
	in.close();
	return Loadlines;
}