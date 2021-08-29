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
//+++ Purpose: implement Abaqus' inp mesh file importer
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.07.21  add more abaqusIO funs concerning other inp keylines

#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include "Mesh/MeshIOBase.h"
#include "Mesh/MeshType.h"

#include "Utils/MessagePrinter.h"
#include "Utils/StringUtils.h"

class AbaqusIO:public MeshIOBase{
public:
    AbaqusIO(){
        _MeshFileName.clear();
        _HasSetMeshFileName=false;
    }

    virtual bool ReadMeshFromFile(Mesh &mesh) override;
    virtual void SetMeshFileName(string filename)override{_MeshFileName=filename;_HasSetMeshFileName=true;}
    virtual string GetMeshFileName()const override{return _MeshFileName;}

	// for Amp/boundary/load information   =====XY======
	//inline string GetAbaqusFileName()const { return _MeshFileName; }   //
	//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > GetAmpsTimeFromInp();
	vector<std::pair<std::vector<std::string>, std::vector<double>> >  GetAmpsTimeFromInp();
	vector<std::pair<string, vector<int>> > GetBoundaryDofFromInp();
	vector<std::pair<pair<string, string>, vector<double>> > GetCLoadFromInp();
	vector<std::pair<vector<string>, double> > GetDsLoadFromInp();
	vector<std::pair<vector<string>, vector<double>> > GetDloadsFromInp();
	vector<std::pair<string, std::vector<string>> > GetLoadStrListFromInp();
	struct material_info {
		string matname;
		double density;
		double Youngs_modulus;
		double poisson;
		double damp_alpha = 0.0;
		double damp_beta = 0.0;
	};//Mate;
	vector<material_info> Mate;
	vector<material_info> GetMateFromInp();
	vector<std::pair<pair<string, string>, vector<double>> > Cloads;
	vector<std::pair<vector<string>, double> > Dsloads;//string:amp,setname,type   double:value,directions
	vector<std::pair<vector<string>, vector<double>> > Dloads;//string:amp,setname,type   double:value,directions
	vector<std::pair<string, std::vector<string>> > Loadlines;

	int GetBeamSectionNumFromInp() const;
	int GetSolidSectionNumFromInp() const;
	int GetShellSectionNumFromInp() const;
	vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> > GetBeamSectionESetsandMatAdFromInp() const;

	vector<pair<pair<string, string>, double>> GetSolidSectionESetsandMatAFromInp() const;
	vector<pair<pair<string, string>, double>> GetShellSectionESetsandMatTFromInp() const;

	struct step_info {
		//vector<int> static_run;
		//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > AmpsTime;	//amp name/type/value-table
		vector<std::pair<std::vector<std::string>, std::vector<double>> > AmpsTime;
		vector<std::pair<string, vector<int>> > boundaryDoF;
		//vector<std::pair<string, vector<double>> > CLoad;
		//vector<int> output;
	}Step;		//History


private:
    int GetElmtNodesNumFromInp() const;
    int GetElmtDimFromInp() const;
	int GetSpaceDimFromInp() const;
	int GetElmtOrderFromInp()const;
    int GetElmtVTKCellTypeFromInp()const;
    MeshType GetElmtMeshTypeFromInp()const;
    string GetElmtMeshTypeNameFromInp()const;

    int GetSubElmtNodesNumFromInp() const;
	int GetSubElmtDimFromInp() const;//原代码多了FromGmshInp
    int GetSubElmtOrderFromInpElmt()const;
    MeshType GetSubElmtMeshTypeFromInp()const;

    int GetSubSubElmtNodesNumFromInp() const;
    MeshType GetSubSubElmtMeshTypeFromInp()const;


    int GetNodesNumFromInp()const;
    int GetElmtsNumFromInp()const;
    int GetNsetsNumFromInp()const;
    int GetElsetsNumFromInp()const;
    int GetSurfacesNumFromInp()const;
    int GetSurfaceElmtsNumFromInp()const;
    int GetSurfaceEdgeIDViaSurfaceNameFromInp(string surfacesetname)const;
    vector<int> GetSurfaceElmtIDViaSurfaceNameFromInp(string surfacesetname)const;
	//vector<int> GetElmtIDViaElSetNameFromInp(string elsetname)const;
	//vector<int> GetNodeIDViaNSetNameFromInp(string nsetname)const;


    ifstream _in;

    int _nMaxDim=-1,_nMinDim=4;//mesh dimension
	int _nSpaceDim = -1;	   //space dimension
    int _nPhysicGroups=0;
    int _nNodeSetPhysicalGroups=0;
    int _nNodes=0,_nElmts=0;
    int _nBulkElmts=0,_nSurfaceElmts=0,_nLineElmts=0;

    double _Xmax,_Xmin,_Ymax,_Ymin,_Zmax,_Zmin;

    int _nNodesPerBulkElmt=-1;
    int _nNodesPerLineElmt=0;
    int _nNodesPerSurfaceElmt=0;
    int _nOrder=1;

	int _BeamSectionNum;
	int _SolidSectionNum;
	int _ShellSectionNum;

};