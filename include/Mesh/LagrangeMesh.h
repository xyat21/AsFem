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
//+++ Date   : 2020.06.28
//+++ Purpose: Define bulk mesh class for the bulk domain(except 
//+++          the interface domain)
//+++          This class offers the functions to get as well as
//+++          create the mesh and the related geometric info
//+++          In order to import mesh from external file, one should
//+++          call the meshio class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.07.21  add funs to read more inp keylines info.

#pragma once

#include<iostream>
#include<iomanip>
#include<vector>
#include <array>
#include<cstdio>
#include<cmath>
#include<set>
#include<numeric>
#include<algorithm>
#include<fstream>

#include "petsc.h"
#include "Utils/MessagePrinter.h"

#include "Mesh/MeshType.h"
#include "Mesh/Nodes.h"

using namespace std;


class LagrangeMesh{
public:
    LagrangeMesh();

    bool CreateLagrangeMesh();
    void SaveLagrangeMesh(string inputfilename="") const;
    //************************************************************
    //*** for the basic settings
    //************************************************************
    void SetBulkMeshDim(const int &ndim) {_nMaxDim=ndim;_nMinDim=ndim-1;}
    void SetBulkMeshMaxDim(const int &ndim) {_nMaxDim=ndim;}
    void SetBulkMeshMinDim(const int &ndim) {_nMinDim=ndim;}
	void SetSpaceDim(const int &nspacedim) { _nSpaceDim = nspacedim; }//
	void SetBulkMeshMeshOrder(const int &p){_nOrder=p;}
    //*** for element numbers 
    void SetBulkMeshNx(const int &nx){_Nx=nx;}
    void SetBulkMeshNy(const int &ny){_Ny=ny;}
    void SetBulkMeshNz(const int &nz){_Nz=nz;}
    void SetBulkMeshNodesNum(const int &n){_nNodes=n;}
    void SetBulkMeshNodesNumPerBulkElmt(const int &n){_nNodesPerBulkElmt=n;}
    void SetBulkMeshNodesNumPerSurfaceElmt(const int &n){_nNodesPerSurfaceElmt=n;}
    void SetBulkMeshNodesNumPerLineElmt(const int &n){_nNodesPerLineElmt=n;}
    void SetBulkMeshElmtsNum(const int &n){_nElmts=n;}
    void SetBulkMeshBulkElmtsNum(const int &n){_nBulkElmts=n;}
    void SetBulkMeshSurfaceElmtsNum(const int &n){_nSurfaceElmts=n;}
    void SetBulkMeshLineElmtsNum(const int &n){_nLineElmts=n;}
    //*** for element geometry settings
    void SetBulkMeshXmin(const double &xmin){_Xmin=xmin;}
    void SetBulkMeshXmax(const double &xmax){_Xmax=xmax;}
    void SetBulkMeshYmin(const double &ymin){_Ymin=ymin;}
    void SetBulkMeshYmax(const double &ymax){_Ymax=ymax;}
    void SetBulkMeshZmin(const double &zmin){_Zmin=zmin;}
    void SetBulkMeshZmax(const double &zmax){_Zmax=zmax;}
    //*** for mesh type setting
    void SetBulkMeshMeshTypeName(string meshname);
    void SetBulkMeshMeshType(const MeshType &type){_BulkMeshType=type;}
    void SetBulkMeshSurfaceMeshType(const MeshType &type){_SurfaceMeshType=type;}
    void SetBulkMeshLineMeshType(const MeshType &type){_LineMeshType=type;}
    //*** for elmt volume settings
    void SetBulkMeshIthElmtVolume(const int &i,const double &volume){_ElmtVolume[i-1]=volume;}
	void SetBulkMeshIthBulkElmtVolume(const int &i,const double &volume){_ElmtVolume[i+_nElmts-_nBulkElmts-1]=volume;}//�����subelmt,�����ǰ�ߺ�bulkelmt.
	//void SetBulkMeshIthBulkElmtMass(const int &i, const double &mass) {_ElmtMass[i+_nElmts- _nBulkElmts -1]=mass;}//����ֱ���ýڵ��Ч������ɽڵ���������
	void SetModelIthNodeeqVolume(const int &i, const double &volume) { _NodeeqVolume[i - 1] = volume; }
	void SetTotalVolume(const double &volume) { _TotalVolume = volume; }

	//*** for physical group settings
    void SetBulkMeshPhysicalGroupNums(const int &n){_nPhysicalGroups=n;}
    void SetBulkMeshNodeSetPhysicalGroupNums(const int &n){_nNodeSetPhysicalGroups=n;}

	void SetBeamSectionNum(const int &n) { _BeamSectionNum = n; }
	void SetSolidSectionNum(const int &n) { _SolidSectionNum = n; }
	void SetShellSectionNum(const int &n) { _ShellSectionNum = n; }
	//************************************************************
    //*** for advanced getting functions (allow value modify externally!)
    //************************************************************
    vector<double>&      GetBulkMeshNodeCoordsPtr(){return _NodeCoords;}
    vector<vector<int>>& GetBulkMeshElmtConnPtr(){return _ElmtConn;}
    //vector<double>&      GetBulkMeshElmtVolumePtr(){return _ElmtVolume;}
	vector<double>&      GetBulkMeshElmtVolumeListPtr() { return _ElmtVolume; }
	vector<double>&      GetBulkMeshNodeeqVolumeListPtr() { return _NodeeqVolume; }

    vector<int>&         GetBulkMeshElmtVTKCellTypeListPtr(){return _ElmtVTKCellTypeList;}
    vector<int>&         GetBulkMeshElmtPhyIDListPtr(){return _ElmtPhyIDList;}
    vector<int>&         GetBulkMeshElmtDimListPtr(){return _ElmtDimList;}
    vector<MeshType>&    GetBulkMeshElmtMeshTypeListPtr(){return _ElmtMeshTypeList;}
    //*** for physical groups
    vector<string>&                   GetBulkMeshPhysicalGroupNameListPtr(){return _PhysicalGroupNameList;}
    vector<int>&                      GetBulkMeshPhysicalGroupIDListPtr(){return _PhysicalGroupIDList;}
    vector<int>&                      GetBulkMeshPhysicalGroupDimListPtr(){return _PhysicalGroupDimList;}
    vector<pair<string,int>>&         GetBulkMeshPhysicalGroupName2DimListPrt(){return _PhysicalGroupName2DimList;}
    vector<pair<int,string>>&         GetBulkMeshPhysicalGroupID2NameListPtr(){return _PhysicalGroupID2NameList;}
    vector<pair<string,int>>&         GetBulkMeshPhysicalGroupName2IDListPtr(){return _PhysicalGroupName2IDList;}
    vector<pair<string,int>>&         GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr(){return _PhysicalGroupName2NodesNumPerElmtList;}
    vector<pair<string,vector<int>>>& GetBulkMeshPhysicalName2ElmtIDsListPtr(){return _PhysicalName2ElmtIDsList;}
    // for node set physical groups
    vector<string>&                   GetBulkMeshNodeSetPhysicalNameListPtr(){return _NodeSetPhysicalGroupNameList;}
    vector<int>&                      GetBulkMeshNodeSetPhysicalIDListPtr(){return _NodeSetPhysicalGroupIDList;}
    vector<pair<int,string>>&         GetBulkMeshNodeSetPhysicalID2NameListPtr(){return _NodeSetPhysicalGroupID2NameList;}
    vector<pair<string,int>>&         GetBulkMeshNodeSetPhysicalName2IDListPtr(){return _NodeSetPhysicalGroupName2IDList;}
    vector<pair<string,vector<int>>>& GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr(){return _NodeSetPhysicalName2NodeIDsList;}

	// =====XY======
	vector<std::pair<pair<string, string>, vector<double>> >& GetCLoadListPtr() { return _CLoadList; };
	vector<std::pair<vector<string>, double> >& GetDsLoadListPtr() { return _DsLoadList; };
	vector<std::pair<vector<string>, vector<double>> >& GetDLoadsListPtr() { return _DLoadsList; };
	vector<std::pair<string, std::vector<string>> >& GetLoadInfoListPtr() { return _LoadStrList; };
	struct material_info {
		string matname;
		double density;
		double Youngs_modulus;
		double poisson;
		double damp_alpha = 0.0;
		double damp_beta = 0.0;
	};//Mate;
	struct step_info {
		//vector<int> static_run;
		//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > _AmpsTimeList;
		vector<std::pair<std::vector<std::string>, std::vector<double>> > _AmpsTimeList;
		vector<std::pair<string, vector<int>> > _BoundaryDoFList;
		//vector<std::pair<string, vector<double>> > CLoad;
		//vector<int> output;
	}Step;		//History
	vector<material_info>& GetMateListPtr() { return _MateList; };
	//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> >& GetAmpsTimeListPtr() { return Step._AmpsTimeList; };
	vector<std::pair<std::vector<std::string>, std::vector<double>> >& GetAmpsTimeListPtr() { return Step._AmpsTimeList; };
	vector<std::pair<string, vector<int>> >& GetBoundaryDoFListPtr() { return Step._BoundaryDoFList; };
	//vector<pair<string, pair<double, double>> >& GetBeamSectionESetsandMatAdListPtr() { return _BeamSectionESetsandMatAdList; };
	vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> >& GetBeamSectionESetsandMatAdListPtr() { return _BeamSectionESetsandMatAdList; };
	vector<pair<pair<string, string>, double>>& GetSolidSectionESetsandMatAListPtr() { return _SolidSectionESetsandMatAList; };
	vector<pair<pair<string, string>, double>>& GetShellSectionESetsandMatTListPtr() { return _ShellSectionESetsandMatTList; };

	inline int GetBeamSectionNum()const { return _BeamSectionNum; }
	inline int GetSolidSectionNum()const { return _SolidSectionNum; }
	inline int GetShellSectionNum()const { return _ShellSectionNum; }

    //************************************************************
    //*** for the basic getting functions
    //************************************************************
    //*** for dimention of the bulk mesh
    inline int GetBulkMeshDim()const{return _nMaxDim;}
    inline int GetBulkMeshMinDim()const{return _nMinDim;}
	inline int GetSpaceDim()const { return _nSpaceDim; }
	//*** for nodes num of the bulk mesh
    inline int GetBulkMeshNodesNum()const{return _nNodes;}
    inline int GetBulkMeshNodesNumPerBulkElmt()const{return _nNodesPerBulkElmt;}
    inline int GetBulkMeshNodesNumPerSurfaceElmt()const{return _nNodesPerSurfaceElmt;}
    inline int GetBulkMeshNodesNumPerLineElmt()const{return _nNodesPerLineElmt;}
    //*** for elmts num of the bulk mesh
    inline int GetBulkMeshNx()const{return _Nx;}
    inline int GetBulkMeshNy()const{return _Ny;}
    inline int GetBulkMeshNz()const{return _Nz;}
    inline int GetBulkMeshElmtsNum()const{return _nElmts;}
    inline int GetBulkMeshBulkElmtsNum()const{return _nBulkElmts;}
    inline int GetBulkMeshSurfaceElmtsNum()const{return _nSurfaceElmts;}
    inline int GetBulkMeshLineElmtsNum()const{return _nLineElmts;}

	//*** for model total volume
	inline double GetTotalVolume()const { return _TotalVolume; }
	//*** for mesh order
    inline int GetBulkMeshOrder()const{return _nOrder;}
    //*** for mesh type
    inline MeshType GetBulkMeshBulkElmtType()const{return _BulkMeshType;}
    inline MeshType GetBulkMeshSurfaceElmtType()const{return _SurfaceMeshType;}
    inline MeshType GetBulkMeshLineElmtType()const{return _LineMeshType;}
    inline string   GetBulkMeshBulkElmtTypeName()const{return _BulkMeshTypeName;}
    //*** for node's coordinate
    inline double GetBulkMeshIthNodeJthCoord(const int &i,const int &j)const{return _NodeCoords[(i-1)*3+j-1];}
    inline void   GetBulkMeshIthNodeCoords(const int &i,double (&coords)[4])const{
        coords[0]=0.0;
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
    }
    inline void   GetBulkMeshIthNodeCoords(const int &i,vector<double> &coords)const{
        coords[0]=0.0;
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
    }
    inline vector<double> GetBulkMeshIthNodeCoords(const int &i)const{
        vector<double> coords(4,0.0);
        coords[0]=0.0;
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
        return coords;
    }
    inline void GetBulkMeshIthElmtNodes(const int &e,Nodes &nodes)const{
        for(int i=1;i<=GetBulkMeshIthElmtNodesNum(e);i++){
            nodes(i,0)=0.0;
            nodes(i,1)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthElmtJthNodeID(e,i),1);
            nodes(i,2)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthElmtJthNodeID(e,i),2);
            nodes(i,3)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthElmtJthNodeID(e,i),3);
        }
    }
    inline void GetBulkMeshIthBulkElmtNodes(const int &e,Nodes &nodes)const{
        for(int i=1;i<=GetBulkMeshIthBulkElmtNodesNum(e);i++){
            nodes(i,0)=0.0;
            nodes(i,1)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthBulkElmtJthNodeID(e,i),1);
            nodes(i,2)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthBulkElmtJthNodeID(e,i),2);
            nodes(i,3)=GetBulkMeshIthNodeJthCoord(GetBulkMeshIthBulkElmtJthNodeID(e,i),3);
        }
    }
    //*** for elmt connectivity information
    inline int  GetBulkMeshIthElmtDim(const int &i)const{return _ElmtDimList[i-1];}
    inline int  GetBulkMeshIthElmtPhyID(const int &i)const{return _ElmtPhyIDList[i-1];}

    inline int  GetBulkMeshIthElmtVTKCellType(const int &i)const{return _ElmtVTKCellTypeList[i-1];}
    inline int  GetBulkMeshIthBulkElmtVTKCellType(const int &i)const{return _ElmtVTKCellTypeList[i+_nElmts-_nBulkElmts-1];}
    
    inline int  GetBulkMeshIthElmtNodesNum(const int &i)const{return _ElmtConn[i-1][0];}//��i��Ԫ�Ľڵ���
    inline int  GetBulkMeshIthBulkElmtNodesNum(const int &i)const{return _ElmtConn[i+_nElmts-_nBulkElmts-1][0];}
    
    inline int  GetBulkMeshIthElmtJthNodeID(const int &i,const int &j)const{return _ElmtConn[i-1][j];}//���ص�i����Ԫ��j�ڵ�
    inline int  GetBulkMeshIthBulkElmtJthNodeID(const int &i,const int &j)const{return _ElmtConn[i+_nElmts-_nBulkElmts-1][j];}//���ص�i��bulk��Ԫ��j�ڵ��
    inline void GetBulkMeshIthElmtNodeIDs(const int &i,vector<int> &elConn)const{
        for(int j=1;j<=_ElmtConn[i-1][0];j++) elConn[j-1]=_ElmtConn[i-1][j];
    }//�����i����Ԫ�����ӽڵ��vector<int>��elConn[0��������-1]
    inline vector<int> GetBulkMeshIthElmtNodesID(const int &i)const{
        vector<int> conn(_ElmtConn[i-1][0],0);
        for(int j=1;j<=_ElmtConn[i-1][0];j++) conn[j-1]=_ElmtConn[i-1][j];
        return conn;
    }//���ص�i����Ԫ�����ӽڵ��vector<int>
    inline int GetBulkMeshIthElmtIDViaPhyName(string phyname,const int &id)const{
        for(const auto &it:_PhysicalName2ElmtIDsList){
            if(it.first==phyname){
                return it.second[id-1];
            }
        }
        return -1;
    }
    inline int GetBulkMeshIthElmtNodesNumViaPhyName(const string phyname,const int &id)const{//ͨ���������ֻ�ȡ��i��Ԫ�����ӽڵ���
        int e;
        for(const auto &it:_PhysicalName2ElmtIDsList){
            if(it.first==phyname){
                e=it.second[id-1];
                return _ElmtConn[e-1][0];
            }
        }
        return -1;
    }
    inline vector<int> GetBulkMeshIthElmtNodeIDs(const int &i)const{//���ص�i����Ԫ�����ӽڵ��vector<int>
        vector<int> temp(_ElmtConn[i-1][0],0);
        for(int j=1;j<=_ElmtConn[i-1][0];j++) temp[j-1]=_ElmtConn[i-1][j];
        return temp;
    }
    inline void GetBulkMeshIthBulkElmtConn(const int &e,vector<int> &conn)const{//�����e��bulk��Ԫ�����ӽڵ��vector<int>��conn[0��������-1]
        for(int i=1;i<=GetBulkMeshIthBulkElmtNodesNum(e);i++){
            conn[i-1]=GetBulkMeshIthBulkElmtJthNodeID(e,i);
        }
    }
    //*** for physical group information
    inline int GetBulkMeshNodeSetPhysicalGroupNum()const{return _nNodeSetPhysicalGroups;}
    inline int GetBulkMeshPhysicalGroupNum()const{return _nPhysicalGroups;}
    inline int GetBulkMeshIthPhysicalID(const int &i)const{return _PhysicalGroupIDList[i-1];}
    inline int GetBulkMeshIthPhysicalDim(const int &i)const{return _PhysicalGroupDimList[i-1];}
    inline string GetBulkMeshIthPhysicalName(const int &i)const{return _PhysicalGroupNameList[i-1];}
    inline int GetBulkMeshDimViaPhyName(string phyname)const{
        for(int i=0;i<static_cast<int>(_PhysicalGroupNameList.size());i++){
            if(_PhysicalGroupNameList[i]==phyname){
                return _PhysicalGroupDimList[i];
            }
        }
        return -1;
    }
    inline string GetBulkMeshPhysicalNameViaPhyID(const int &phyid)const{
        for(const auto &it:_PhysicalGroupID2NameList){
            if(it.first==phyid){
                return it.second;
            }
        }
        return "";
    }
    inline int GetBulkMeshPhysicalIDViaName(string phyname)const{
        for(const auto &it:_PhysicalGroupName2IDList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return -1;
    }
    inline int GetBulkMeshNodesNumPerElmtViaPhysicalName(string phyname)const{
        for(const auto &it:_PhysicalGroupName2NodesNumPerElmtList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return -1;
    }
    inline vector<int> GetBulkMeshElmtIDsViaPhysicalName(string phyname)const{
        for(const auto &it:_PhysicalName2ElmtIDsList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return vector<int>(0);
    }
    inline int GetBulkMeshElmtsNumViaPhysicalName(string phyname)const{
        for(const auto &it:_PhysicalName2ElmtIDsList){
            if(it.first==phyname){
                return static_cast<int>(it.second.size());
            }
        }
        return 0;
    }
    inline vector<int> GetBulkMeshNodeIDsViaPhysicalName(string phyname)const{
        for(const auto &it:_NodeSetPhysicalName2NodeIDsList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return vector<int>(0);
    }
    inline int GetBulkMeshNodeIDsNumViaPhysicalName(string phyname)const{
        for(const auto &it:_NodeSetPhysicalName2NodeIDsList){
            if(it.first==phyname){
                return static_cast<int>(it.second.size());
            }
        }
        return 0;
    }
    inline int GetBulkMeshIthNodeIDViaPhyName(const string phyname,const int &id)const{
        for(const auto &it:_NodeSetPhysicalName2NodeIDsList){
            if(it.first==phyname){
                return it.second[id-1];
            }
        }
        return -1;
    }
    inline int GetBulkMeshIthNodeSetPhysicalID(const int &i)const{
        return _NodeSetPhysicalGroupIDList[i-1];
    }
    inline string GetBulkMeshIthNodeSetPhysicalName(const int &i)const{
        return _NodeSetPhysicalGroupNameList[i-1];
    }


    //************************************************************
    //*** for mesh information printer
    //************************************************************
    void PrintBulkMeshInfo()const;
    void PrintBulkMeshInfoDetails()const;

private:
    bool Create1DLagrangeMesh();
    bool Create2DLagrangeMesh();
    bool Create3DLagrangeMesh();

protected:
    //************************************************************
    //*** for the basic information of our mesh
    //************************************************************
    bool _IsMeshCreated=false;
    int _nNodes,_nElmts;
    int _nNodesPerBulkElmt,_nNodesPerSurfaceElmt,_nNodesPerLineElmt;
    int _nBulkElmts,_nSurfaceElmts,_nLineElmts;
    int _nMaxDim,_nMinDim, _nSpaceDim;
    int _Nx,_Ny,_Nz;
    double _Xmax,_Xmin,_Ymax,_Ymin,_Zmin,_Zmax;
    int _nOrder;
    MeshType _BulkMeshType,_SurfaceMeshType,_LineMeshType;
    vector<double>      _NodeCoords;// store all the nodes/controlpts' coordinate
    vector<vector<int>> _ElmtConn;  // store all the element(bulk+surface+line+node elements)
    vector<double>      _ElmtVolume;// it could be: volume(3D), area(2D) or length(1D)
	vector<double>      _NodeeqVolume;//=_ElmtVolume/NodesnumPerelmt and accumulate

    vector<int>         _ElmtVTKCellTypeList;
    vector<int>         _ElmtPhyIDList;
    vector<int>         _ElmtDimList;
    vector<MeshType>    _ElmtMeshTypeList;
    int                 _BulkElmtVTKCellType;
    string              _BulkMeshTypeName;
    double              _TotalVolume;

    //************************************************************
    //*** for the basic physical group information
    //************************************************************
    int _nPhysicalGroups;
    vector<string>                   _PhysicalGroupNameList;
    vector<int>                      _PhysicalGroupIDList;
    vector<int>                      _PhysicalGroupDimList;
    vector<pair<string,int>>         _PhysicalGroupName2DimList;
    vector<pair<int,string>>         _PhysicalGroupID2NameList;
    vector<pair<string,int>>         _PhysicalGroupName2IDList;
    vector<pair<string,int>>         _PhysicalGroupName2NodesNumPerElmtList;
    vector<pair<string,vector<int>>> _PhysicalName2ElmtIDsList;
    //*** for nodal physical group
    int _nNodeSetPhysicalGroups;
    vector<string>                   _NodeSetPhysicalGroupNameList;
    vector<int>                      _NodeSetPhysicalGroupIDList;
    vector<pair<int,string>>         _NodeSetPhysicalGroupID2NameList;
    vector<pair<string,int>>         _NodeSetPhysicalGroupName2IDList;
    vector<pair<string,vector<int>>> _NodeSetPhysicalName2NodeIDsList;

	//************************************************************
	//*** for the Amp/boundary/load information   =====XY======
	//************************************************************
	//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > _AmpsTimeList;//amp name/type/value-table
	//vector<std::pair<string, vector<int>> > _BoundaryDofList;
	vector<std::pair<pair<string, string>, vector<double>> > _CLoadList;
	vector<std::pair<vector<string>, double> > _DsLoadList;//string:amp,setname,type   double:value,directions
	vector<std::pair<vector<string>, vector<double>> > _DLoadsList;//string:amp,setname,type   double:value,directions
	vector<std::pair<string, std::vector<string>> > _LoadStrList;
	//struct material_info {
	//	string matname;
	//	double density;
	//	double Youngs_modulus;
	//	double poisson;
	//	double damp_alpha = 0.0;
	//	double damp_beta = 0.0;
	//};//Mate;
	vector<material_info> _MateList;

	int _BeamSectionNum;
	int _SolidSectionNum;
	int _ShellSectionNum;
	//vector<pair<string, pair<double, double>> > _BeamSectionESetsandMatAdList;
	vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> > _BeamSectionESetsandMatAdList;

	vector<pair<pair<string, string>, double>> _SolidSectionESetsandMatAList;
	vector<pair<pair<string, string>, double>> _ShellSectionESetsandMatTList;

	//struct step_info {
	//	//vector<int> static_run;
	//	vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > _AmpsTimeList;	
	//	vector<std::pair<string, vector<int>> > _BoundaryDoFList;
	//	//vector<std::pair<string, vector<double>> > CLoad;
	//	//vector<int> output;
	//}Step;		//History

};