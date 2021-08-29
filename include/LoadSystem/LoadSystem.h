//****************************************************************

//****************************************************************
#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

//******************************************
//*** for VFIFPM own header
//******************************************
#include "Utils/MessagePrinter.h"
#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"
//#include "DofHandler/DofHandler.h" // this line must be comment out to get rid of circular include issue from DofHandler class !!!
#include "FE/FE.h"
#include "FESystem/FECalcType.h"
#include "AmpSystem/AmpBlock.h"

#include "LoadSystem/LoadBlock.h"
#include "LoadSystem/LoadType.h"

using namespace std;

class Mesh;
class DofHandler;

class LoadSystem {
public:
	LoadSystem();
	void AddLoadBlock2List(LoadBlock &loadblock);
	void InitLoadSystem(const Mesh &mesh);

	inline int GetLoadBlockNums()const { return _nLoadBlocks; }
	inline LoadBlock GetIthLoadBlock(const int &i)const { return _LoadBlockList[i - 1]; }

	//**************************************************************
	//*** Apply load conditions
	//**************************************************************
	void ApplyLoad(const Mesh &mesh, const DofHandler &dofHandler, FE &fe, const FECalcType &calctype, const double &t, Vec &U, Vec &RHS);//根据需要改输入
	void ApplyInitialLoad(const Mesh &mesh, const DofHandler &dofHandler, const double &t, Vec &U);//根据需要改输入

	void PrintLoadSystemInfo()const;

private:
	//**************************************************************
	//*** some basic get functions
	//**************************************************************
	inline LoadType GetIthLoadBlockLoadType(const int &i)const { return _LoadBlockList[i - 1]._LoadType; }
	inline string GetIthLoadBlockSetName(const int &i)const { return _LoadBlockList[i - 1]._LoadSetName; }
	inline vector<double> GetIthLoadBlockLoadValue(const int &i)const { return _LoadBlockList[i - 1]._LoadValues; }

	//**************************************************************
	//*** for different nodal/elemental type load conditions
	//**************************************************************
	void ApplyCload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &RHS);
	void ApplyDload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &RHS);
	void ApplyDsload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &RHS);

	//**************************************************************
	//*** for different sub type load conditions
	void ApplyPressure(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &ampvalue, const LoadBlock &loadblock, Vec &RHS);
	void ApplySurfTraction(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &ampvalue, const LoadBlock &loadblock, Vec &RHS);


	//**************************************************************
	//*** for other general load conditions
	//**************************************************************
	void RunLoadLibs(const LoadType bctype, const FECalcType &calctype, const Vector3d &normals);
	void User1LC(const FECalcType &calctype, const Vector3d &normals);

private:
	int _nLoadBlocks;
	vector<LoadBlock> _LoadBlockList;

private:
	double _xi, _eta, _JxW;

	double _xs[3][3];
	Vector3d _normals;
	PetscMPIInt _rank, _size;

};

