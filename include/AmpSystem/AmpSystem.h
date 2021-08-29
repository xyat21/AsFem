//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.10
//+++ Purpose: Define the amplitude curve system in VFIFPM
//+++          Here we can apply:
//+++               1)  ...
//+++               2) 
//+++               3)  as well as user defined amp (uamp)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>


//******************************************
//*** for VFIFPM own header
//******************************************
#include "Utils/MessagePrinter.h"

#include "AmpSystem/AmpBlock.h"
#include "AmpSystem/AmpType.h"
//
//#include "Mesh/Nodes.h"
//#include "Mesh/Mesh.h"
//// #include "DofHandler/DofHandler.h" // this line must be comment out to get rid of circular include issue from DofHandler class !!!
//#include "FE/FE.h"
//#include "FESystem/FECalcType.h"

#include "Utils/Vector3d.h"

using namespace std;

//class Mesh;
//class DofHandler;

class AmpSystem {
public:
	AmpSystem();
	void AddAmpBlock2List(AmpBlock &ampblock);

	//void InitAmpSystem(const Mesh &mesh);

	inline int GetAmpBlockNums()const { return _nAmpBlocks; }
	inline AmpBlock GetIthAmpBlock(const int &i)const { return _AmpBlockList[i - 1]; }

	//bool CheckAppliedAmpNameIsValid(const Mesh &mesh);
	//**************************************************************
	//*** add basic settings
	//**************************************************************
	void SetAmpScaleFactor(const double &factor) { _ScaleFactor = factor; }


	//**************************************************************
	//*** Compute amplitude value
	//**************************************************************
	void ComputeiAmp(const int &ith, const double &t, double &ampvalue);//计算第i个幅值曲线的t时刻幅值。
	void ComputeInitialAmp(const double &t, double &ampvalue);

	void PrintAmpSystemInfo()const;

private:
	//**************************************************************
	//*** some basic get functions
	//**************************************************************
	inline AmpType GetIthAmpBlockAmpType(const int &i)const { return _AmpBlockList[i - 1]._AmpType; }
	inline vector<string> GetIthAmpBlockNameVec(const int &i)const { return _AmpBlockList[i - 1]._AmpNameList; }
	inline vector<double> GetIthAmpBlockAmpParams(const int &i)const { return _AmpBlockList[i - 1]._AmpParams; }
	inline string GetIthAmpBlockName(const int &i)const { return _AmpBlockList[i - 1]._AmpBlockName; }
	int GetAmpIndexFromName(const string &ampname) {
		int index = 0;
		int tempNo= 0;
		for (auto it : _AmpBlockList) {
			index++;
			if (ampname == it._AmpBlockName) { tempNo=_AmpBlockList[index - 1]._AmpIndex; }
		}
		return tempNo;//补充 [如果没有匹配值] 则返回0，默认是const type类型
	}

	//**************************************************************
	//*** for different amplitude curves compute
	//**************************************************************
	void ComputeTabularAmp(const double &t, const double &smoothdata, const vector<double> &ampparams, double &ampvalue);
	void ComputeEqualSpaceAmp(const double &t, const double &smoothdata, const double &fixInterval, const vector<double> &ampparams, double &ampvalue);//NOTE:inp has two different ways to fill !!
	void ComputePeriodicAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputeModulatedAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputeDecayAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputeSolutionBasedAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputeSmoothAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputeSpectrumAmp(const double &t, const vector<double> &ampparams, double &ampvalue);
	void ComputePSDAmp(const double &t, const vector<double> &ampparams, double &ampvalue);

	//**************************************************************
	//*** for other general amplitude curves
	//**************************************************************


private:
	int _nAmpBlocks;
	vector<AmpBlock> _AmpBlockList;

private:
	double _ScaleFactor;
	//int _nBCDim, _nDim, _nBulkDim, _nNodesPerBCElmt;
	PetscMPIInt _rank, _size;

	////*******************************
	////*** for boundary integration
	////*******************************
	//double _xi, _eta, _JxW;
	//Nodes _elNodes;
	//double _xs[3][3], _dist;

	//Vector3d _normals, _gpGradU, _gpCoord;
	//double _gpU;
	//double _localR, _localK;
	//Vec _Useq;
	//VecScatter _scatteru;

};
















