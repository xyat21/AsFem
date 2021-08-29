//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <string>
#include <vector>

#include "Utils/Vector3d.h"
#include "LoadSystem/LoadType.h"
#include "AmpSystem/AmpType.h"

using namespace std;


class LoadBlock {//
public:
	LoadBlock() {
		_LoadBlockName.clear();
		_LoadkeyName.clear();
		_LoadTypeName.clear();
		_LoadType = LoadType::NULLLC;
		_LoadSetName.clear();
		_Partunit.nodesNo.clear();
		_Partunit.elesNo.clear();
		_Partunit.sDof.clear();
		_Partunit.values.clear();
		_DofsName.clear();
		//_LoadAmpName.clear();
		_LoadAmpName = "instantaneous";
		_LoadIndex=0;
		//_DofID = -1;
		_LoadValues.clear();
		_direction.setZero();
		_IsFollower = false;
		_IsConstArea = false;
		_IsOPnew = false;
	}

	//keyword info line
	int            _LoadIndex;//accumulated nums 
	string         _LoadBlockName;//"loadtype+i" Cload1 ...
	//LoadOP		   _OPType;//只有两种  不如用bool控制 但是可以没有控制项...
	bool		   _IsOPnew;//真则删除所有现有CLoads并重新定义新载荷，假则保留先前定义，然后修改或继续定义新CLoads
	string		   _LoadAmpName;
	bool           _IsFollower;
	bool           _IsConstArea;

	//obtained from data line
	string         _LoadSetName;
	struct PartUnit
	{
		vector<int> nodesNo;//如果不是setname，比如用函数AnalyticalField施加会分配计算 力大小
		vector<int> elesNo;
		vector<string> sDof;// 自由度分1-6，GRAV/体力分BX-BZ转换为1-3（注意仅能据此确定载荷类型）
		vector<double> values;

	}_Partunit;

	// judged from abrr. name
	string         _LoadkeyName;
	string         _LoadTypeName;
	LoadType       _LoadType;

	//string         _DofName;
	vector<string> _DofsName;//一句keyword下的数据有多行，如果是set[一句keyword下只可能有一个set] 则可以多行合并
	//int            _DofID;
	//double         _LoadValue;
	vector<double> _LoadValues;

	// other parameters
	Vector3d _direction;//即使是逐点模式（函数计算的）,各点体力方向均相同

	void Init() {
		_LoadBlockName.clear();
		_LoadkeyName.clear();
		_LoadTypeName.clear();
		_LoadType = LoadType::NULLLC;
		_LoadSetName.clear();
		_Partunit.nodesNo.clear();
		_Partunit.elesNo.clear();
		_Partunit.sDof.clear();
		_Partunit.values.clear();
		_DofsName.clear();
		//_LoadAmpName.clear();
		_LoadAmpName = "instantaneous";
		_LoadIndex = 0;
		//_DofID = -1;
		_LoadValues.clear();
		_direction.setZero();
		_IsFollower = false;
		_IsConstArea = false;
		_IsOPnew = false;
	}
};














