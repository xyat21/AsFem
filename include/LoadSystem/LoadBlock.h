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
	//LoadOP		   _OPType;//ֻ������  ������bool���� ���ǿ���û�п�����...
	bool		   _IsOPnew;//����ɾ����������CLoads�����¶������غɣ���������ǰ���壬Ȼ���޸Ļ����������CLoads
	string		   _LoadAmpName;
	bool           _IsFollower;
	bool           _IsConstArea;

	//obtained from data line
	string         _LoadSetName;
	struct PartUnit
	{
		vector<int> nodesNo;//�������setname�������ú���AnalyticalFieldʩ�ӻ������� ����С
		vector<int> elesNo;
		vector<string> sDof;// ���ɶȷ�1-6��GRAV/������BX-BZת��Ϊ1-3��ע����ܾݴ�ȷ���غ����ͣ�
		vector<double> values;

	}_Partunit;

	// judged from abrr. name
	string         _LoadkeyName;
	string         _LoadTypeName;
	LoadType       _LoadType;

	//string         _DofName;
	vector<string> _DofsName;//һ��keyword�µ������ж��У������set[һ��keyword��ֻ������һ��set] ����Զ��кϲ�
	//int            _DofID;
	//double         _LoadValue;
	vector<double> _LoadValues;

	// other parameters
	Vector3d _direction;//��ʹ�����ģʽ����������ģ�,���������������ͬ

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














