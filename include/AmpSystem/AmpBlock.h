//****************************************************************

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <string>
#include <vector>


#include "AmpSystem/AmpType.h"

using namespace std;


class AmpBlock {
public:
	AmpBlock() {
		_AmpBlockName.clear();
		_AmpTypeName.clear();
		_AmpType = AmpType::CONSTAMP;
		_AmpTab.clear();
		_AmpParams.clear();
		_AmpIndex=0;
		_Smoothdata = 0.0;
		_FixInterval = 0.0;
	}

	int			   _AmpIndex;
	string         _AmpBlockName;//name=
	string         _AmpTypeName;//amptype as string "..."
	AmpType        _AmpType;//DEFINITION=

	vector<std::pair<double, double> > _AmpTab;//tabular (time,ampvalue)
	vector<double> _AmpParams;//不用上面格式，每种需要定义不同格式，不如都存储然后根据type去识别所需参数格式
	//double ComputeAt(const AmpType &amptype, const double &dt, vector<std::pair<double, double> > &amptab);
	//double ComputeAt(const AmpType &amptype, const double &dt);
	vector<string> _AmpNameList;// it could be either an element set or a node set  没用
	double		   _Smoothdata;
	double		   _FixInterval;

	void Init() {
		_AmpBlockName.clear();
		_AmpTypeName.clear();
		_AmpType = AmpType::CONSTAMP;
		_AmpTab.clear();
		_AmpParams.clear();
		_AmpIndex=0;
		_Smoothdata=0.0;
		_FixInterval=0.0;
	}

};
