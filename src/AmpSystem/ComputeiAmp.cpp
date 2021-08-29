//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//#include <algorithm>    // std::copy
#include "AmpSystem/AmpSystem.h"

void AmpSystem::ComputeiAmp(const int &ith, const double &t, double &ampvalue) {

	//vector<string> ampnamelist;
	vector<double> ampparams;
	double		   smoothdata;
	double		   fixInterval;

	AmpBlock it;
	if (ith==0) {
		it._AmpParams.push_back(0.0);
		it._AmpParams.push_back(1.0);
		it._AmpParams.push_back(1.0);
		it._AmpParams.push_back(1.0);
		it._AmpBlockName = "instantaneous";
		it._AmpTypeName  = "CONSTAMP";
		it._AmpType = AmpType::CONSTAMP;
	}
	else {
		it = _AmpBlockList[ith - 1];
	}	

	//for (auto it : _AmpBlockList) {
		//ampnamelist = it._AmpNameList;
		ampparams = it._AmpParams;
		smoothdata= it._Smoothdata;
		fixInterval= it._FixInterval;
		if (it._AmpType == AmpType::TABULAR) {
			ComputeTabularAmp(t, smoothdata, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::EQUALSPACED) {
			ComputeEqualSpaceAmp(t, smoothdata, fixInterval, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::PERIODIC) {
			ComputePeriodicAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::MODULATED) {
			ComputeModulatedAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::DECAY) {
			ComputeDecayAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::SOLUTIONBASED) {
			ComputeSolutionBasedAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::SMOOTHSTEP) {
			ComputeSmoothAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::SPECTRUM) {
			ComputeSpectrumAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::PSD) {
			ComputePSDAmp(t, ampparams, ampvalue);
		}
		else if (it._AmpType == AmpType::CONSTAMP) { 
			ampvalue = 1.0;
			continue;
		}
		else {
			//User Amp
			MessagePrinter::PrintErrorTxt("unsupported Amplitude type in current AmpType, please check your code");
			MessagePrinter::AsFem_Exit();
		}

	//}
}

void AmpSystem::ComputeSolutionBasedAmp(const double &t, const vector<double> &ampparams, double &ampvalue) { continue; }
void AmpSystem::ComputeSpectrumAmp(const double &t, const vector<double> &ampparams, double &ampvalue) { continue; }
void AmpSystem::ComputePSDAmp(const double &t, const vector<double> &ampparams, double &ampvalue) { continue; }

//void LinearInterpolation(vector<double> &amppnts, const double &t, double &ampvalue) {
//	//amppnts=(x1,y1,x2,y2)  range:x1<=t<=x2
//	if ((amppnts[2] - amppnts[0]) != 0.0) {
//		ampvalue = amppnts[1] + (amppnts[3] - amppnts[1]) / (amppnts[2] - amppnts[0])*(t - amppnts[0]);
//	}
//}
//
//void AmpSystem::ComputeTabularAmp(const double &t, const double &smoothdata, const vector<double> &ampparams, double &ampvalue) {
//	//ampparams[0+2n]=time;	ampparams[1+2n]=Amp
//	int polyn;
//	vector<double> twopoints(4);
//	polyn= ampparams.size()/2-1;
//	for (int n = 1; n <= polyn; ++n) {
//		if ((t>= ampparams[2*(n-1)])&&(t < ampparams[2*n])) {
//			twopoints.clear();
//			std::copy(ampparams.begin()+ 2 * (n - 1), ampparams.begin() + 2 * (n + 1), twopoints.begin()); //Çø¼ä¿½±´
//			LinearInterpolation(twopoints, t, ampvalue);//ampvalue can be calculated.
//			break;
//		}	
//		else if (t >= ampparams[2*polyn]) {//*(ampparams.end()-2)
//			ampvalue = ampparams.back();
//		}
//	}
//
//}


