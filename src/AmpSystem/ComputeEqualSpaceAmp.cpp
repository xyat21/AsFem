//****************************************************************
//* This file is part of the VFFPM framework
//* Vector Form Finite Particle Method program (VFIFPM)
//* All rights reserved, Xiaoyuan @ CopyRight 2021
//* https://github.com/XiaoyuanGao/LDFSI.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the Amplitude curve in VFFPM
//+++          here we define the Equal Space Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void LinearInterpolation(vector<double> &amppnts, const double &t, double &ampvalue) {
	//amppnts=(x1,y1,x2,y2)  range:x1<=t<=x2
	if ((amppnts[2] - amppnts[0]) != 0.0) {
		ampvalue = amppnts[1] + (amppnts[3] - amppnts[1]) / (amppnts[2] - amppnts[0])*(t - amppnts[0]);
	}
}

void AmpSystem::ComputeEqualSpaceAmp(const double &t, const double &smoothdata, const double &fixInterval, const vector<double> &ampparams, double &ampvalue) {
	//default tn=n*fixInterval;		ampparams[n]=Amp
	int pn;
	double tn,tn1;
	vector<double> twopoints(4);
	pn = ampparams.size();
	for (int n = 1; n <= pn-1; ++n) {
		tn = (n-1)*fixInterval;
		tn1 = n * fixInterval;
		if ((t >= tn) && (t < tn1)) {
			twopoints.clear();
			twopoints.push_back(tn);
			twopoints.push_back(ampparams[n-1]);
			twopoints.push_back(tn1);
			twopoints.push_back(ampparams[n]);
			LinearInterpolation(twopoints, t, ampvalue);//ampvalue can be calculated.
			break;
		}
		else if (t >= (pn-1) * fixInterval) {
			ampvalue = ampparams.back();
		}
	}

}

