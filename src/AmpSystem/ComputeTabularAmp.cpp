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
//+++          here we define the Tabular Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void LinearInterpolation(vector<double> &amppnts, const double &t, double &ampvalue) {
	//amppnts=(x1,y1,x2,y2)  range:x1<=t<=x2
	if ((amppnts[2] - amppnts[0]) != 0.0) {
		ampvalue = amppnts[1] + (amppnts[3] - amppnts[1]) / (amppnts[2] - amppnts[0])*(t - amppnts[0]);
	}
}

void AmpSystem::ComputeTabularAmp(const double &t, const double &smoothdata, const vector<double> &ampparams, double &ampvalue) {
	//ampparams[0+2n]=time;	ampparams[1+2n]=Amp
	int polyn;
	vector<double> twopoints(4);
	polyn = ampparams.size() / 2 - 1;
	ampvalue = 0.0;
	for (int n = 1; n <= polyn; ++n) {
		if ((t >= ampparams[2 * (n - 1)]) && (t < ampparams[2 * n])) {
			twopoints.clear();
			std::copy(ampparams.begin() + 2 * (n - 1), ampparams.begin() + 2 * (n + 1), twopoints.begin()); //Çø¼ä¿½±´
			LinearInterpolation(twopoints, t, ampvalue);//ampvalue can be calculated.
			break;
		}
		else if (t >= ampparams[2 * polyn]) {//*(ampparams.end()-2)
			ampvalue = ampparams.back();
		}
	}

}

