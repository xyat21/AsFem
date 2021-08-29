//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the Amplitude curve in VFFPM
//+++          here we define the Smooth Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void SmoothInterpolation(vector<double> &amppnts, const double &t, double &ampvalue) {
	//amppnts=(x1,y1,x2,y2)  range:x1<=t<=x2
	double kexi;
	kexi = (t- amppnts[0]) / (amppnts[2] - amppnts[0]);
	ampvalue = amppnts[1] + (amppnts[3] - amppnts[1]) *pow(kexi,3.0)*(10-15*kexi+6*kexi*kexi);
}

void AmpSystem::ComputeSmoothAmp(const double &t, const vector<double> &ampparams, double &ampvalue) {//http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch01abk09.html#usb-kws-mamplitude
	//ampparams[0+2n]=time;	ampparams[1+2n]=Amp
	int polyn;
	vector<double> twopoints(4);
	polyn = ampparams.size() / 2 - 1;
	ampvalue = 0.0;
	for (int n = 1; n <= polyn; ++n) {
		if ((t >= ampparams[2 * (n - 1)]) && (t < ampparams[2 * n])) {
			twopoints.clear();
			std::copy(ampparams.begin() + 2 * (n - 1), ampparams.begin() + 2 * (n + 1), twopoints.begin()); //Çø¼ä¿½±´
			SmoothInterpolation(twopoints, t, ampvalue);//ampvalue can be calculated.
			break;
		}
		else if (t < ampparams[0]) {//
			ampvalue = ampparams[1];//.front()
		}
		else if (t >= ampparams[2 * polyn]) {//*(ampparams.end()-2)
			ampvalue = ampparams.back();
		}
	}
}
