//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the Amplitude curve in VFFPM
//+++          here we define the Periodic Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void AmpSystem::ComputePeriodicAmp(const double &t, const vector<double> &ampparams, double &ampvalue) {//http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch01abk09.html#usb-kws-mamplitude
	int N = (int) ampparams[0];//the number of terms in the Fourier series.  An Bn
	double w = ampparams[1];//the circular frequency, in radians per time.
	double t0 = ampparams[2];//the starting time.
	double A0 = ampparams[3];//the constant term in the Fourier series.
	double ai;
	ampvalue = 0.0;
	//for ampparams[i>3] refers to An,Bn data.
	if (t<t0) {
		ampvalue = A0;
	}
	else {
		ampvalue = A0;
		for (int n = 1; n <= N;++n) {
			ai = (ampparams[2+2*n])*cos(n*w*(t-t0)) +(ampparams[3 + 2 * n])*sin(n*w*(t - t0));
			ampvalue += ai;
		}
	}
}

