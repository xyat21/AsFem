//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the Amplitude curve in VFFPM
//+++          here we define the Modulated Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void AmpSystem::ComputeModulatedAmp(const double &t, const vector<double> &ampparams, double &ampvalue) {//http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch01abk09.html#usb-kws-mamplitude
	double A0 = ampparams[0];
	double A = ampparams[1];
	double t0 = ampparams[2];//the starting time.
	double w1 = ampparams[3];
	double w2 = ampparams[4];
	double ai=0.0;
	ampvalue = 0.0;

	//if (t<=t0) {
	//	ampvalue = A0;
	//}
	//else {
	//	ampvalue = A0;
	//	ai = A*sin(w1*(t-t0)) *sin(w2*(t - t0));
	//	ampvalue += ai;
	//}

	if (t > t0) { ai = A * sin(w1*(t - t0)) *sin(w2*(t - t0)); }
	ampvalue = A0 + ai;
}
