//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the Amplitude curve in VFFPM
//+++          here we define the Decay Amp
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void AmpSystem::ComputeDecayAmp(const double &t, const vector<double> &ampparams, double &ampvalue) {//http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch01abk09.html#usb-kws-mamplitude
	double A0 = ampparams[0];//the constant term.
	double A = ampparams[1];//the coefficient of the exponential function.
	double t0 = ampparams[2];//the starting time.
	double td = ampparams[3];//the decay time of the exponential function.
	double ai=0.0;
	ampvalue = 0.0;

	if (t >= t0) { ai = A * exp(-(t - t0)/td); }
	ampvalue = A0 + ai;
}
