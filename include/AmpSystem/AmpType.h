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
//+++          For example: Tablular Amp, Periodic Amp, Smooth Amp, and so on
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class AmpType {//http://127.0.0.1:2080/v6.13/books/usi/default.htm?startat=pt06ch57.html#usi-amp  //http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch01abk09.html#usb-kws-mamplitude
	CONSTAMP,	   //http://127.0.0.1:2080/v6.13/books/usi/default.htm?startat=pt06ch57s03.html
	TABULAR,
	EQUALSPACED,
	PERIODIC,//
	MODULATED,//
	DECAY,
	SOLUTIONBASED,
	SMOOTHSTEP,
	SPECTRUM,
	PSD,
	//...
	USER1A,
	USER2A,
	USER3A,
	USER4A,
	USER5A,
	USER6A,
	USER7A,
	USER8A,
	USER9A,
	USER10A
};














