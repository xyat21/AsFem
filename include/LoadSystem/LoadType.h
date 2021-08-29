//****************************************************************
//* This file is part of the VFFPM framework
//* Vector Form Finite Particle Method program (VFFPM)
//* All rights reserved, Xiaoyuan @ CopyRight 2021
//* https://github.com/XiaoyuanGao/LDFSI.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.09
//+++ Purpose: Define the loading condition type in VFFPM
//+++          For example: Concentrated Loads, Pressure Loads, Body Loads, and so on
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once
//Cload		1-6					[concentrated load, moment]						[unit: F, FL]
//Dsload	TRSHR/TRVEC/P		[pressure, shear traction,general traction ]	[unit: F/L2, ...] /former is vec before project to local direc, latter is global direc./
//Dload		GRAV/BX BY BZ		[gravity, body force]							[unit: F/L3]
enum class LoadType {//http://127.0.0.1:2080/v6.13/books/usi/default.htm?startat=pt03ch16s09.html#usi-lbi-loadeditors
	NULLLC,
	CONCENTRATEDLC,
	MOMENTLC,
	PRESSURELC,//http://127.0.0.1:2080/v6.13/books/usb/default.htm?startat=pt07ch34s04aus122.html “‘œ¬ «DLOAD	//Data lines to define distributed surface pressures:  1.Surface name. 2.Distributed load type label P, PNU, SP, or VP.  3.Reference load magnitude, which can be modified by using the *AMPLITUDE option.
	SHELLEDGELC,//http://127.0.0.1:2080/v6.13/books/key/default.htm?startat=ch04abk31.html#usb-kws-hdload
	SURFACETRACTIONLC,
	//NODALFORCELC,
	//NODALFLUXLC,
	PIPEPRESSURELC,
	BODYLC,
	LINELC,
	GRAVITYLC,
	//...
	USER1LC,
	USER2LC,
	USER3LC,
	USER4LC,
	USER5LC,
	USER6LC,
	USER7LC,
	USER8LC,
	USER9LC,
	USER10LC
};

enum class LoadOP {
	MOD,
	NEW
};














