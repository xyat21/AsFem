//****************************************************************

//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.14
//+++ Purpose: print out some basic info for LoadSystem class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "LoadSystem/LoadSystem.h"

void LoadSystem::PrintLoadSystemInfo()const{
    MessagePrinter::PrintNormalTxt("Boundary conditions information summary:");
    const int len=68;
    char buff[len];
    string str;
    for(auto it:_LoadBlockList){
        snprintf(buff,len," +Load block information:");
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   load index numbr %3d, load block name = [%40s]", it._LoadIndex,it._LoadBlockName.c_str());
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        if(it._IsFollower){
            snprintf(buff,len,"   load type name  = %15s (Load Follower)",it._LoadTypeName.c_str());
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        else{
            snprintf(buff,len,"   load type name  = %15s",it._LoadTypeName.c_str());
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        //*
		if (it._LoadAmpName.size()>0) {
			snprintf(buff, len, "   load Amplitude name  = %15s", it._LoadAmpName.c_str());
			str = buff;
			MessagePrinter::PrintNormalTxt(str);
		}
		if (it._LoadSetName.size() > 0) {
			snprintf(buff, len, "   load set name and first line = %15s, %10s, %14.6e", it._LoadSetName.c_str(), it._DofsName[0].c_str(), it._LoadValue);
			str = buff;
			MessagePrinter::PrintNormalTxt(str);
		}
		else if (it._Partunit.values.size() > 0) {
			if (it._Partunit.nodesNo.size() > 0) {
				snprintf(buff, len, "   load first line node and value = %5d, %5s, %14.6e", it._Partunit.nodesNo[0], it._Partunit.sDof[0].c_str(), it._Partunit.values[0]);
				str = buff;
				MessagePrinter::PrintNormalTxt(str);
			}
			else if (it._Partunit.elesNo.size() > 0) {
				snprintf(buff, len, "   load first line element and value = %5d, %5s, %14.6e", it._Partunit.elesNo[0], it._Partunit.sDof[0].c_str(), it._Partunit.values[0]);
				str = buff;
				MessagePrinter::PrintNormalTxt(str);
			}
		}
		if (it._direction.norm > 0) {
			snprintf(buff, len, "   load abbr. and value vec  = %10s, %14.6e, %10.8e, %10.8e, %10.8e", it._DofsName[0].c_str(), it._LoadValue,it._direction[0], it._direction[1], it._direction[2]);
			str = buff;
			MessagePrinter::PrintNormalTxt(str);
		}

    }
    MessagePrinter::PrintDashLine();
}