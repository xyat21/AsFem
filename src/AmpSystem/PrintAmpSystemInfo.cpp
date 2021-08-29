//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.10
//+++ Purpose: print out some basic info for AmpSystem class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

void AmpSystem::PrintAmpSystemInfo()const{
    MessagePrinter::PrintNormalTxt("Amplitude curves information summary:");
    const int len=68;
    char buff[len];
    string str;
    for(auto it:_AmpBlockList){
        snprintf(buff,len," +Amplitude block information:");
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   amplitude block name = [%40s]",it._AmpBlockName.c_str());
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        if(it._Smoothdata >0.0){
            snprintf(buff,len,"   amplitude type name  = %15s, smooth data = %14.6e ", it._AmpTypeName.c_str(),it._Smoothdata);
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        else{
            snprintf(buff,len,"   amplitude type name  = %15s",it._AmpTypeName.c_str());
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }

        //*
		snprintf(buff, len, "   Amplitude definition data:");
		str = buff;
		MessagePrinter::PrintNormalTxt(str);
		//*
		str = "   ";
		for(auto ampdata:it._AmpParams){
		    str+= to_string(ampdata)+",  ";
		}
		int lastN = 4;
		str.substr(str.size() - lastN);//去掉最后的",  "一个尾缀	//s = s.substr(0, s.length()-3);		//str.erase (str.end()-n, str.end());
		MessagePrinter::PrintNormalTxt(str);


		//snprintf(buff,len,"   dof name            = %15s, dof index=%2d",it._DofName.c_str(),it._DofID);
  //      str=buff;
  //      MessagePrinter::PrintNormalTxt(str);
  //      //*
  //      snprintf(buff,len,"   boundary value      = %14.6e",it._BCValue);
  //      str=buff;
  //      MessagePrinter::PrintNormalTxt(str);

        ////*
        //str="   amplitude name       =";
        //for(auto ampname:it._AmpNameList){
        //    str+=ampname+" ";
        //}
        //MessagePrinter::PrintNormalTxt(str);
    }
    MessagePrinter::PrintDashLine();
}