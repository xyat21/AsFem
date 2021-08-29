//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.10
//+++ Purpose: Define the amplitude curve system in VFIFPM
//+++          Here we can apply:
//+++               1)  ...
//+++               2) 
//+++               3)  as well as user defined amp (uamp)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "AmpSystem/AmpSystem.h"

AmpSystem::AmpSystem() {
	_nAmpBlocks = 0;
	_AmpBlockList.clear();

	_ScaleFactor = 1.0;

}

//************************************
void AmpSystem::AddAmpBlock2List(AmpBlock &ampblock) {
    string msg;
    if(_AmpBlockList.size()<1){
        _AmpBlockList.push_back(ampblock);
        _nAmpBlocks=int(_AmpBlockList.size());
    }
    else{
        bool NotInList=true;
        for(const auto &it:_AmpBlockList){
            if(it._AmpBlockName==ampblock._AmpBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            _AmpBlockList.push_back(ampblock);
            _nAmpBlocks=static_cast<int>(_AmpBlockList.size());
        }
        else{
            msg="duplicated ["+ampblock._AmpBlockName+"] in the [amps] sub block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}

void AmpSystem::ComputeInitialAmp(const double &t, double &ampvalue) {
    _ScaleFactor=1.0;
	//create as ramp 0-1 curve below

}