//****************************************************************

//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.14
//+++ Purpose: here we can apply the different types of load
//+++          condition, and, once more, we only need to focus on
//+++          the calculation on each ... !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "LoadSystem/LoadSystem.h"

void LoadSystem::RunLoadLibs(const LoadType loadtype,const FECalcType &calctype,
	const Vector3d &normals){
    switch (loadtype)
    {
    case LoadType::USER1LC:
        User1LC(calctype, normals);
        break;


    default:
        MessagePrinter::PrintErrorTxt("unsupported load condition type");
        MessagePrinter::AsFem_Exit();
        break;
    }
}