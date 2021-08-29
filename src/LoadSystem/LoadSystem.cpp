//****************************************************************

//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.14
//+++ Purpose: Define the load condition system in VFIFPM
//+++          Here we can apply:
//+++               1) cload, i.e. concentrated force, moment ...
//+++               2) dload, i.e. gravity, body force
//+++               3) dsload as well as user-defined-load (ulc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "LoadSystem/LoadSystem.h"

LoadSystem::LoadSystem(){
    _nLoadBlocks=0;
    _LoadBlockList.clear();

    _xi=0.0;_eta=0.0;_JxW=0.0;
    
    _xs[0][0]=0.0;_xs[0][1]=0.0;_xs[0][2]=0.0;
    _xs[1][0]=0.0;_xs[1][1]=0.0;_xs[1][2]=0.0;
    _xs[2][0]=0.0;_xs[2][1]=0.0;_xs[2][2]=0.0;

    _normals=0.0;
}

//************************************
void LoadSystem::AddLoadBlock2List(LoadBlock &loadblock){
    string msg;
    if(_LoadBlockList.size()<1){
        _LoadBlockList.push_back(loadblock);
        _nLoadBlocks=int(_LoadBlockList.size());
    }
    else{
        bool NotInList=true;
        for(const auto &it:_LoadBlockList){
            if(it._LoadBlockName==loadblock._LoadBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            _LoadBlockList.push_back(loadblock);
            _nLoadBlocks=static_cast<int>(_LoadBlockList.size());
        }
        else{
            msg="duplicated ["+loadblock._LoadBlockName+"] in the [loads] sub block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}

void LoadSystem::InitLoadSystem(const Mesh &mesh){
    //_elNodes.InitNodes(mesh.GetBulkMeshNodesNumPerBulkElmt());

    _xi=0.0;_eta=0.0;_JxW=0.0;
    
    _xs[0][0]=0.0;_xs[0][1]=0.0;_xs[0][2]=0.0;
    _xs[1][0]=0.0;_xs[1][1]=0.0;_xs[1][2]=0.0;
    _xs[2][0]=0.0;_xs[2][1]=0.0;_xs[2][2]=0.0;

    _normals=0.0;
}