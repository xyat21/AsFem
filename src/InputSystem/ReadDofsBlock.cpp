//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.01
//+++ Purpose: This function can read the [dofs] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.04.16  add abaqus inp dof keylines info.

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsBlock(ifstream &in, string str, int &linenum, Mesh &mesh, DofHandler &dofHandler) {//Mesh是加的
	// dof block format:
    // [dofs]
    //   name=c1 c2		[c1 c2是指自由度代号吗？]
    // [end]
    char buff[55];
    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    // str=StrToLower(str);
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            // str=StrToLower(str);
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos){
			if (str.find("name=fromAbaqus") != string::npos ||
				str.find("NAME=fromABAQUS") != string::npos) {

				PetscPrintf(PETSC_COMM_WORLD, "*** Read [dofs] block from abaqus file !!!                            ***\n");
				HasName = true;
				//获取namelist传回，然后继续***********************************
				//===============================补充成读ABAQUS代码====================================================================
				//abaqus的BC信息保存形式为Step.boundaryDoF.push_back(make_pair(bcset, bcdof)); 从其中提取bcset字符串
				vector<std::pair<string, vector<int>> > _BoudaryDofFromInp;
				//_BoudaryDofFromInp=Mesh::GetBoundaryDofFromInp(Mesh::GetAbaqusFileName());
				_BoudaryDofFromInp = mesh.GetBoundaryDoFListPtr();
				//vector<string> bcsets;
				vector<string> bcdofs;
				for (auto it : _BoudaryDofFromInp) {
					for (auto iv : it.second) {
						bcdofs.push_back(to_string(iv));//模型约束涉及到的所有dofs
					}
					//bcdofs.push_back(to_string(it.second));//first is bcset name, second is the constrained dof nums.
				}
				sort(bcdofs.begin(), bcdofs.end());//非降序排序
				bcdofs.erase(unique(bcdofs.begin(), bcdofs.end()), bcdofs.end());//把重复的元素删掉，后面的位置向前移
				//vector<std::pair<vector<string>, double> > Mesh::GetDsLoadFromInp                  vector<std::pair<vector<string>, vector<double>> > Mesh::GetDloadsFromInp
				vector<std::pair<vector<string>, double> > _DsloadFromInp;
				_DsloadFromInp = mesh.GetDsLoadListPtr();
				if (_DsloadFromInp.size() > 0) {
					bcdofs.push_back(_DsloadFromInp[0].first[2]);//P 增加压力载荷类型当一个自由度		后续考虑其他类型载荷、目前只将pressure类放入了dof
				}

				namelist = bcdofs;
				bcdofs.clear();
				if (namelist.size() > 0) {
					if (StringUtils::IsUniqueStrVec(namelist)) {
						HasName = true;
						dofHandler.AddDofNameFromStrVec(namelist);
					}
					else {
						snprintf(buff, 55, "line-%d has some errors", linenum);
						MessagePrinter::PrintErrorTxt(string(buff));
						MessagePrinter::PrintErrorTxt(" duplicated dof name from abaqus in [dofs] block !!!");
						MessagePrinter::AsFem_Exit();
					}
				}
				else {
					HasName = false;
					snprintf(buff, 55, "line-%d has some errors", linenum);
					MessagePrinter::PrintErrorTxt(string(buff));
					MessagePrinter::PrintErrorTxt(" no dofsets from ABAQUS found in [dofs] block, name=ABAQUS so inp information ... is expected");
					MessagePrinter::AsFem_Exit();
				}
				//===↑ 与下面读取namelist后操作相同
			}
			else {
				int i=str.find_first_of('=');
				if(str.size()<=5){
					snprintf(buff,55,"line-%d has some errors",linenum);
					MessagePrinter::PrintErrorTxt(string(buff));
					MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
					MessagePrinter::AsFem_Exit();
				}
				string substr=str.substr(i+1,str.length());
				namelist=StringUtils::SplitStr(substr,' ');
				if(namelist.size()>0){
					if(StringUtils::IsUniqueStrVec(namelist)){
						HasName=true;
						dofHandler.AddDofNameFromStrVec(namelist);
					}
					else{
						snprintf(buff,55,"line-%d has some errors",linenum);
						MessagePrinter::PrintErrorTxt(string(buff));
						MessagePrinter::PrintErrorTxt(" duplicated dof name in the [dofs] block");
						MessagePrinter::AsFem_Exit();
					}
				}
				else{
					HasName=false;

					snprintf(buff,55,"line-%d has some errors",linenum);
					MessagePrinter::PrintErrorTxt(string(buff));
					MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
					MessagePrinter::AsFem_Exit();
				}            
			}
		}
        else if(str.find("[]")!=string::npos){
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt(" the bracket pair is not complete in the [dofs] block");
            MessagePrinter::AsFem_Exit();
        }
        else{
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt(" unknown option in the [dofs] block");
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
    }
    if(!HasName){
        MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
        MessagePrinter::AsFem_Exit();
    }
    return HasName;
}