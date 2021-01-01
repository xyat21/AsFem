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
//+++ Date   : 2020.07.12
//+++ Purpose: This function can read the [qpoint] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadQPointBlock(ifstream &in,string str,int &linenum,FE &fe){
    // format:
    // [qpoint]
    //   type=gauss[gausslobatto]
    //   order=2
    //   bcorder=1
    // [end]

    bool HasType;
    vector<double> numbers;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    str=StringUtils::StrToLower(str);
    numbers.clear();

    HasType=false;
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StringUtils::StrToLower(str);
            continue;
        }
        str=StringUtils::RemoveStrSpace(str);
        
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            if((substr.find("gauss")!=string::npos||
               substr.find("GAUSS")!=string::npos)&&substr.size()==5){
                   fe.SetQPointType(QPointType::GAUSSLEGENDRE);
                   HasType=true;
            }
            else if((substr.find("gausslobatto")!=string::npos||
                    substr.find("GAUSSLOBATTO")!=string::npos)&&substr.size()==12){
                   fe.SetQPointType(QPointType::GAUSSLOBATTO);
                   HasType=true;
            }
            else{
                HasType=false;
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("unsupported gauss point type in the [qpoint] block, type=gauss[gausslobatto] is expected in the [qpoint] block",false);
                MessagePrinter::AsFem_Exit();
            }
        }
        else if((str.find("order=")!=string::npos||
                 str.find("ORDER=")!=string::npos)&&
                 str.find("bcorder")==string::npos&&
                 str.find("BCORDER")==string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("order= can not be found in the [qpoint] block, 'order=integer' is expected in the [qpoint] block",false);
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(int(numbers[0])<0||int(numbers[0])>7){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid order value in the [qpoint] block, 'order=integer' is expected in the [qpoint] block",false);
                    MessagePrinter::AsFem_Exit();
                }
                fe.SetBulkQpOrder(int(numbers[0]));
                fe.SetBCQpOrder(int(numbers[0]));
            }
        }
        else if(str.find("bcorder=")!=string::npos||
                str.find("BCORDER=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("bcorder= can not be found in the [qpoint] block, 'bcorder=integer' is expected in the [qpoint] block",false);
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(int(numbers[0])<0||int(numbers[0])>7){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid bcorder value in the [qpoint] block, 'bcorder=integer' is expected in the [qpoint] block",false);
                    MessagePrinter::AsFem_Exit();
                }
                fe.SetBCQpOrder(int(numbers[0]));
            }
        }
        else if(str.find("[]")!=string::npos){
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("the block bracket pair is not complete in the [qpoint] block, please check your input file",false);
            MessagePrinter::AsFem_Exit();
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("unknown option in the [qpoint] block, please check your input file",false);
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
    }

    if(!HasType){
        fe.SetQPointType(QPointType::GAUSSLEGENDRE);
        HasType=true;
    }
    
    return HasType;
}