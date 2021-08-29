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
//+++ Date   : 2020.12.26
//+++ Purpose: Define equation system in AsFem, here you can access
//+++          K matrix and Residual of our system equations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.02  add equation supplymenary variables

#pragma once

#include <iostream>

#include "petsc.h"

#include "DofHandler/DofHandler.h"

using namespace std;

class EquationSystem{
public:
    EquationSystem();

    void InitEquationSystem(const int &ndofs,const int &maxrownnz);
    void CreateSparsityPattern(DofHandler &dofHandler);

    void ReleaseMem();

public:
    Mat _AMATRIX;
	Vec _RHS;//right hand total force.
	Vec _RINT;//internal force
	Vec _RF;//constrain induced reaction force
	Vec _REXT;//external force
	Vec _Mass;
	Vec _Cdamp;
private:
    int _nDofs;
};