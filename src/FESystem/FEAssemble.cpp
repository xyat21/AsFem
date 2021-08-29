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
//+++ Date   : 2020.11.30
//+++ Purpose: Implement the different assemble functions for different
//+++          cases, i.e., residual assemble and jacobian assemble
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.25  add secA for shell sumR

#include "FESystem/FESystem.h"

void FESystem::AssembleSubResidualToLocalResidual(const int &ndofspernode,const int &dofs,const int &iInd,
                                            const VectorXd &subR,VectorXd &localR){
    // TODO: we need a better assemble algorithm for complex coupling case !!!
    for(int i=1;i<=dofs;i++){
        localR((iInd-1)*ndofspernode+i)+=subR(i);//单元iInd节点的i自由度的值 i重复[因为不同高斯点在同一节点作用/单元属于多sub block?]时累加  
												 //按总自由度编号存储转换为节点对应自由度编号存储 本质上还是1-dofs
    }
}
//***********************************
//void FESystem::AccumulateLocalResidual(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
//                                const VectorXd &localR,vector<double> &sumR){
void FESystem::AccumulateLocalResidual(const int &dofs, const vector<double> &dofsactiveflag, const double &JxW, const double &secA,
	const VectorXd &localR,vector<double> &sumR){
    for(int i=1;i<=dofs;i++){
        //sumR[i-1]+=localR(i)*JxW*dofsactiveflag[i-1];//没计shell 厚度
		sumR[i-1]+= localR(i)*secA*JxW;//只考虑内力，因此不管约束自由度导致合力为0，在其他地方单独计算支反力，然后合力使其为0.
	}
}
//*****************************************
void FESystem::AssembleLocalResidualToGlobalResidual(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &residual,Vec &rhs){
    VecSetValues(rhs,ndofs,dofindex.data(),residual.data(),ADD_VALUES);
}
//*************************************************************
void FESystem::AssembleSubJacobianToLocalJacobian(const int &ndofspernode,
                                            const int &iInd,const int &jInd,
                                            const MatrixXd &subK,MatrixXd &localK){
    for(int i=1;i<=ndofspernode;i++){
        for(int j=1;j<=ndofspernode;j++){
            localK((iInd-1)*ndofspernode+i,(jInd-1)*ndofspernode+j)+=subK(i,j);
        }
    }
}
void FESystem::AccumulateLocalJacobian(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                const MatrixXd &localK,vector<double> &sumK){
    for(int i=1;i<=dofs;i++){
        if(dofsactiveflag[i-1]>0.0){
            for(int j=1;j<=dofs;j++){
                sumK[(i-1)*dofs+j-1]+=localK(i,j)*JxW;
                if(localK(i,j)*JxW>_MaxKMatrixValue) _MaxKMatrixValue=localK(i,j)*JxW;
            }
        }
    }
}
void FESystem::AssembleLocalJacobianToGlobalJacobian(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &jacobian,Mat &K){
    MatSetValues(K,ndofs,dofindex.data(),ndofs,dofindex.data(),jacobian.data(),ADD_VALUES);
}
//**********************************************************************
void FESystem::AssembleSubHistToLocal(const int &e,const int &ngp,const int &gpInd,const Materials &mate,SolutionSystem &solutionSystem){
    solutionSystem._ScalarMaterials[(e-1)*ngp+gpInd-1]=mate.ScalarMaterials;
    solutionSystem._VectorMaterials[(e-1)*ngp+gpInd-1]=mate.VectorMaterials;
    solutionSystem._Rank2TensorMaterials[(e-1)*ngp+gpInd-1]=mate.Rank2Materials;
    solutionSystem._Rank4TensorMaterials[(e-1)*ngp+gpInd-1]=mate.Rank4Materials;
}
void FESystem::AssembleLocalHistToGlobal(const int &e,const int &ngp,SolutionSystem &solutionSystem){
    for(int gpInd=0;gpInd<ngp;gpInd++){
        solutionSystem._ScalarMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._ScalarMaterials[(e-1)*ngp+gpInd];
        solutionSystem._VectorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._VectorMaterials[(e-1)*ngp+gpInd];
        solutionSystem._Rank2TensorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._Rank2TensorMaterials[(e-1)*ngp+gpInd];
        solutionSystem._Rank4TensorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._Rank4TensorMaterials[(e-1)*ngp+gpInd];
    }
}