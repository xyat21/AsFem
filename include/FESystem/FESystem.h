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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.27  add funs to achieve virtual inverse motion

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

#include "petsc.h"


//*************************************
//*** For AsFem's own header file
//*************************************
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCSystem/BCSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "FE/FE.h"
#include "FE/ShapeFun.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

using namespace std;

class FESystem{
public:
    FESystem();
    void InitBulkFESystem(const Mesh &mesh,
                    const DofHandler &dofHandler,
                    FE &fe,
                    SolutionSystem &solution);

public:
    void SetHistNumPerGPoint(const int &n){_nHist=n;}
    void SetProjNumPerNode(const int &n) {_nProj=n;}
    void SetKMatrixFactor(const double &factor){_KMatrixFactor=factor;}
    inline double GetKMatrixFactor() const{return _KMatrixFactor;}
    void ResetMaxAMatrixValue(){_MaxKMatrixValue=-1.0e6;}
    void SetMaxAMatrixValue(const double &val) {_MaxKMatrixValue=val;}
    inline double GetMaxAMatrixValue()const {return _MaxKMatrixValue;}
    inline double GetBulkVolume() const {return _BulkVolumes;}
	void SetMinLcs(const double &Lcs) { _Lcs = Lcs; }
	inline double GetMinLcs()const { return _Lcs; }

    // for FEM simulation related functions
    void FormBulkFE(const FECalcType &calctype,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                SolutionSystem &solutionSystem,
				Mat &AMATRIX, Vec &MASS, Vec &CDAMP, Vec &RHS);
    
	//void Initdt(Mesh &mesh, ElmtSystem &elmtSystem, MateSystem &mateSystem);//为了通过单元面积 尺寸计算显式临界步长
	PetscReal Initdtcr(Mesh &mesh);
	//PetscReal Dtcr(Mesh &mesh);
	void CalcRF(DofHandler &dofHandler, Vec &REXT, Vec &RINT, Vec &RHS, Vec &RF);//update reaction force, xxxand recalculate RHSxxx.
	void CalcreRF(Vec &RF, Vec &reRF);
	void CalcUtt(Vec &M, Vec &RHS, Vec &Utt);
	void InitUold(const double &dt, Vec &CDAMP, Vec &Uold, Vec &U, Vec &V, Vec &Utt);
	void SolveUnew(const double &dt, Vec &CDAMP, Vec &Uold, Vec &U, Vec &Unew, Vec &Utt);

	//void SolveUd(Mesh &mesh, const DofHandler &dofHandler, SolutionSystem &solutionSystem, Vec &Nd);// inverse motion to get deformation
	void CalctranU(const int &nNodes, const int &nDim, const int &nDofs, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Ub, vector<Vector3d> &tranU);
	void CalcOutPlaneRotate(const int &nNodes, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Ub, vector<Vector3d> &Xai1, vector<Vector3d> &Xbi1, vector<Vector3d> &eai1, vector<Vector3d> &ebi1, vector<double> &Xai1norm, vector<double> &Xbi1norm, Vector3d &thetaout, Vector3d &nout, vector<Vector3d> &Utout);
	void CalcInPlaneRotate(const int &nNodes, const int &nDim, vector<Vector3d> &Xai1, vector<Vector3d> &Xbi1, Vector3d &na, Vector3d &thetain);
	//void CalcRotateVec(Vector3d &thetaout, Vector3d &thetain, Vector3d &nrotate, double &totalrotate);
	void SolveUdeform(const int &nNodes, vector<Vector3d> &tranU, vector<Vector3d> &Xbi1, Vector3d &nrotate, double &totalrotate, vector<Vector3d> &Ud);
	void SolveUlocal(const int &nNodes, const int &nDim, const int &nDofs, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Xbi1, vector<Vector3d> &Ud, Vector3d &na, Vector3d &nrotate, vector<Vector3d> &Xlocal, vector<Vector3d> &Ulocal, vector<double> &elUlocal);
	void RenewRcomponts(const int &nDofsPerNode, const vector<Vector3d> &Xlocal, vector<double> &_R);
	void TransformLocalRToGlobalR(const int &nNodes, const int &nDim, const int &nDofs, const MatrixXd &Atheta, const MatrixXd &Ql, const double &totalrotate, vector<double> &_R);
	//void TransformLocalDletaSigmaToGlobal(const int &nNodes, const int &nDim, const int &nDofs, const MatrixXd &Ql, Rank2MateType &Rank2Materials);
	void TransformLocalDletaSigmaToGlobal(const MatrixXd &Ql, Rank2MateType &Rank2Materials);
	void RenewReacforceViaSect(const double &secA, map<string, double> &gpProj);


private:
    //*********************************************************
    //*** assemble residual to local and global one
    //*********************************************************
    void AssembleSubResidualToLocalResidual(const int &ndofspernode,const int &dofs,const int &iInd,
                                            const VectorXd &subR,VectorXd &localR);
	//void AccumulateLocalResidual(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
	//                             const VectorXd &localR,vector<double> &sumR);
	void AccumulateLocalResidual(const int &dofs, const vector<double> &dofsactiveflag, const double &JxW, const double &secA,
								 const VectorXd &localR,vector<double> &sumR);
    void AssembleLocalResidualToGlobalResidual(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &residual,Vec &rhs);

    //*********************************************************
    //*** assemble jacobian to local and global one
    //*********************************************************
    void AssembleSubJacobianToLocalJacobian(const int &ndofspernode,
                                            const int &iInd,const int &jInd,
                                            const MatrixXd &subK,MatrixXd &localK);
    void AccumulateLocalJacobian(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                 const MatrixXd &localK,vector<double> &sumK);
    void AssembleLocalJacobianToGlobalJacobian(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &jacobian,Mat &K);

    void AssembleLocalToGlobal(const int &isw,const int &ndofs,vector<int> &elDofs,
                               vector<double> &localK,vector<double> &localR,
                               Mat &AMATRIX,Vec &RHS);
    //*********************************************************
    //*** for projection
    //*********************************************************
    void AssembleLocalProjectionToGlobal(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                         const map<string,double> &ProjVariables,
                                         const ScalarMateType &ScalarMate,
                                         const VectorMateType &VectorMate,
                                         const Rank2MateType &Rank2Mate,
                                         const Rank4MateType &Rank4Mate,
                                         SolutionSystem &solutionSystem);

    void AssembleLocalProjVariable2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                          const int &nProj,vector<string> ProjNameVec,const map<string,double> &elProj,Vec &ProjVec);

    void AssembleLocalProjScalarMate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                            const int &nProj,vector<string> ProjNameVec,const ScalarMateType &ScalarMate,Vec &ProjVec);

    void AssembleLocalProjVectorMate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                            const int &nProj,vector<string> ProjNameVec,const VectorMateType &VectorMate,Vec &ProjVec);

    void AssembleLocalProjRank2Mate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                            const int &nProj,vector<string> ProjNameVec,const Rank2MateType &Rank2Mate,Vec &ProjVec);

    void AssembleLocalProjRank4Mate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                            const int &nProj,vector<string> ProjNameVec,const Rank4MateType &Rank4Mate,Vec &ProjVec);

    void Projection(const int &nTotalNodes,SolutionSystem &solutionSystem);


    //*********************************************************
    //*** for history variables
    //*********************************************************
    void AssembleSubHistToLocal(const int &e,const int &ngp,const int &gpInd,const Materials &mate,SolutionSystem &solutionSystem);
    void AssembleLocalHistToGlobal(const int &e,const int &ngp,SolutionSystem &solutionSystem);
    

public:
    void PrintFESystemInfo() const;



private:
    double _BulkVolumes=0.0;
    MatrixXd _localK;    //used in uel, the size is the maximum dofs per element
	VectorXd _localR, _localM, _localC;    //used in uel, the size is the maximum dofs per element
	MatrixXd _subK; // used in each sub element, the size is the maximum dofs per node
	VectorXd _subR, _subM, _subC; // used in each sub element, the size is the maximum dofs per node
	vector<double> _K, _R, _M, _C;//used in assemble

	vector<double> _Lc;//vec for charateristic length of all elements
	double _Lcs;// minimum characteristic length

    Nodes _elNodes;
    vector<int> _elConn,_elDofs;
    vector<double> _elDofsActiveFlag;
	vector<double> _elU, _elUo, _elV, _elRF, _elUtt;
    vector<double> _elUold,_elVold;
	vector<double> _gpU, _gpV, _gpUtt;

    vector<double> _gpUOld,_gpVOld;
    vector<double> _gpHist,_gpHistOld;
    map<string,double> _gpProj;
	vector<Vector3d> _gpGradU, _gpGradV, _gpGradUtt;

    vector<Vector3d> _gpGradUOld,_gpGradVOld;
    vector<double> _MaterialValues;
    Vector3d _gpCoord;
    int _nHist,_nProj,_nGPoints;
    double _MaxKMatrixValue=-1.0e9,_KMatrixFactor=0.1;

    ElmtType elmttype;
    MateType matetype;
    vector<int> localDofIndex;
    int mateindex;

	// XY for inverse deformation calculation
	MeshType emeshtype;
	vector<Vector3d> Xe, Ub, Ua, Ud, Ulocal, Xlocal;//nodes coords,new disp., last disp.  pure deformation U after inverse motion.  local U.  local X accordinates.
	vector<Vector3d> tranU, Utout, Utoi;//relative vector for other nodes to 1st node ->translate U ->+out plane rotate U ->+in plane rotate U 
	vector<Vector3d> Xai1, Xbi1, eai1, ebi1;
	vector<double> Xai1norm, Xbi1norm, elUlocal;
	Vector3d na, nb;//a b 时刻的12 13边长向量外积，法向量
	Vector3d thetaout, thetain;
	MatrixXd Atheta, Ql;
	double totalrotate;
	Vector3d nout, nin, nrotate;

private:
    //************************************
    //*** For PETSc related vairables
    PetscMPIInt _rank,_size;
	VecScatter _scatteru, _scatteruo, _scatterv, _scatterutt, _scatterrf, _scatterproj, _scatterhist, _scatterhistold, _scatteruold, _scattervold;
    Vec _Useq, _Uoseq,_Uoldseq;// this can contain the ghost node from other processor
    Vec _Vseq,_Voldseq, _U_ttseq, _RFseq;
    Vec _ProjSeq, _HistSeq, _HistOldSeq;
};