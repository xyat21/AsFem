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
//+++ Date   : 2020.12.30
//+++ Purpose: here we call the TS+SNES solver from PETSc to solve
//+++          our nonlinear system equation R(x,xdot,t)->0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"


//***************************************************************
//*** here we define a monitor to print out the iteration info
//***************************************************************
PetscErrorCode MySNESMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
    char buff[70];
    string str;
    TSAppCtx *user=(TSAppCtx*)ctx;
    user->iters=iters;
    user->rnorm=rnorm;
    if(iters==0){
        SNESGetSolutionNorm(snes,&user->dunorm);
    }
    else{
        SNESGetUpdateNorm(snes,&user->dunorm);
    }
    user->enorm=rnorm*user->dunorm;
    if(iters==0){
        user->rnorm0=rnorm;
        user->dunorm0=user->dunorm;
        user->enorm0=user->enorm;
    }
    if(user->IsDepDebug){
        snprintf(buff,70,"  SNES solver:iters=%3d,|R|=%11.4e,|dU|=%11.4e,dt=%7.2e",iters,rnorm,user->dunorm,user->dt);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }
    return 0;
}
//**********************************
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;
    char buff[68];
    string str;
    double dt;

    TSGetTimeStep(ts,&dt);

    user->time=time;
    user->step=step;
    user->dt=dt;
    // update previous solution
    VecCopy(user->_solutionSystem._Unew,user->_solutionSystem._Uold);
    VecCopy(user->_solutionSystem._V,user->_solutionSystem._Vold);
    // update current solution
    VecCopy(U,user->_solutionSystem._Unew);
    
    snprintf(buff,68,"Time step=%8d, time=%13.5e, dt=%13.5e",step,time,dt);
	//if (step < 0) return 0; /* step of -1 indicates an interpolated solution  XY*/
	str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(!user->IsDepDebug){
        snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",user->iters+1,user->rnorm0,user->rnorm);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }


    if(user->_fectrlinfo.IsProjection){
        user->_feSystem.FormBulkFE(FECalcType::Projection,time,dt,user->_fectrlinfo.ctan,
                                   user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,user->_mateSystem,
                                   //user->_solutionSystem,user->_equationSystem._AMATRIX,user->_equationSystem._RHS);
								   user->_solutionSystem, user->_equationSystem._AMATRIX, user->_equationSystem._Mass, user->_equationSystem._Cdamp, user->_equationSystem._RHS);

    }
    if(step%user->_outputSystem.GetIntervalNum()==0){
        user->_outputSystem.WriteResultToFile(step,user->_mesh,user->_dofHandler,user->_solutionSystem);
        user->_outputSystem.WriteResultToPVDFile(time,user->_outputSystem.GetOutputFileName());
        MessagePrinter::PrintNormalTxt("Write result to "+user->_outputSystem.GetOutputFileName());
        MessagePrinter::PrintDashLine();

    }
    if(step%user->_postprocess.GetOutputIntervalNum()==0){
        user->_postprocess.RunPostprocess(time,user->_mesh,user->_dofHandler,user->_fe,user->_solutionSystem);
    }


    // update history variable
    user->_feSystem.FormBulkFE(FECalcType::UpdateHistoryVariable,time,dt,user->_fectrlinfo.ctan,
                               user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,user->_mateSystem,
                               //user->_solutionSystem,user->_equationSystem._AMATRIX,user->_equationSystem._RHS);
							   user->_solutionSystem, user->_equationSystem._AMATRIX, user->_equationSystem._Mass, user->_equationSystem._Cdamp, user->_equationSystem._RHS);


    if(user->IsAdaptive&&step>=1){
        if(user->iters<=user->OptiIters){
            dt=user->dt*user->GrowthFactor;
            if(dt>user->DtMax) dt=user->DtMax;
        }
        else{
            dt=user->dt*user->CutbackFactor;
            if(dt<user->DtMin) dt=user->DtMin;
        }
        TSSetTimeStep(ts,dt);
    }

    return 0;
}
//***************************************************************
//*** for our Residual and Jacobian calculation
//***************************************************************
PetscErrorCode ComputeIResidual(TS ts,PetscReal t,Vec U,Vec V,Vec RHS,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;

    TSGetTimeStep(ts,&user->dt);
    // apply the initial dirichlet boundary condition
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);

    VecCopy(U,user->_solutionSystem._Unew);
    VecCopy(V,user->_solutionSystem._V);

    user->_feSystem.FormBulkFE(FECalcType::ComputeResidual,
                        t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        user->_solutionSystem,
                        //user->_equationSystem._AMATRIX,RHS);
						user->_equationSystem._AMATRIX, user->_equationSystem._Mass, user->_equationSystem._Cdamp, RHS);

    
    user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                    FECalcType::ComputeResidual,t,user->_fectrlinfo.ctan,U,
                    user->_equationSystem._AMATRIX,RHS);
    
    return 0;
}
//******************************************************
PetscErrorCode ComputeIJacobian(TS ts,PetscReal t,Vec U,Vec V,PetscReal s,Mat A,Mat B,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;
    int i;

    TSGetTimeStep(ts,&user->_fectrlinfo.dt);
    TSGetTimeStep(ts,&user->dt);
    user->_feSystem.ResetMaxAMatrixValue();// we reset the penalty factor
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);

    VecCopy(U,user->_solutionSystem._Unew);
    VecCopy(V,user->_solutionSystem._V);

    user->_fectrlinfo.ctan[0]=1.0;
    user->_fectrlinfo.ctan[1]=s;// dUdot/dU		//matrix dF/dU + s*dF/dU_t两项d的系数？  to be the Jacobian of F(t,U,W+a*U) where F(t,U,U_t) = 0 is the DAE to be solved. 
								//https://www.mcs.anl.gov/petsc/petsc-dev/docs/manualpages/TS/TSSetIJacobian.html

    user->_feSystem.FormBulkFE(FECalcType::ComputeJacobian,
                        t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        user->_solutionSystem,
                        //A,user->_equationSystem._RHS);
						A, user->_equationSystem._Mass, user->_equationSystem._Cdamp, user->_equationSystem._RHS);

    
    if(user->_feSystem.GetMaxAMatrixValue()>1.0e12){
        user->_bcSystem.SetBCPenaltyFactor(1.0e20);
    }
    else if(user->_feSystem.GetMaxAMatrixValue()>1.0e6&&user->_feSystem.GetMaxAMatrixValue()<=1.0e12){
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);
    }
    else if(user->_feSystem.GetMaxAMatrixValue()>1.0e3&&user->_feSystem.GetMaxAMatrixValue()<=1.0e6){
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e12);
    }
    else{
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e16);
    }

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                    FECalcType::ComputeJacobian,t,user->_fectrlinfo.ctan,U,
                    A,user->_equationSystem._RHS);

    MatGetSize(B,&i,&i);

    return 0;
}

//*************************************************************************
bool TimeStepping::Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
			BCSystem &bcSystem, ICSystem &icSystem, AmpSystem &ampSystem, LoadSystem &loadSystem,
			SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            OutputSystem &outputSystem,
            Postprocess &postprocessSystem,
            FEControlInfo &fectrlinfo){

    _appctx=TSAppCtx{mesh,dofHandler,
				   bcSystem,icSystem,ampSystem, loadSystem,
                   elmtSystem,mateSystem,
                   solutionSystem,equationSystem,
                   fe,feSystem,
                   outputSystem,
                   postprocessSystem,
                   fectrlinfo,
                   //*****************
                   0.0,0.0,
                   0.0,0.0,
                   0.0,0.0,
                   0,
                   fectrlinfo.IsDebug,fectrlinfo.IsDepDebug,
                   0.0,0.0,0,
                    //********************************
                    _Adaptive,
                    _OptIters,
                    _GrowthFactor,_CutBackFactor,
                    _DtMin,_DtMax
                   };
    

	_appctx._icSystem.ApplyIC(_appctx._mesh, _appctx._dofHandler, _appctx._solutionSystem._U);//后面有需要则更改可引入Vnew和Uttnew初始条件.用t时刻，不需要t+1时刻.
	_appctx._bcSystem.ApplyInitialBC(_appctx._mesh, _appctx._dofHandler, 0.0, _appctx._solutionSystem._U);//这两步只更新了BC 与U相关的初始值,insert模式[故对同一节点与IC不能同时存在]
    //_appctx._icSystem.ApplyIC(_appctx._mesh,_appctx._dofHandler,_appctx._solutionSystem._Unew);
    //_appctx._bcSystem.ApplyInitialBC(_appctx._mesh,_appctx._dofHandler,0.0,_appctx._solutionSystem._Unew);
    _appctx._feSystem.FormBulkFE(FECalcType::InitHistoryVariable,_appctx._fectrlinfo.t,_appctx._fectrlinfo.dt,_appctx._fectrlinfo.ctan,
                                 _appctx._mesh,_appctx._dofHandler,_appctx._fe,_appctx._elmtSystem,_appctx._mateSystem,
                                 _appctx._solutionSystem,
                                 //_appctx._equationSystem._AMATRIX,_appctx._equationSystem._RHS);
								_appctx._equationSystem._AMATRIX, _appctx._equationSystem._Mass, _appctx._equationSystem._Cdamp, _appctx._equationSystem._RINT);//计算单元体积 装配各变量(这里是hist相关)矩阵
								//初始化 M C总体矩阵 其他选项再计算单元内力RHS。  粒子法不考虑单元K矩阵.

	//init Mass vec and Damp vec with known values.

	////initial critical time step with 0.9 reduced factor
	PetscReal dt_cr0;
	dt_cr0 = _appctx._feSystem.Initdtcr(_appctx._mesh);//简化下面过程
	_appctx.dt = dt_cr0;//


	//apply external forces: 节点集中力Pi 单元上集中力分配sum(Pij)  单元上均布力分配sum(Pik) 单元体力m*gravity 随体载荷
	_appctx._loadSystem.ApplyLoad(_appctx._mesh, _appctx._dofHandler, _appctx._fe, FECalcType::ComputeResidual, _appctx._fectrlinfo.t, _appctx._solutionSystem._U, _appctx._equationSystem._REXT);
	//RHS is total force, here last account the constraint force, to make constrainted dof to be 0. and return Cforce
	_appctx._feSystem.CalcRF(_appctx._dofHandler, _appctx._equationSystem._REXT, _appctx._equationSystem._RINT, _appctx._equationSystem._RHS, _appctx._equationSystem._RF);
	_appctx._feSystem.CalcreRF(_appctx._equationSystem._RF, _appctx._solutionSystem._reRF);
	_appctx._feSystem.CalcUtt(_appctx._equationSystem._Mass, _appctx._equationSystem._RHS, _appctx._solutionSystem._Utt); 要补角加速度、惯性矩 = > 质量矩阵4 - 6;
	//init solution t-1时刻
	_appctx._feSystem.InitUold(_appctx.dt, _appctx._equationSystem._Cdamp, _appctx._solutionSystem._Uold, _appctx._solutionSystem._U, _appctx._solutionSystem._V, _appctx._solutionSystem._Utt);

	for (_appctx.time = 0.0; _appctx.time < _appctx._fectrlinfo.t; _appctx.time += _appctx.dt) {
		_appctx._feSystem.SolveUnew(_appctx.dt, _appctx._equationSystem._Cdamp, _appctx._solutionSystem._Uold, _appctx._solutionSystem._U, _appctx._solutionSystem._Unew, _appctx._solutionSystem._Utt);
		//逆运动 消除平动和转动的刚体位移，得到变形位移Ud，计算单元内力

		//_appctx._feSystem.SolveUd();		

		_appctx._feSystem.FormBulkFE(FECalcType::ComputeResidual,
			_appctx._fectrlinfo.t, _appctx._fectrlinfo.dt, _appctx._fectrlinfo.ctan,
			_appctx._mesh, _appctx._dofHandler, _appctx._fe,
			_appctx._elmtSystem, _appctx._mateSystem,
			_appctx._solutionSystem,
			_appctx._equationSystem._AMATRIX, _appctx._equationSystem._Mass, _appctx._equationSystem._Cdamp, _appctx._equationSystem._RINT);

		VecCopy(_appctx._solutionSystem._U, _appctx._solutionSystem._Uold);
		VecCopy(_appctx._solutionSystem._Unew, _appctx._solutionSystem._U);


	}






    TSSetIFunction(_ts,_appctx._equationSystem._RHS,ComputeIResidual,&_appctx);
    TSSetIJacobian(_ts,_appctx._equationSystem._AMATRIX,_appctx._equationSystem._AMATRIX,ComputeIJacobian,&_appctx);
    
    TSMonitorSet(_ts,MyTSMonitor,&_appctx,NULL);
    SNESMonitorSet(_snes,MySNESMonitor,&_appctx,0);

    SNESSetForceIteration(_snes,PETSC_TRUE);

    SNESSetFromOptions(_snes);

    TSSetFromOptions(_ts);

    //***************************************************
    //*** before solve,we need to write some basic info
    //*** to pvd file
    //***************************************************
    _appctx._outputSystem.WritePVDFileHeader();
    TSSolve(_ts,_appctx._solutionSystem._Unew);
    _appctx._outputSystem.WritePVDFileEnd();

    return true;
}
