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
//+++ Purpose: Loop all the bulk elements, where we can calculate
//+++          residual, k matrix as well as projection quantities
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.28  add some new funs to implement the virtual inverse motion.

#include "FESystem/FESystem.h"

void FESystem::FormBulkFE(const FECalcType &calctype,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                SolutionSystem &solutionSystem,
				Mat &AMATRIX, Vec &MASS, Vec &CDAMP, Vec &RHS){
    
	PetscReal damp_factor, inodeMass, elsect, rho;
	vector<double> nodeeqVolume, nodeeqMass;
	vector<vector<double>> nodalDofFlag;
	PetscInt nActiveDofs;
	char buff[70];
	string str;
	string str2;
	char buff2[70];

	if(calctype==FECalcType::ComputeResidual){
        VecSet(RHS,0.0);

		nodeeqVolume = mesh.GetBulkMeshNodeeqVolumeListPtr();
		if (mesh.GetBeamSectionNum() > 0) {//这里都默认只有一个截面属性[0]，梁认为是矩形rect,面积是a*b
			elsect = mesh.GetBeamSectionESetsandMatAdListPtr()[0].second.first[0] * mesh.GetBeamSectionESetsandMatAdListPtr()[0].second.first[1];
		}
		else if (mesh.GetSolidSectionNum() > 0) {//截面积
			//elsect = mesh.GetSolidSectionESetsandMatAListPtr()[0].second;//体单元没有截面积啊？
			elsect = 1.0;
		}
		else if (mesh.GetShellSectionNum() > 0) {//壳厚度
			elsect = mesh.GetShellSectionESetsandMatTListPtr()[0].second;
		}
	}
    else if(calctype==FECalcType::ComputeJacobian){
        MatZeroEntries(AMATRIX);
    }

    else if(calctype==FECalcType::Projection){
        VecSet(solutionSystem._Proj,0.0);
        VecSet(solutionSystem._ProjScalarMate,0.0);
        VecSet(solutionSystem._ProjVectorMate,0.0);
        VecSet(solutionSystem._ProjRank2Mate,0.0);
        VecSet(solutionSystem._ProjRank4Mate,0.0);
    }
    else if(calctype==FECalcType::UpdateHistoryVariable||
            calctype==FECalcType::InitHistoryVariable||
            calctype==FECalcType::InitMaterialAndProjection){
        // do nothing!
		// below XY added
		if (calctype == FECalcType::InitHistoryVariable) {//第一次初始化变量时 把质量矩阵M/C计算好 nnodes*ndofs
			damp_factor = mesh.GetMateListPtr()[0].damp_alpha;
			VecSet(MASS, 0.0);
			//VecSet(CDAMP, 0.0);
			VecSet(CDAMP, damp_factor);

			nodeeqVolume = mesh.GetBulkMeshNodeeqVolumeListPtr();
			if (mesh.GetBeamSectionNum() > 0) {//这里都默认只有一个截面属性[0]，梁认为是矩形rect,面积是a*b
				elsect = mesh.GetBeamSectionESetsandMatAdListPtr()[0].second.first[0] * mesh.GetBeamSectionESetsandMatAdListPtr()[0].second.first[1];
			}
			else if (mesh.GetSolidSectionNum() > 0) {//截面积
				//elsect = mesh.GetSolidSectionESetsandMatAListPtr()[0].second;//体单元没有截面积啊？
				elsect = 1.0;

			}
			else if (mesh.GetShellSectionNum() > 0) {//壳厚度
				elsect = mesh.GetShellSectionESetsandMatTListPtr()[0].second;
			}
			rho = mesh.GetMateListPtr()[0].density;//只考虑一种材料
			//nodeeqMass= rho * elsect * nodeeqVolume;
		}

    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in FormBulkFE, please check your code");
        MessagePrinter::AsFem_Exit();
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);

    // we can get the correct value on the ghosted node!
    // please keep in mind, we will always use Unew and V in SNES and TS !!!
    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

	VecScatterCreateToAll(solutionSystem._U, &_scatteruo, &_Uoseq);// XY
	VecScatterBegin(_scatteru, solutionSystem._U, _Uoseq, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(_scatteru, solutionSystem._U, _Uoseq, INSERT_VALUES, SCATTER_FORWARD);

	VecScatterCreateToAll(solutionSystem._Utt, &_scatterutt, &_U_ttseq);
	VecScatterBegin(_scatterutt, solutionSystem._Utt, _U_ttseq, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(_scatterutt, solutionSystem._Utt, _U_ttseq, INSERT_VALUES, SCATTER_FORWARD);

	VecScatterCreateToAll(solutionSystem._HistOld, &_scatterhistold, &_HistOldSeq);
	VecScatterBegin(_scatterhistold, solutionSystem._HistOld, _HistOldSeq, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(_scatterhistold, solutionSystem._HistOld, _HistOldSeq, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterCreateToAll(solutionSystem._V,&_scatterv,&_Vseq);
    VecScatterBegin(_scatterv,solutionSystem._V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterv,solutionSystem._V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);


    // for the disp and velocity in the previous step
    VecScatterCreateToAll(solutionSystem._Uold,&_scatteruold,&_Uoldseq);
    VecScatterBegin(_scatteruold,solutionSystem._Uold,_Uoldseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteruold,solutionSystem._Uold,_Uoldseq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(solutionSystem._Vold,&_scattervold,&_Voldseq);
    VecScatterBegin(_scattervold,solutionSystem._Vold,_Voldseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scattervold,solutionSystem._Vold,_Voldseq,INSERT_VALUES,SCATTER_FORWARD);

    int rankne=mesh.GetBulkMeshBulkElmtsNum()/_size;
    int eStart=_rank*rankne;
    int eEnd=(_rank+1)*rankne;
    if(_rank==_size-1) eEnd=mesh.GetBulkMeshBulkElmtsNum();

    PetscInt nDofs,nNodes,nDofsPerNode,nDofsPerSubElmt,e;
    PetscInt i,j,jj;
    PetscInt nDim,gpInd;
    PetscReal xi,eta,zeta,w,JxW,DetJac,elVolume;
    nDim=mesh.GetDim();

    _BulkVolumes=0.0;
    for(int ee=eStart;ee<eEnd;++ee){
        e=ee+1;
        mesh.GetBulkMeshIthBulkElmtNodes(e,_elNodes);
        mesh.GetBulkMeshIthBulkElmtConn(e,_elConn);
		emeshtype = mesh.GetBulkMeshBulkElmtType();// XY 
        dofHandler.GetIthBulkElmtDofIndex0(e,_elDofs,_elDofsActiveFlag);//vector<vector<double>> _NodalDofFlag
        nDofs=dofHandler.GetIthBulkElmtDofsNum(e);
        nNodes=mesh.GetBulkMeshIthBulkElmtNodesNum(e);
        nDofsPerNode=nDofs/nNodes;

        // for the disp and velocity in current time step
        VecGetValues(_Useq,nDofs,_elDofs.data(),_elU.data());
		VecGetValues(_Uoseq, nDofs, _elDofs.data(), _elUo.data());// XY
		VecGetValues(_U_ttseq, nDofs, _elDofs.data(), _elUtt.data());
        VecGetValues(_Vseq,nDofs,_elDofs.data(),_elV.data());
        // for the disp and velocity in the previous time step
        VecGetValues(_Uoldseq,nDofs,_elDofs.data(),_elUold.data());
        VecGetValues(_Voldseq,nDofs,_elDofs.data(),_elVold.data());
  
		//XY 
		//扣除平移[相对1st节点]
		CalctranU(nNodes, nDim, nDofs, Xe, Ua, Ub, tranU);
		//扣除面外和面内转动

		//根据网格类型计算.
		if (((emeshtype == MeshType::EDGE2) || (emeshtype == MeshType::TRUSS2) || (emeshtype == MeshType::BEAM2) || (emeshtype == MeshType::PIPE2)) && (nNodes == 2)) {
			//需要的是单元类型meshtype->rod  tri quad ...
			for (i = 1; i <= nNodes; i++) {

			}
		}
		else if ((emeshtype == MeshType::TRI3) && (nNodes == 3)) {
			//面外法向转动theta_op   out-plane rotate

			for (i = 1; i <= nNodes; i++) {


			}

			//面内平行转动theta_ip   in-plane rotate



		}
		else if ((emeshtype == MeshType::QUAD4) && (nNodes == 4)) {

		}
		else if ((emeshtype == MeshType::TET4) && (nNodes == 4)) {

		}
		else if ((emeshtype == MeshType::HEX8) && (nNodes == 8)) {

		}
		else {
			snprintf(buff2, 70, "+++ Unsupported MeshType with nodes =%2d for pure deformation +++!", nNodes);//
			str2 = buff2;
			MessagePrinter::PrintErrorTxt(str2);
			MessagePrinter::AsFem_Exit();
		}

        if(calctype==FECalcType::ComputeResidual){
            fill(_R.begin(),_R.end(),0.0);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            fill(_K.begin(),_K.end(),0.0);
        }
        else if(calctype==FECalcType::Projection){
            for(auto &it:_gpProj) it.second=0.0;
        }
		else if (calctype == FECalcType::InitHistoryVariable || calctype == FECalcType::UpdateHistoryVariable) {//this else by XY add
			fill(_gpHist.begin(), _gpHist.end(), 0.0);
			if (calctype == FECalcType::InitHistoryVariable) {
				fill(_M.begin(), _M.end(), 0.0);
				damp_factor = mesh.GetMateListPtr()[0].damp_alpha;
				fill(_C.begin(), _C.end(), damp_factor);//实际上求解阻尼可以不是这个材料的阻尼
				//_localM.setZero(); _localC.setZero(); 这应该在节点循环里  实际这些都没用 组装成Mass Cdamp就行，以后不变
			}
		}

        xi=0.0;eta=0.0;zeta=0.0;DetJac=1.0;w=1.0;
        elVolume=0.0;
        for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
            // init all the local K&R array/matrix
            // get local history(old) value on each gauss point
            if(calctype!=FECalcType::InitMaterialAndProjection){
                mateSystem.GetMaterialsOldPtr().ScalarMaterials=solutionSystem._ScalarMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetMaterialsOldPtr().VectorMaterials=solutionSystem._VectorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetMaterialsOldPtr().Rank2Materials=solutionSystem._Rank2TensorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetMaterialsOldPtr().Rank4Materials=solutionSystem._Rank4TensorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
            }
            // calculate the current shape funs on each gauss point
            if(nDim==1){
                w =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                fe._BulkShp.Calc(xi,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            else if(nDim==2){
                w  =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                eta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                fe._BulkShp.Calc(xi,eta,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            else if(nDim==3){
                w   =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi  =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                eta =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                zeta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,3);
                fe._BulkShp.Calc(xi,eta,zeta,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            JxW=w*DetJac;//w refers weight高斯点权重系数
            elVolume+=1.0*JxW;//这是计算单元体积(一维长度，二维面积，三维体积)吗,高斯点权重求和。
            // calculate the coordinate of current gauss point
            _gpCoord(1)=0.0;_gpCoord(2)=0.0;_gpCoord(3)=0.0;
            for(i=1;i<=nNodes;++i){
                _gpCoord(1)+=_elNodes(i,1)*fe._BulkShp.shape_value(i);
                _gpCoord(2)+=_elNodes(i,2)*fe._BulkShp.shape_value(i);
                _gpCoord(3)+=_elNodes(i,3)*fe._BulkShp.shape_value(i);
            }
            
            
            if(calctype==FECalcType::ComputeResidual){
                _localR.setZero();
            }
            else if(calctype==FECalcType::ComputeJacobian){
                _localK.setZero();
            }
            else if(calctype==FECalcType::Projection){
                for(auto &it:_gpProj) it.second=0.0;
            }
            else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
                for(auto &it:solutionSystem._ScalarMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._VectorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._Rank2TensorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._Rank4TensorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
            }
            // now we do the loop for local element, *local element could have multiple contributors according
            // to your model, i.e. one element (or one domain) can be assigned by multiple [elmt] sub block in your input file !!!
            for(int ielmt=1;ielmt<=static_cast<int>(dofHandler.GetIthElmtElmtMateTypePair(e).size());ielmt++){//比如铺层材料?  //其实是该单元处还有次级单元，如梁单元，共用节点
                elmttype=dofHandler.GetIthElmtJthKernelElmtType(e,ielmt);
                matetype=dofHandler.GetIthElmtJthKernelMateType(e,ielmt);
                localDofIndex=dofHandler.GetIthBulkElmtJthKernelDofIndex(e,ielmt);
                mateindex=dofHandler.GetIthBulkElmtJthKernelMateIndex(e,ielmt);
                nDofsPerSubElmt=static_cast<int>(localDofIndex.size());

                // now we calculate the local dofs and their derivatives
                // *this is only the local one, which means, i.e., if current element use dofs=u v
                // then we only calculate u v and their derivatives on each gauss point,
                // then for the next loop, the same element may use 'dofs=u v w', then we will calculate
                // u v w and their derivatives, and so on!!!
                // In short, here we dont offer the localK and localR for the whole element, instead, the quantities
                // of a single gauss point according to each sub [elmt] block
                for(j=1;j<=nDofsPerSubElmt;j++){
                    // !!!: the index starts from 1, not 0, please following the same way in your UEL !!!
                    _gpU[j]=0.0;_gpV[j]=0.0; _gpUtt[j] = 0.0; _gpUOld[j]=0.0;_gpVOld[j]=0.0;
                    _gpGradU[j](1)=0.0;_gpGradU[j](2)=0.0;_gpGradU[j](3)=0.0;
                    _gpGradUOld[j](1)=0.0;_gpGradUOld[j](2)=0.0;_gpGradUOld[j](3)=0.0;
                    _gpGradV[j](1)=0.0;_gpGradV[j](2)=0.0;_gpGradV[j](3)=0.0;
                    _gpGradVOld[j](1)=0.0;_gpGradVOld[j](2)=0.0;_gpGradVOld[j](3)=0.0;
					_gpGradUtt[j](1) = 0.0; _gpGradUtt[j](2) = 0.0; _gpGradUtt[j](3) = 0.0;
                    jj=localDofIndex[j-1];
                    for(i=1;i<=nNodes;++i){
                        _gpU[j]+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);
                        _gpUOld[j]+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);

                        _gpV[j]+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);
                        _gpVOld[j]+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);

						_gpUtt[j] += _elUtt[(i - 1)*nDofsPerNode + jj - 1] * fe._BulkShp.shape_value(i);// XY

                        _gpGradU[j](1)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradU[j](2)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradU[j](3)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradUOld[j](1)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradUOld[j](2)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradUOld[j](3)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradV[j](1)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradV[j](2)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradV[j](3)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradVOld[j](1)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradVOld[j](2)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradVOld[j](3)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

						_gpGradUtt[j](1) += _elUtt[(i - 1)*nDofsPerNode + jj - 1] * fe._BulkShp.shape_grad(i)(1);// XY
						_gpGradUtt[j](2) += _elUtt[(i - 1)*nDofsPerNode + jj - 1] * fe._BulkShp.shape_grad(i)(2);
						_gpGradUtt[j](3) += _elUtt[(i - 1)*nDofsPerNode + jj - 1] * fe._BulkShp.shape_grad(i)(3);

                    }
                }

                if(calctype==FECalcType::ComputeResidual){
                    _subR.setZero();
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    _subK.setZero();
                }
                //*****************************************************
                //*** For user material calculation(UMAT)
                //*****************************************************
                if(calctype==FECalcType::InitHistoryVariable){
                    mateSystem.InitBulkMateLibs(matetype,mateindex,nDim,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV);
                }
                else{
                    mateSystem.RunBulkMateLibs(matetype,mateindex,nDim,t,dt,_gpCoord,_gpU,_gpUOld,_gpV,_gpVOld,
                                               _gpGradU,_gpGradUOld,_gpGradV,_gpGradVOld);//查V elV gpV等 都加入Utt量计算
                }
                //*****************************************************
                //*** For user element calculation(UEL)
                //*****************************************************
                if(calctype==FECalcType::ComputeResidual){
                    for(i=1;i<=nNodes;i++){
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,nDofsPerSubElmt,t,dt,ctan,
                            _gpCoord,_gpU,_gpUOld,_gpV,_gpVOld,_gpGradU,_gpGradUOld,_gpGradV,_gpGradVOld,
                            fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(i),// for Residual, we only need test fun
                            fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(i),
                            mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),
                            _gpProj,_subK,_subR);
                        AssembleSubResidualToLocalResidual(nDofsPerNode,nDofsPerSubElmt,i,_subR,_localR);
                    }
                }
                else if(calctype==FECalcType::ComputeJacobian){//XY all the ComputeJacobian are not needed, later to comment away!
                    for(i=1;i<=nNodes;i++){
                        for(j=1;j<=nNodes;j++){
                            elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,nDofsPerSubElmt,t,dt,ctan,
                            _gpCoord,_gpU,_gpUOld,_gpV,_gpVOld,_gpGradU,_gpGradUOld,_gpGradV,_gpGradVOld,
                            fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(j),
                            fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(j),
                            mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),
                            _gpProj,_subK,_subR);
                            AssembleSubJacobianToLocalJacobian(nDofsPerNode,i,j,_subK,_localK);
                        }
                    }
                }
                else if(calctype==FECalcType::Projection){
                    for(i=1;i<=nNodes;i++){
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,nDofsPerSubElmt,t,dt,ctan,
                        _gpCoord,_gpU,_gpUOld,_gpV,_gpVOld,_gpGradU,_gpGradUOld,_gpGradV,_gpGradVOld,
                        fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(i),// for Residual, we only need test fun
                        fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(i),
                        mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),
                        _gpProj,_subK,_subR);
                    }
                    // here we should not assemble the local projection, because the JxW should not be accumulated
                    // inside the element-loop, but the gpProj should be.
                    // therefore, each sub element should use its own place of gpProj, in short, the gpProj is shared
                    // between different elements
                }
            }//=====> end-of-sub-element-loop

            //***********************************************
            //*** accumulate all the local contribution inside gauss loop
            if(calctype==FECalcType::ComputeResidual){
                //AccumulateLocalResidual(nDofs,_elDofsActiveFlag,JxW,_localR,_R);
				AccumulateLocalResidual(nDofs, _elDofsActiveFlag, JxW, elsect, _localR, _R);//shell element with t to get Rint,通过各高斯点累加   ==这里R改为内力不是合力 还有其他影响吗？如节点力插值
				//这里根据约束条件计算支反力Rconst

            }
            else if(calctype==FECalcType::ComputeJacobian){
                //AccumulateLocalJacobian(nDofs,_elDofsActiveFlag,JxW,_localK,_K);
            }
            else if(calctype==FECalcType::Projection){
				// △_sigma_local to global via Q^T, then somewhere later with sigma at time to rotate R together
				TransformLocalDletaSigmaToGlobal(Ql, mateSystem.GetRank2MatePtr());// tb时刻虚拟构形下的全局坐标系下的应力增量.
				RenewReacforceViaSect(elsect, _gpProj);
				AssembleLocalProjectionToGlobal(nNodes, JxW, fe._BulkShp, _gpProj,//proj里含有reaction force,要用截面secA
												mateSystem.GetScalarMatePtr(),
                                                mateSystem.GetVectorMatePtr(),
                                                mateSystem.GetRank2MatePtr(),
                                                mateSystem.GetRank4MatePtr(),
                                                solutionSystem);
            }
            else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
                AssembleSubHistToLocal(e,fe._BulkQPoint.GetQpPointsNum(),gpInd,mateSystem.GetMaterialsPtr(),solutionSystem);
            }
        }//----->end of gauss point loop
		if (calctype == FECalcType::InitHistoryVariable) {
			mesh.SetBulkMeshIthBulkElmtVolume(e, elVolume);//自己编的在AbaqusIOReadMesh文件里计算了单元体积.... 这里重复了
			_BulkVolumes += elVolume;
		}

        
        if(calctype==FECalcType::ComputeResidual){
			//单元内力需要进一步通过合力和力矩平衡确定简化掉的自由度上的分力.
			//三角形薄膜单元内力分量更新
			RenewRcomponts(nDofsPerNode, Xlocal, _R);
			TransformLocalRToGlobalR(nNodes, nDim, nDofs, Atheta, Ql, totalrotate, _R);
			AssembleLocalResidualToGlobalResidual(nDofs,_elDofs,_R,RHS);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            //AssembleLocalJacobianToGlobalJacobian(nDofs,_elDofs,_K,AMATRIX);
        }
        else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
            AssembleLocalHistToGlobal(e,fe._BulkQPoint.GetQpPointsNum(),solutionSystem);
        }
    }//------>end of element loop

    //********************************************************************
    //*** finish all the final assemble for different matrix and array
    //********************************************************************
    if(calctype==FECalcType::ComputeResidual){
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
    }
    else if(calctype==FECalcType::Projection){
        Projection(mesh.GetBulkMeshNodesNum(),solutionSystem);
    }
	else if (calctype == FECalcType::InitHistoryVariable || calctype == FECalcType::UpdateHistoryVariable) {
		VecAssemblyBegin(solutionSystem._Hist);
		VecAssemblyEnd(solutionSystem._Hist);
		VecCopy(solutionSystem._Hist, solutionSystem._HistOld);

		if (calctype == FECalcType::InitHistoryVariable) {
			nodalDofFlag = dofHandler.GetNodalDofFlag();
			if (nodalDofFlag.size() != nodeeqVolume.size()) {
				MessagePrinter::PrintStars();
				MessagePrinter::PrintErrorTxt("  nodalDofFlag.size()!= nodeeqVolume.size() ! ");
				//MessagePrinter::AsFem_Exit();//测试node数与等效节点质量数相等否
			}

			nActiveDofs = 0;
			for (i = 0; i < nodalDofFlag.size(); i++) {//total nNodes
				inodeMass = rho * elsect*nodeeqVolume[i];
				这里4 - 6dof应该是惯性矩计算;
				for (j = 0; j < nodalDofFlag[i].size(); j++) {//nDofsPerNode
					if (nodalDofFlag[i][j] > -10.0) {//只考虑有activedof>0.0节点质量  不合理，得避免m分母为0,改为>-10，所有mass都考虑
						//MASS[nActiveDofs] = inodeMass;// nodeeqMass[i];
						VecSetValue(MASS, nActiveDofs, inodeMass, INSERT_VALUES);
						nActiveDofs += 1;
					}
				}
			}

			VecAssemblyBegin(MASS);
			VecAssemblyEnd(MASS);
			VecAssemblyBegin(CDAMP);
			VecAssemblyEnd(CDAMP);
		}
	}



    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecScatterDestroy(&_scatteruold);
    VecScatterDestroy(&_scatterv);
    VecScatterDestroy(&_scattervold);
	VecScatterDestroy(&_scatterutt);
	//VecScatterDestroy(&_scatterproj);
	//VecScatterDestroy(&_scatterhist);
	//VecScatterDestroy(&_scatterhistold);

    VecDestroy(&_Useq);
    VecDestroy(&_Uoldseq);

    VecDestroy(&_Vseq);
    VecDestroy(&_Voldseq);

	VecDestroy(&_U_ttseq);
	//VecDestroy(&_ProjSeq);
	//VecDestroy(&_HistSeq);
	//VecDestroy(&_HistOldSeq);
}

PetscReal FESystem::Initdtcr(Mesh &mesh) {

	//initial critical time step with 0.9 reduced factor
	vector<double> Lc, Ae;
	PetscReal dt_cr0, ltemp, li, E, rho, nu, cspeed, fdamp, damp;
	MeshType meshtype;
	//计算单元最小长度  目前考虑bar/line tri quad中的一种
	Ae = mesh.GetBulkMeshElmtVolumeListPtr();
	E = mesh.GetMateListPtr()[0].Youngs_modulus;
	rho = mesh.GetMateListPtr()[0].density;
	nu = mesh.GetMateListPtr()[0].poisson;
	damp = mesh.GetMateListPtr()[0].damp_alpha;
	fdamp = sqrt(damp*damp + 1) - damp;
	Lc.clear();
	vector<int>	nID;
	int eID, ni;
	vector<Vector3d> elNcoords;
	int _nNodesPerBulkElmt = mesh.GetBulkMeshNodesNumPerBulkElmt;
	nID.resize(_nNodesPerBulkElmt, 0);
	elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
	meshtype = mesh.GetBulkMeshBulkElmtType;
	for (int e = 1; e <= mesh.GetBulkMeshElmtsNum; ++e) {
		eID = mesh.GetBulkMeshIthElmtPhyID[e];//单元号
		for (ni = 1; ni <= _nNodesPerBulkElmt; ni++) {//bar 1,2 quad4 ni=1,2,3,4
			nID[ni - 1] = mesh.GetBulkMeshElmtConnPtr()[eID - 1][ni];
			for (int k = 1; k <= 3; k++) {//如果是二维则z=0
				elNcoords[ni - 1](k) = mesh.GetBulkMeshNodeCoordsPtr()[(nID[ni - 1] - 1) * 3 + k - 1];//elmt各点坐标xyz
			}
		}
		if (meshtype == MeshType::EDGE2 || meshtype == MeshType::BEAM2 || meshtype == MeshType::TRUSS2 || meshtype == MeshType::PIPE2) {
			li = (elNcoords[1] - elNcoords[0]).norm();
			ltemp = li;
			cspeed = sqrt(E / rho);
			//Lc.push_back(ltemp);
		}
		else if (meshtype == MeshType::TRI3) {
			ltemp = 0.0;
			li = (elNcoords[2] - elNcoords[0]).norm();
			for (ni = 1; ni <= _nNodesPerBulkElmt - 1; ni++) {
				ltemp = max(ltemp, (elNcoords[ni] - elNcoords[ni - 1]).norm());
			}
			ltemp = max(ltemp, li);
			cspeed = sqrt(E / (rho*(1 - nu * nu)));
			ltemp = (1.0 + 1.0)*Ae[e - 1] / ltemp;
		}
		else if (meshtype == MeshType::QUAD4) {
			ltemp = 0.0;
			li = (elNcoords[3] - elNcoords[0]).norm();
			for (ni = 1; ni <= _nNodesPerBulkElmt - 1; ni++) {
				ltemp = max(ltemp, (elNcoords[ni] - elNcoords[ni - 1]).norm());
			}
			ltemp = max(ltemp, li);
			cspeed = sqrt(E / (rho*(1 - nu * nu)));
			ltemp = (1.0 + 0.0)*Ae[e - 1] / ltemp;//LS-DYNA manual 19.3
		}
		else if (meshtype == MeshType::TET4) {}
		else if (meshtype == MeshType::HEX8) {}
		Lc.push_back(ltemp / cspeed);
	}
	dt_cr0 = *min_element(Lc.begin(), Lc.end());//逐单元法 还有其他大时间步的方法;后面可以类似ABAQUS输出10个单元最小步长
	dt_cr0 = dt_cr0 * 0.9*fdamp;//初始缩减系数0.9

	return dt_cr0;
}

void FESystem::CalcRF(DofHandler &dofHandler, Vec &REXT, Vec &RINT, Vec &RHS, Vec &RF){
	VecSet(RHS, 0);
	//VecAXPY(RHS, 1.0, REXT);
	//VecAXPY(RHS, 1.0, RINT);//compute total force without RF<--constrained force.
	//VecAXPBYPCZ(RHS, 1.0, 1.0, 0.0, REXT, RINT);
	VecWAXPY(RHS, 1.0, REXT, RINT);//compute total force without RF<--constrained force.
	int i,j,nDofs;
	vector<vector<double>> nodalDofFlag;
	nodalDofFlag = dofHandler.GetNodalDofFlag();
	VecZeroEntries(RF);// VecSet(RF,0);
	VecCopy(RHS, RF);
	nDofs = 0;
	for (i = 0; i < nodalDofFlag.size(); i++) {//total nNodes
		for (j = 0; j < nodalDofFlag[i].size(); j++) {//nDofsPerNode
			if (nodalDofFlag[i][j] > 0.0) {//只考虑有activedof支反力
				VecSetValue(RF, nDofs, 0.0, INSERT_VALUES);
				nDofs += 1;
			}
			else if (nodalDofFlag[i][j] == 0.0) {			
				VecSetValue(RHS, nDofs, 0.0, INSERT_VALUES);
				nDofs += 1;
			}			
		}
	}
	VecScale(RF, -1.0);
	//VecAXPY(RHS, 1.0, RF);//上面已经把约束算上，deactivedof的力为0.

	VecAssemblyBegin(RHS);
	VecAssemblyEnd(RHS);
	VecAssemblyBegin(RF);
	VecAssemblyEnd(RF);
}

void FESystem::CalcreRF(Vec &RF, Vec &reRF) {
	VecSet(reRF, 0);
	VecAXPY(reRF, -1.0, RF);//BC beared force==total force at constraint points.
}

void FESystem::CalcUtt(Vec &M, Vec &RHS, Vec &Utt) {
	VecPointwiseDivide(Utt, RHS, M);//Computes the componentwise division w = x/y.
}

void FESystem::InitUold(const double &dt, Vec &CDAMP, Vec &Uold, Vec &U, Vec &V, Vec &Utt) {
	double c1, psi;
	//psi = CDAMP[0];//不行
	//PetscInt ix[2] = {0,1};
	VecAssemblyBegin(CDAMP);
	VecAssemblyEnd(CDAMP);
	VecGetValues(CDAMP, 1, 0, &psi);//目前没找到直接获取一个值的函数,开始用的psi[2]...
	c1 = 1.0 + 0.5*dt*psi;
	VecCopy(U, Uold);
	VecAXPBYPCZ(Uold, -c1*dt, 0.5*dt*dt, 1.0, V, Utt);
}

void FESystem::SolveUnew(const double &dt, Vec &CDAMP, Vec &Uold, Vec &U, Vec &Unew, Vec &Utt) {
	double c1,c2, psi;
	VecAssemblyBegin(CDAMP);
	VecAssemblyEnd(CDAMP);
	VecGetValues(CDAMP, 1, 0, &psi);//目前没找到直接获取一个值的函数,开始用的psi[2]...
	c1 = 1.0 + 0.5*dt*psi;
	c2 = 1.0 - 0.5*dt*psi;
	VecCopy(U, Unew);
	VecAXPBYPCZ(Unew, -c2/c1, dt*dt/c1, 2.0/c1, Uold, Utt);
}

void FESystem::CalctranU(const int &nNodes, const int &nDim, const int &nDofs, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Ub, vector<Vector3d> &tranU) {
	tranU.resize(nNodes, 0.0);
	Xe.resize(nNodes, 0.0);
	Ua.resize(nNodes, 0.0);
	Ub.resize(nNodes, 0.0);

	//tranU[0] = 0.0; tranU[1] = 0.0; tranU[2] = 0.0; //node1 xyz
	//Xe[0](1) = _elNodes(1, 1); 
	//Xe[0](2) = _elNodes(1, 2); 
	//Xe[0](3) = _elNodes(1, 3);
	//Ua[0](1) = _elUo[0]; 
	//Ua[0](2) = _elUo[1]; 
	//Ua[0](3) = _elUo[2];
	//Ub[0](1) = _elU[0]; 
	//Ub[0](2) = _elU[1]; 
	//Ub[0](3) = _elU[2];
	//tranU[0](1) = 0.0; tranU[0](2) = 0.0; tranU[0](3) = 0.0; //node1 xyz
	for (int ie = 0; ie < nNodes; ie++) {
		//Xe[ie](1) = _elNodes(ie + 1, 1);
		//Xe[ie](2) = _elNodes(ie + 1, 2);
		//Xe[ie](3) = _elNodes(ie + 1, 3);

		//Ua[ie](1) = _elUo[nDofs*ie];
		//Ua[ie](2) = _elUo[nDofs*ie + 1];
		//Ua[ie](3) = _elUo[nDofs*ie + 2];
		//Ub[ie](1) = _elU[nDofs*ie];
		//Ub[ie](2) = _elU[nDofs*ie + 1];
		//Ub[ie](3) = _elU[nDofs*ie + 2];

		//tranU[ie](1) = _elU[nDofs*ie] - _elU[0];
		//tranU[ie](2) = _elU[nDofs*ie + 1] - _elU[1];
		//tranU[ie](3) = _elU[nDofs*ie + 2] - _elU[2];

		for (int xi = 1; xi <= nDim; xi++) {//note: 1d point element will error!
			Xe[ie](xi) = _elNodes(ie + 1, xi);
			Ua[ie](xi) = _elUo[nDofs*ie + xi - 1];
			Ub[ie](xi) = _elU[nDofs*ie + xi - 1];
			tranU[ie](xi) = (_elU[nDofs*ie + xi - 1] - _elU[xi - 1]) - (_elUo[nDofs*ie + xi - 1]- _elUo[xi - 1]);//b时刻相对于a时刻 相对位移增量
		}
	}
}

void FESystem::CalcOutPlaneRotate(const int &nNodes, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Ub, vector<Vector3d> &Xai1, vector<Vector3d> &Xbi1, vector<Vector3d> &eai1, vector<Vector3d> &ebi1, vector<double> &Xai1norm, vector<double> &Xbi1norm, Vector3d &thetaout, Vector3d &nout, vector<Vector3d> &Utout) {
	// XY	
	double theta1=0.0;

	Xai1.resize(nNodes, 0.0);//Xa11=0  为了编号对应，保留11编号.
	Xbi1.resize(nNodes, 0.0);
	Xai1norm.resize(nNodes, 0.0);
	Xbi1norm.resize(nNodes, 0.0);
	eai1.resize(nNodes, 0.0);//ea11=0
	ebi1.resize(nNodes, 0.0);
	for (int ie = 1; ie < nNodes; ie++) {//start from 2nd node
		Xai1[ie] = (Xe[ie]+ Ua[ie])-(Xe[0] + Ua[0]);
		Xai1norm[ie] = Xai1[ie].norm();
		eai1[ie] = Xai1[ie]/ Xai1norm[ie];

		Xbi1[ie] = (Xe[ie] + Ub[ie]) - (Xe[0] + Ub[0]);
		Xbi1norm[ie] = Xbi1[ie].norm();
		ebi1[ie] = Xbi1[ie] / Xbi1norm[ie];
	}
	na = eai1[1] ^ eai1[2];
	nb = ebi1[1] ^ ebi1[2];
	theta1 = asin((na^nb).norm()/(na.norm()*nb.norm()));
	nout = (na^nb)/((na^nb).norm());
	thetaout = nout*theta1;	
}

void FESystem::CalcInPlaneRotate(const int &nNodes, const int &nDim, vector<Vector3d> &Xai1, vector<Vector3d> &Xbi1, Vector3d &na, Vector3d &thetain){
	Vector3d Xac(0.0), Xbc(0.0);
	vector<Vector3d> Xaic, Xbic, Xpic;//ta tb各点相对于C坐标    (b->a)投影点相对于Ca坐标
	vector<double> thetai(nNodes,0.0);
	double theta2=0.0;

	Xaic.resize(nNodes, 0.0);
	Xbic.resize(nNodes, 0.0);
	Xpic.resize(nNodes, 0.0);
	for (int xi = 1; xi <= nDim; xi++) {
		for (int ie = 0; ie < nNodes; ie++) {//相对于1st node的形心C坐标
			Xac(xi) += 1.0/nNodes *(Xai1[ie](xi));
			//Xac(2) += 1.0 / 3.0*;
			//Xac(3) += 1.0 / 3.0*;
			Xbc(xi) += 1.0/nNodes *(Xbi1[ie](xi));
		}
	}

	for (int ie = 0; ie < nNodes; ie++) {//各点相对于C坐标
		Xaic[ie] = Xai1[ie]-Xac;
		Xbic[ie] = Xbi1[ie] - Xbc;
		Xpic[ie] = Xbic[ie]-na*(Xbic[ie].ODot(na));//投影点
		thetai[ie] = asin(na*((Xaic[ie] ^ Xpic[ie])/(Xaic[ie] ^Xpic[ie]).norm()));
	}
	theta2 = std::accumulate(std::begin(thetai), std::end(thetai), 0.0);
	theta2 = theta2 / nNodes;
	thetain = na * theta2;
}

void FESystem::SolveUdeform(const int &nNodes, vector<Vector3d> &tranU, vector<Vector3d> &Xbi1, Vector3d &nrotate, double &totalrotate, vector<Vector3d> &Ud) {
	Ud.resize(nNodes, 0.0);
	MatrixXd Atheta(3,3,0.0);
	Atheta(1, 2)= -nrotate(3);
	Atheta(1, 3) = nrotate(2);
	Atheta(2, 1) = nrotate(3);
	Atheta(2, 3) =-nrotate(1);
	Atheta(3, 1) =-nrotate(2);
	Atheta(3, 2) = nrotate(1);

	//vector<VectorXd> XXbi1;
	//XXbi1.resize(nNodes, 0.0);
	//for (int vi = 0; vi < nNodes; vi++) {// transverse from Vector3d to VectorXd, for MatrixXd operation, but later reverse action again needed, it's better to add MatrixXd*Vector3d operation.
	//	XXbi1[vi].Resize(3, 0.0);
	//	for (int vii = 0; vii < 3; vii++) {//可以用copy
	//		XXbi1[vi](vii) = Xbi1[vi](vii);
	//	}
	//}

	for (int ie = 1; ie < nNodes; ie++) {
		Ud[ie] = tranU[ie]+ (Atheta*Atheta*(1-cos(-totalrotate))+Atheta*sin(-totalrotate))*Xbi1[ie];
	}

}

void FESystem::SolveUlocal(const int &nNodes, const int &nDim, const int &nDofs, vector<Vector3d> &Xe, vector<Vector3d> &Ua, vector<Vector3d> &Xbi1, vector<Vector3d> &Ud, Vector3d &na, Vector3d &nrotate, vector<Vector3d> &Xlocal, vector<Vector3d> &Ulocal, vector<double> &elUlocal) {
	MatrixXd Ql(nDim, nDim,0.0);
	Xlocal.resize(nNodes, 0.0);
	Ulocal.resize(nNodes, 0.0);
	elUlocal.resize(nNodes*nDofs, 0.0);
	Vector3d ex, ey, ez,tempex(0.0);
	ez = na;
	if (Ud[1].norm()>1.0e-15) {
		ex = Ud[1] / (Ud[1].norm());//(2．76)公式后面解释当节点2变形很小，直接用1a2a边为ex方向
	}
	else {
		tempex = Ud[1]+ (Xe[1] + Ua[1]) - (Xe[0] + Ua[0]);
		ex = tempex / (tempex.norm());//Xa2=(Xe[1] + Ua[1])       Xa1=(Xe[0] + Ua[0])
	}
	ey = (ez^ex) / (ez^ex).norm;
	for (int dim = 1; dim <= nDim; dim++) {//坐标转换矩阵
		Ql(nDim,dim)=ez(dim);//ez
		Ql(nDim-1, dim) = ey(dim);//nDim=3, 3d才有
		Ql(1, dim) = ex(dim);
	}
	for (int i = 0; i < nNodes; i++) {
		Xlocal[i] = Ql * Xbi1[i];//(2．77)
		Ulocal[i] = Ql * Ud[i];//(2．78)
	}
	//薄膜三角形单元 u1=v1=w1=v2=w2=0，w3=0. 如果e2设在2节点方向上，只剩三个独立非零位移量(u2,u3,v3)
	//避免系统误差，强行赋值
	Ulocal[1](2) = 0.0; Ulocal[1](3)=0.0; Ulocal[2](3) = 0.0;

	for (int ie = 0; ie < nNodes; ie++) {
		for (int xi = 1; xi <= nDim; xi++) {//note: 1d point element will error!
			//Ub[ie](xi) = _elU[nDofs*ie + xi - 1];
			elUlocal[nDofs*ie + xi - 1] = Ulocal[ie](xi);//反操作，将转换的位移1-3dof还原到类似原_elU向量里.  没赋转角值4-6dof
		}
	}


}

void FESystem::RenewRcomponts(const int &nDofsPerNode, const vector<Vector3d> &Xlocal, vector<double> &_R) {
	double x2, x3, y2, y3, f2x, f2y, f3x, f3y;
	double f1x, f1y, f2y;
	x2 = Xlocal[1](1);
	x3 = Xlocal[2](1);
	y2 = Xlocal[1](2);
	y3 = Xlocal[2](2);
	f2x = _R[1 * nDofsPerNode];
	f3x = _R[2 * nDofsPerNode];
	f3y = _R[2 * nDofsPerNode+1];

	f2y = (f3x*y3+f2x*y2-f3y*x3)/x2;	//	sumM1=0
	f1x = -(f2x+f3x);					//	Fx=0
	f1y = -(f2y+f3y);					//	Fy=0

	_R[0] = f1x;
	_R[1] = f1y;
	_R[1*nDofsPerNode+1] = f2y;
}

void FESystem::RenewReacforceViaSect(const double &secA, map<string, double> &gpProj) {
	gpProj["reacforce_x"] *= secA;
	gpProj["reacforce_y"] *= secA;
	gpProj["reacforce_z"] *= secA;
}

void FESystem::TransformLocalRToGlobalR(const int &nNodes, const int &nDim, const int &nDofs, const MatrixXd &Atheta, const MatrixXd &Ql, const double &totalrotate, vector<double> &_R) {
	MatrixXd Qt;
	vector<Vector3d> elR, elRbar;
	elR.resize(nNodes,0.0);
	elRbar.resize(nNodes, 0.0);
	for (int ie = 0; ie < nNodes; ie++) {
		for (int xi = 1; xi <= nDim; xi++) {//note: 1d point element will error!
			elR[ie](xi) = _R[nDofs*ie+xi-1];//ndofspernode
		}
	}
	Qt=Ql.Transpose();
	for (int ie = 1; ie < nNodes; ie++) {
		elRbar[ie] = Qt*elR[ie];
		elRbar[ie] = (Atheta*Atheta*(1 - cos(totalrotate)) + Atheta * sin(totalrotate))*elRbar[ie];
	}

	for (int ie = 0; ie < nNodes; ie++) {
		for (int xi = 1; xi <= nDim; xi++) {
			_R[nDofs*ie+xi-1] = elR[ie](xi);
		}
	}

}

//void FESystem::TransformLocalDletaSigmaToGlobal(const int &nNodes, const int &nDim, const int &nDofs, const MatrixXd &Ql, Rank2MateType &Rank2Materials) {
void FESystem::TransformLocalDletaSigmaToGlobal(const MatrixXd &Ql, Rank2MateType &Rank2Materials) {
	//MatrixXd Qt;
	//Qt = Ql.Transpose();
	Rank2Materials.at("stress") = Rank2Materials.at("stress").TransposeClampMultiply(Ql);
}
