//****************************************************************

//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.14
//+++ Purpose: here we apply the load conditions we defined in
//+++          each [loads] sub blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "LoadSystem/LoadSystem.h"
#include "DofHandler/DofHandler.h"

void LoadSystem::ApplyLoad(const Mesh &mesh, const DofHandler &dofHandler, FE &fe, const FECalcType &calctype, const double &t, Vec &U, Vec &REXT){
	//Amp t ele-node-dof value matrix
	//�����ȼ����ampvalue ����ampvalue������Ҫampblock��t

	for (auto it : _LoadBlockList) {
		if (it._LoadkeyName == "Cload") {//Cload		1-6
			ApplyCload(mesh, dofHandler, calctype, t, loadblock, REXT);
		}
		else if (it._LoadkeyName == "Dload") {//Dload		GRAV/BX BY BZ
			ApplyDload(mesh, dofHandler, calctype, t, loadblock, REXT);
		}
		else if (it._LoadkeyName == "Dsload") {//Dsload	TRSHR/TRVEC/P
			ApplyDsload(mesh, dofHandler, calctype, t, loadblock, REXT);
		}
	}
}
//****************************************************
void LoadSystem::ApplyInitialLoad(const Mesh &mesh, const DofHandler &dofHandler, const double &t, Vec &U){

}

void LoadSystem::ApplyCload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &REXT) {
	string _LoadAmpName, _LoadSetName;

	vector<int> nodesNo;
	vector<string> sDof;
	vector<double> cvalues;
	int iamp;
	double ampvalue, loadvalue;
	vector<pair<string, vector<int>>> _NodeSetName2NodeIDsList;
	vector<int>    _Nodeids;
	int			   dofID, _IthNodeJthDofIndex;

	vector<string> _DofsName;//"1-6"
	vector<double> _LoadValues;
	bool           _IsFollower;

	_LoadAmpName = loadblock._LoadAmpName;
	_LoadSetName = loadblock._LoadSetName;
	_DofsName    = loadblock._DofsName;
	_LoadValues  = loadblock._LoadValues;

	if (loadblock._Partunit.nodesNo.size() > 0) {
		nodesNo = loadblock._Partunit.nodesNo;
		sDof = loadblock._Partunit.sDof;
		cvalues = loadblock._Partunit.values;
	}
	_IsFollower  = loadblock._IsFollower;
	_normals = 0.0;

	iamp = GetAmpIndexFromName(_LoadAmpName);
	ComputeiAmp(iamp, t, ampvalue);
	if (_LoadSetName.size() > 0 && nodesNo.size() == 0) {//nodeset input
		_NodeSetName2IDsList = mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr();
		for (int i = 0; i < _NodeSetName2IDList.size(), ++i) {
			if (_LoadSetName == _NodeSetName2IDList[i].first) {
				_Nodeids = _NodeSetName2IDList[i].second;
			}
		}
	}
	else if (_LoadSetName.size() == 0 && nodesNo.size() > 0) {//dicretize node input
		_Nodeids.clear();
		_Nodeids = nodesNo;
		_DofsName.clear();
		_DofsName = sDof;
		_LoadValues.clear();
		_LoadValues = cvalues;
	}

	int rankne, eStart, eEnd;
	int e, i, iInd;
	MPI_Comm_size(PETSC_COMM_WORLD, &_size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);


	rankne = _Nodeids.size() / _size;
	eStart = _rank * rankne;
	eEnd = (_rank + 1)*rankne;
	if (_rank == _size - 1) eEnd = _Nodeids.size();
	for (e = eStart; e < eEnd; ++e) {
		if (_LoadSetName.size() > 0 && nodesNo.size() == 0) {
			for (int s = 0; s < _DofsName.size(), ++s) {
				dofID = atoi(_DofsName[s].c_str());
				iInd = dofHandler.GetIthNodeJthDofIndex(_Nodeids[e], dofID) - 1;
				loadvalue = _LoadValues[s] * ampvalue;
				if (_IsFollower) {
					//�ҳ��ڵ�ת���ǣ�normal
					//����ά��ά ��Ԫ�����й� �ο�applyBC



				}//follow�����ת�����κ��xyz�����ɶȷ���

				if (calctype == FECalcType::ComputeResidual) {
					VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);//�Ƿ�follow��Ӱ��װ��λ��
				}
			}
		}
		else if (_LoadSetName.size() == 0 && nodesNo.size() > 0) {
			dofID = atoi(_DofsName[e].c_str());
			iInd = dofHandler.GetIthNodeJthDofIndex(_Nodeids[e], dofID) - 1;
			loadvalue = _LoadValues[e] * ampvalue;
			if (_IsFollower) {
				//�ҳ��ڵ�ת���ǣ�normal
				//����ά��ά ��Ԫ�����й� �ο�applyBC
				//ͬ��


			}//follow�����ת�����κ��xyz�����ɶȷ���

			if (calctype == FECalcType::ComputeResidual) {
				VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);//�Ƿ�follow��Ӱ��װ��λ��
			}
		}
	}
}

void LoadSystem::ApplyDload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &REXT) {
	string _LoadAmpName, _LoadSetName;

	vector<int> elesNo;
	vector<string> sDof;
	vector<double> dvalues;
	Vector3d _direction;
	int iamp;
	double ampvalue, loadvalue;
	double erho, evol;
	vector<pair<string, vector<int>>> _ElSetName2NodeIDsList;
	vector<int>    _Elmtids;
	int			   _nNodesPerElmt, _IthElmtJthNodeID, dofID;
	vector<double> _ElmtsVolume;

	vector<string> _DofsName;//"GRAV BX BY BZ"����
	vector<double> _LoadValues;
	//bool           _IsFollower;

	_LoadAmpName = loadblock._LoadAmpName;
	_LoadSetName = loadblock._LoadSetName;
	_DofsName = loadblock._DofsName;
	_LoadValues = loadblock._LoadValues;//����ֻ��һ�У���������������
	_direction = loadblock._direction;

	if (loadblock._Partunit.elesNo.size() > 0) {
		elesNo = loadblock._Partunit.elesNo;
		sDof = loadblock._Partunit.sDof;
		dvalues = loadblock._Partunit.values;
	}
	//_IsFollower = loadblock._IsFollower;
	//_normals = 0.0;

	iamp = GetAmpIndexFromName(_LoadAmpName);
	ComputeiAmp(iamp, t, ampvalue);
	if (_LoadSetName.size() > 0 && elesNo.size() == 0) {//eleset input
		_ElSetName2IDsList = mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr();
		for (int i = 0; i < _ElSetName2IDList.size(), ++i) {
			if (_LoadSetName == _ElSetName2IDList[i].first) {
				_Elmtids = _ElSetName2IDList[i].second;
			}
		}
		//_Elmtids = GetBulkMeshElmtIDsViaPhysicalName(_LoadSetName);
	}
	else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {//dicretize element input
		_Elmtids.clear();
		_Elmtids = elesNo;
		_DofsName.clear();
		_DofsName = sDof;
		_LoadValues.clear();
		_LoadValues = dvalues;
	}

	int rankne, eStart, eEnd;
	int e, i, j, k, iInd;//iInd ��_IthNodeJthDofIndex, 
	vector<int> iconn;
	MPI_Comm_size(PETSC_COMM_WORLD, &_size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

	rankne = _Elmtids.size() / _size;
	eStart = _rank * rankne;
	eEnd = (_rank + 1)*rankne;
	if (_rank == _size - 1) eEnd = _Elmtids.size();
	_nNodesPerElmt = mesh.GetBulkMeshNodesNumPerBulkElmt;//һ�ֵ�Ԫ
	_ElmtsVolume = mesh.GetBulkMeshElmtVolumeListPtr();
	erho = mesh.GetMateListPtr()[0].density;//һ�ֲ���
	for (e = eStart; e < eEnd; ++e) {
		evol = _ElmtsVolume[_Elmtids[e]-1];
		iconn = mesh.GetBulkMeshIthElmtNodesID(_Elmtids[e]);
		if (_LoadSetName.size() > 0 && elesNo.size() == 0) {
			//_nNodesPerElmt = GetBulkMeshIthElmtNodesNumViaPhyName(_LoadSetName,e+1);
			for(i = 1; i <= _nNodesPerElmt; ++i){			
				j = iconn[i-1];
				//_nNodesPerElmt, _IthElmtJthNodeID, _IthNodeJthDofIndex, dofID;
				for (k = 1; k <= 3; ++k) {//xyz����  ���ɻ�ȡdim���
					if (_DofsName.size() == 1) {//gravity //GRAV
						loadvalue = erho*evol*_LoadValues[0]/ _nNodesPerElmt * _direction(k)* ampvalue;//��������/�ڵ���
					}
					else if (_DofsName.size() >= 2) {//body force //Ҫת��BX -> 1  ��άû��BZ ->3
						loadvalue = evol*_LoadValues[k-1] / _nNodesPerElmt* ampvalue;//��������/�ڵ���
					}
					dofID = k;
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
		else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {//��㸳ֵ
			if (_DofsName[0].find("GRAV") != string::npos) {//_direction.normsq()>0.0
				for (i = 1; i <= _nNodesPerElmt; ++i) {
					j = iconn[i - 1];
					for (k = 1; k <= 3; ++k) {//xyz����
						loadvalue = erho * evol*_LoadValues[e] / _nNodesPerElmt * _direction(k)* ampvalue;//��������/�ڵ���
						dofID = k;
						iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
						if (calctype == FECalcType::ComputeResidual) {
							VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
						}
					}
				}
			}
			else if (_DofsName[0].find("GRAV") == string::npos) {
				loadvalue = evol * _LoadValues[e] / _nNodesPerElmt * ampvalue;//��������/�ڵ���
				for (i = 1; i <= _nNodesPerElmt; ++i) {
					j = iconn[i - 1];
					if (_DofsName[e].find("BX") != string::npos) {
						dofID = 1;
					}
					else if (_DofsName[e].find("BY") != string::npos) {
						dofID = 2;
					}
					else if (_DofsName[e].find("BZ") != string::npos) {
						dofID = 3;
					}
					else {
						MessagePrinter::PrintStars();
						MessagePrinter::PrintErrorTxt("  unexpected dofsname in element "_Elmtids[e]", dof "_DofsName[e]" ! ");
						MessagePrinter::AsFem_Exit();
					}
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
	}
}


void LoadSystem::ApplyDsload(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &t, const LoadBlock &loadblock, Vec &REXT) {
	string _LoadAmpName;
	double ampvalue;
	int iamp;
	vector<string> _DofsName, sDof;

	_LoadAmpName = loadblock._LoadAmpName;
	_DofsName = loadblock._DofsName;
	if (loadblock._Partunit.elesNo.size() > 0) {
		sDof = loadblock._Partunit.sDof;
	}

	iamp = GetAmpIndexFromName(_LoadAmpName);
	ComputeiAmp(iamp, t, ampvalue);

	if(_DofsName[0].find("P") != string::npos || sDof[0].find("P") != string::npos){
		ApplyPressure(mesh, dofHandler, calctype, ampvalue, loadblock, REXT);
	}
	else if (_DofsName[0].find("TR") != string::npos || sDof[0].find("TR") != string::npos) {
		ApplySurfTraction(mesh, dofHandler, calctype, ampvalue, loadblock, REXT);
	}
	//else if (_DofsName[0].find("TRVEC") != string::npos || sDof[0].find("TRVEC") != string::npos) {
	//	ApplyTRVEC(mesh, dofHandler, calctype, ampvalue, loadblock, REXT);
	//}
	else {
		MessagePrinter::PrintStars();
		MessagePrinter::PrintErrorTxt("  unexpected surface load type "_DofsName[0]","sDof[0]" !! ");
		MessagePrinter::AsFem_Exit();
	}
}

void LoadSystem::ApplyPressure(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &ampvalue, const LoadBlock &loadblock, Vec &REXT) {
	//_IsFollower
	bool _IsConstArea;
	string _LoadSetName, _eqSetName;

	vector<int> elesNo;
	vector<double> pvalues;
	//Vector3d _normals;
	double loadvalue;
	double earea, shellt;
	vector<pair<string, vector<int>>> _ElSetName2NodeIDsList;
	vector<int>    _Elmtids;
	int			   _nNodesPerElmt, _IthElmtJthNodeID, dofID;
	vector<double> _ElmtsVolume;
	vector<double> _LoadValues;
	Vector3d _pdirection;
	vector<pair<string, vector<int>>> _PhysicalName2ElmtIDsList;
	vector<pair<pair<string, string>, double>> _ShellESetsandMatTList;
	bool idseqset = true;

	_normals = 0.0;
	_LoadSetName = loadblock._LoadSetName;
	_LoadValues = loadblock._LoadValues;//��setֻ��һ�У���������ж���
	_IsConstArea = loadblock._IsConstArea;

	if (loadblock._Partunit.elesNo.size() > 0) {
		elesNo = loadblock._Partunit.elesNo;
		pvalues = loadblock._Partunit.values;
	}

	_PhysicalName2ElmtIDsList= mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr();
	if (_LoadSetName.size() > 0 && elesNo.size() == 0) {//eleset input
		_Elmtids = mesh.GetBulkMeshElmtIDsViaPhysicalName(_LoadSetName);
	}
	else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {//dicretize element input
		_Elmtids.clear();
		_Elmtids = elesNo;
		_LoadValues.clear();
		_LoadValues = pvalues;
	}
	int rankne, eStart, eEnd;
	int e, i, j, k, iInd;//iInd ��_IthNodeJthDofIndex, 
	vector<int> iconn;
	MPI_Comm_size(PETSC_COMM_WORLD, &_size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

	rankne = _Elmtids.size() / _size;
	eStart = _rank * rankne;
	eEnd = (_rank + 1)*rankne;
	if (_rank == _size - 1) eEnd = _Elmtids.size();
	_nNodesPerElmt = mesh.GetBulkMeshNodesNumPerBulkElmt;//һ�ֵ�Ԫ
	_ElmtsVolume = mesh.GetBulkMeshElmtVolumeListPtr();
	//���ڿ�/Ĥ��Ԫ ���ﲹ���ȡ��Ԫ���t �������������� [����������Ԫ��pressure��λ������������]
	_ShellESetsandMatTList= mesh.GetShellSectionESetsandMatTListPtr();
	if (_LoadSetName.size() > 0 && elesNo.size() == 0) {
		for (auto it : _ShellESetsandMatTList) {
			if (it[0] == _LoadSetName) {
				shellt = it[2];
				break;
			}
		}
	}
	else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {
		//ͨ����Ԫ ȷ��set��Ȼ��ͬ�ϻ�ȡ���
		for (auto it : _PhysicalName2ElmtIDsList) {
			if (_Elmtids.size() == it.second.size()) {
				for (int ie = 0; ie < _Elmtids.size();++ie) {
					if (_Elmtids[ie]!= it.second[ie]) {
						idseqset = false;
					}
					if ((ie == _Elmtids.size()-1)&&(idseqset==true)) {
						_eqSetName = it.first;
						break;
					}
				}
			}
		}
		if (_eqSetName.size() > 0) {
			for (auto it : _ShellESetsandMatTList) {
				if (it[0] == _eqSetName) {
					shellt = it[2];
					break;
				}
			}
		}
		else {
			MessagePrinter::PrintStars();
			MessagePrinter::PrintErrorTxt("  NOT find corresponding loadset for discritized elements Pressure load !! ");
			MessagePrinter::AsFem_Exit();
		}
	}

	for (e = eStart; e < eEnd; ++e) {
		evol = _ElmtsVolume[_Elmtids[e] - 1];
		if (_IsConstArea) {
			earea = evol / shellt;//��Ԫ����ǰԭʼ��� 
		}
		else {//===========��������������� ����=========================
			earea = evol / shellt;// xxx
			//_normals = 0.0;
			continue;
		}
		iconn = mesh.GetBulkMeshIthElmtNodesID(_Elmtids[e]);
		//_pdirection; ����ѹ����������ȷ��
		if (_LoadSetName.size() > 0 && elesNo.size() == 0) {
			for (i = 1; i <= _nNodesPerElmt; ++i) {
				j = iconn[i - 1];
				for (k = 1; k <= 3; ++k) {//xyz����  ���ɻ�ȡdim���
					loadvalue = earea *_LoadValues[0] / _nNodesPerElmt * _pdirection(k)* ampvalue;//����ѹ��/�ڵ���
					dofID = k;
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
		else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {
			for (i = 1; i <= _nNodesPerElmt; ++i) {
				j = iconn[i - 1];
				for (k = 1; k <= 3; ++k) {//xyz����  ���ɻ�ȡdim���
					loadvalue = earea * _LoadValues[e] / _nNodesPerElmt * _pdirection(k)* ampvalue;//����ѹ��/�ڵ���
					dofID = k;
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
	}
}

void LoadSystem::ApplySurfTraction(const Mesh &mesh, const DofHandler &dofHandler, FECalcType &calctype, const double &ampvalue, const LoadBlock &loadblock, Vec &REXT) {
	//for TRSHR always _IsFollower, but TRVEC can choose _IsFollower.
	bool _IsFollower, _IsConstArea;
	string _LoadSetName, _eqSetName;

	vector<int> elesNo;
	vector<double> tracvalues;
	//Vector3d _normals;
	double loadvalue;
	double earea, shellt;
	vector<pair<string, vector<int>>> _ElSetName2NodeIDsList;
	vector<int>    _Elmtids;
	int			   _nNodesPerElmt, _IthElmtJthNodeID, dofID;
	vector<double> _ElmtsVolume;
	vector<double> _LoadValues;
	Vector3d _tracdirection;
	vector<pair<string, vector<int>>> _PhysicalName2ElmtIDsList;
	vector<pair<pair<string, string>, double>> _ShellESetsandMatTList;
	bool idseqset = true;

	_normals = 0.0;
	_tracdirection = loadblock._direction;
	_IsConstArea = loadblock._IsConstArea;
	if (loadblock._Partunit.elesNo.size() > 0) {
		elesNo = loadblock._Partunit.elesNo;
		tracvalues = loadblock._Partunit.values;
	}
	_PhysicalName2ElmtIDsList = mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr();
	if (_LoadSetName.size() > 0 && elesNo.size() == 0) {//eleset input
		_LoadSetName = loadblock._LoadSetName;
		_LoadValues = loadblock._LoadValues;//��setֻ��һ�У���������ж���
		_Elmtids = mesh.GetBulkMeshElmtIDsViaPhysicalName(_LoadSetName);
	}
	else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {//dicretize element input
		_Elmtids.clear();
		_Elmtids = elesNo;
		_LoadValues.clear();
		_LoadValues = tracvalues;
	}

	int rankne, eStart, eEnd;
	int e, i, j, k, iInd;//iInd ��_IthNodeJthDofIndex, 
	vector<int> iconn;
	MPI_Comm_size(PETSC_COMM_WORLD, &_size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

	rankne = _Elmtids.size() / _size;
	eStart = _rank * rankne;
	eEnd = (_rank + 1)*rankne;
	if (_rank == _size - 1) eEnd = _Elmtids.size();
	_nNodesPerElmt = mesh.GetBulkMeshNodesNumPerBulkElmt;//һ�ֵ�Ԫ
	_ElmtsVolume = mesh.GetBulkMeshElmtVolumeListPtr();
	//���ڿ�/Ĥ��Ԫ ���ﲹ���ȡ��Ԫ���t �������������� [����������Ԫ��traction force��λ������������]
	_ShellESetsandMatTList = mesh.GetShellSectionESetsandMatTListPtr();
	if (_LoadSetName.size() > 0 && elesNo.size() == 0) {
		for (auto it : _ShellESetsandMatTList) {
			if (it[0] == _LoadSetName) {
				shellt = it[2];
				break;
			}
		}
	}
	else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {
		//ͨ����Ԫ ȷ��set��Ȼ��ͬ�ϻ�ȡ���
		for (auto it : _PhysicalName2ElmtIDsList) {
			if (_Elmtids.size() == it.second.size()) {
				for (int ie = 0; ie < _Elmtids.size(); ++ie) {
					if (_Elmtids[ie] != it.second[ie]) {
						idseqset = false;
					}
					if ((ie == _Elmtids.size() - 1) && (idseqset == true)) {
						_eqSetName = it.first;
						break;
					}
				}
			}
		}
		if (_eqSetName.size() > 0) {
			for (auto it : _ShellESetsandMatTList) {
				if (it[0] == _eqSetName) {
					shellt = it[2];
					break;
				}
			}
		}
		else {
			MessagePrinter::PrintStars();
			MessagePrinter::PrintErrorTxt("  NOT find corresponding loadset for discritized elements Traction load !! ");
			MessagePrinter::AsFem_Exit();
		}
	}

	for (e = eStart; e < eEnd; ++e) {
		evol = _ElmtsVolume[_Elmtids[e] - 1];
		if (_IsConstArea) {
			earea = evol / shellt;//��Ԫ����ǰԭʼ��� 
		}
		else {//===========��������������� ����=========================
			earea = evol / shellt;// xxx
			//_normals = 0.0;//_tracdirection
			continue;
		}
		iconn = mesh.GetBulkMeshIthElmtNodesID(_Elmtids[e]);
		//_tracdirection; ����ѹ����������ȷ��
		if (_LoadSetName.size() > 0 && elesNo.size() == 0) {
			for (i = 1; i <= _nNodesPerElmt; ++i) {
				j = iconn[i - 1];
				for (k = 1; k <= 3; ++k) {//xyz����  ���ɻ�ȡdim���
					loadvalue = earea * _LoadValues[0] / _nNodesPerElmt * _tracdirection(k)* ampvalue;//����ѹ��/�ڵ���
					if (_IsFollower) {
						//�ҳ��ڵ�ת���ǣ�normal
						//����ά��ά ��Ԫ�����й� �ο�applyBC
						//ͬ��


					}//follow�����ת�����κ��xyz�����ɶȷ���

					dofID = k;
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
		else if (_LoadSetName.size() == 0 && elesNo.size() > 0) {
			for (i = 1; i <= _nNodesPerElmt; ++i) {
				j = iconn[i - 1];
				for (k = 1; k <= 3; ++k) {//xyz����  ���ɻ�ȡdim���
					loadvalue = earea * _LoadValues[e] / _nNodesPerElmt * _tracdirection(k)* ampvalue;//����ѹ��/�ڵ���
					if (_IsFollower) {
						//�ҳ��ڵ�ת���ǣ�normal
						//����ά��ά ��Ԫ�����й� �ο�applyBC
						//ͬ��


					}//follow�����ת�����κ��xyz�����ɶȷ���
					dofID = k;
					iInd = dofHandler.GetIthNodeJthDofIndex(j, dofID) - 1;
					if (calctype == FECalcType::ComputeResidual) {
						VecSetValue(REXT, iInd, loadvalue, ADD_VALUES);
					}
				}
			}
		}
	}

}
