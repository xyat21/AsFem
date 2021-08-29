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
//+++ Date   : 2021.02.18
//+++ Purpose: read mesh information from abaqus inp file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.07.21  add some funs for reading more abaqus inp info.

#include <numeric>
#include "Mesh/AbaqusIO.h"
#include "Utils/Vector3d.h"

bool AbaqusIO::ReadMeshFromFile(Mesh &mesh){
    string str,substr;
    vector<double> numbers;
    if(!_HasSetMeshFileName){
        MessagePrinter::PrintErrorTxt("can\'t read mesh, the mesh file name has not been set");
        return false;
    }
    _in.close();_in.clear();
    _in.open(_MeshFileName.c_str(),ios::in);
    if(!_in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t read the .inp file(="+_MeshFileName+"), please make sure the file name is correct or your inp file is in the same folder as your input file");
        MessagePrinter::AsFem_Exit();
    }
    //************************************************************
    //*** init all the mesh array and physical information
    //************************************************************
    mesh.GetBulkMeshNodeCoordsPtr().clear();
    mesh.GetBulkMeshElmtConnPtr().clear();
    mesh.GetBulkMeshElmtVolumeListPtr().clear();
	mesh.GetBulkMeshNodeeqVolumeListPtr().clear();

    mesh.GetBulkMeshElmtVTKCellTypeListPtr().clear();
    mesh.GetBulkMeshElmtPhyIDListPtr().clear();
    mesh.GetBulkMeshElmtDimListPtr().clear();
    mesh.GetBulkMeshElmtMeshTypeListPtr().clear();
    //*** for physical groups
    mesh.GetBulkMeshPhysicalGroupNameListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupIDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupDimListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupID2NameListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2IDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().clear();
    mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().clear();
    mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().clear();

	// =====XY======
	mesh.GetCLoadListPtr().clear();
	mesh.GetDsLoadListPtr().clear();
	mesh.GetDLoadsListPtr().clear();
	mesh.GetMateListPtr().clear();
	mesh.GetAmpsTimeListPtr().clear();
	mesh.GetBoundaryDoFListPtr().clear();
	mesh.GetBeamSectionESetsandMatAdListPtr().clear();
	mesh.GetSolidSectionESetsandMatAListPtr().clear();
	mesh.GetShellSectionESetsandMatTListPtr().clear();
	_BeamSectionNum=0;
	_SolidSectionNum=0;
	_ShellSectionNum=0;

	mesh.SetBeamSectionNum(GetBeamSectionNumFromInp());
	mesh.SetSolidSectionNum(GetSolidSectionNumFromInp());
	mesh.SetShellSectionNum(GetShellSectionNumFromInp());

	vector<std::pair<pair<string, string>, vector<double>> > t_CLoadList;
	vector<std::pair<vector<string>, double> > t_DsLoadList;
	vector<std::pair<vector<string>, vector<double>> > t_DLoadsList;
	vector<std::pair<string, std::vector<string>> > t_LoadStrList;
	vector<material_info> t_MateList;
	//vector<pair<string, pair<double, double>> > t_BeamSectionESetsandMatAdList;
	vector<pair<array<string, 3>, pair<array<double, 2>, array<double, 3>>> > t_BeamSectionESetsandMatAdList;
	vector<pair<pair<string, string>, double>> t_SolidSectionESetsandMatAList;
	vector<pair<pair<string, string>, double>> t_ShellSectionESetsandMatTList;
	//vector<std::pair<pair<string, string>, vector<std::pair<double, double> >> > t_AmpsTimeList;
	vector<std::pair<std::vector<std::string>, std::vector<double>> > t_AmpsTimeList;
	vector<std::pair<vector<string>, vector<double> > > t_AmpsTimeBlock;
	vector<std::pair<string, vector<int>> > t_BoundaryDoFList;
	bool endpartFlag = true;

    //*********************************
	mesh.SetSpaceDim(GetSpaceDimFromInp());
    mesh.SetBulkMeshDim(GetElmtDimFromInp());
    mesh.SetBulkMeshMinDim(GetElmtDimFromInp()-1);
    mesh.SetBulkMeshMeshOrder(GetElmtOrderFromInp());
    mesh.SetBulkMeshNodesNum(GetElmtNodesNumFromInp());
    mesh.SetBulkMeshElmtsNum(0);
    mesh.SetBulkMeshBulkElmtsNum(0);
    mesh.SetBulkMeshSurfaceElmtsNum(0);
    mesh.SetBulkMeshLineElmtsNum(0);
    mesh.SetBulkMeshPhysicalGroupNums(0);


    _nMaxDim=-1;_nMinDim=4; _nSpaceDim = -1;
    _nPhysicGroups=0;
    _nNodeSetPhysicalGroups=0;

	bool DimCheckFlag = true;
	_nSpaceDim = GetSpaceDimFromInp();
    _nMaxDim=GetElmtDimFromInp();
    _nMinDim=_nMaxDim-1;
	//如果是1D二维/三维单元，2D三维单元，则空间维度比单元维度大

    int nLowerDimElmts=0;


    _Xmax=-1.0e16;_Xmin=1.0e16;
    _Ymax=-1.0e16;_Ymin=1.0e16;
    _Zmax=-1.0e16;_Zmin=1.0e16;

    map<string,vector<int>> PhysicalName2NodeIDsSet;
	map<string, vector<int>> PhysicalName2ElmtIDsSet;
	vector<pair<string, vector<int>>> SurfaceElsets;

    vector<int> bulkelmtid;
	vector<double> elVol;

    PhysicalName2NodeIDsSet.clear();
	PhysicalName2ElmtIDsSet.clear();
	SurfaceElsets.clear();

    _nNodes=0;_nElmts=0;_nBulkElmts=0;
    while(!_in.eof()){
        // now we start to read inp file
        // remember: inp do not have physical group information explicitly!
        //           so we use a default value for all the bulk elements(physical id=0)
        // We! do! not! allow! multiple! parts! in your inp!!!
        getline(_in,str);
        if(str.find("*Node")!=string::npos&&str.find("Output") == string::npos){
            // start to read the nodes information
            int nodeid;
            double x,y,z;
            _nNodes=GetNodesNumFromInp();
            mesh.GetBulkMeshNodeCoordsPtr().resize(_nNodes*3,0.0);
            mesh.SetBulkMeshNodesNum(_nNodes);
			mesh.GetBulkMeshNodeeqVolumeListPtr().resize(_nNodes, 0.0);
            for(int i=0;i<_nNodes;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str,',');
                if((static_cast<int>(numbers.size())!=_nMaxDim+1)&& DimCheckFlag==true){
					DimCheckFlag = false;
                    MessagePrinter::PrintWarningTxt("incompatible node coordinates with your mesh dimension",false);
                    //MessagePrinter::AsFem_Exit();
                }
                nodeid=static_cast<int>(numbers[0]);
                x=numbers[1];y=numbers[2];z=0.0;
                if(_nSpaceDim ==3) z=numbers[3];//三维空间
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+1-1]=x;
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+2-1]=y;
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+3-1]=z;

                if(x>_Xmax) _Xmax=x;
                if(x<_Xmin) _Xmin=x;
                if(y>_Ymax) _Ymax=y;
                if(y<_Ymin) _Ymin=y;
                if(z>_Zmax) _Zmax=z;
                if(z<_Zmin) _Zmin=z;
            }
			mesh.SetBulkMeshXmax(_Xmax);
			mesh.SetBulkMeshXmin(_Xmin);
			mesh.SetBulkMeshYmax(_Ymax);
			mesh.SetBulkMeshYmin(_Ymin);
			mesh.SetBulkMeshZmax(_Zmax);
			mesh.SetBulkMeshZmin(_Zmin);
		} // end-of-node-block
        else if(str.find("*Element")!=string::npos&&str.find("Output") == string::npos){
            _nBulkElmts=GetElmtsNumFromInp();// this is only the bulk elements number  // 这是所有单元啊···包括了边界的
            _nNodesPerBulkElmt=GetElmtNodesNumFromInp();
            nLowerDimElmts=GetSurfaceElmtsNumFromInp();//仅针对三维单元or 包括shell单元时施加的压力单元？
            nLowerDimElmts=0;// currently, we do not account for surface element!!!=================
            _nElmts=_nBulkElmts+nLowerDimElmts;			//理解下这里 nLowerDimElmts=SurfaceElmtsNum=0
            mesh.SetBulkMeshElmtsNum(_nElmts);
            mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
            mesh.SetBulkMeshNodesNumPerBulkElmt(_nNodesPerBulkElmt);
            bulkelmtid.clear(); 
            if(_nMaxDim==2){
                mesh.SetBulkMeshLineElmtsNum(nLowerDimElmts);
				//mesh.SetBulkMeshLineElmtsNum(0);
                mesh.SetBulkMeshNodesNumPerLineElmt(GetSubElmtNodesNumFromInp());//二维平面单元的子单元是线单元，但注意有二维坐标的线单元的子单元得到是点[所以这句不能用于beam类单元]
                mesh.SetBulkMeshLineMeshType(GetSubElmtMeshTypeFromInp());
            }else if(_nMaxDim==3){
                mesh.SetBulkMeshSurfaceElmtsNum(nLowerDimElmts);
                mesh.SetBulkMeshNodesNumPerSurfaceElmt(GetSubElmtNodesNumFromInp());
                mesh.SetBulkMeshNodesNumPerLineElmt(GetSubSubElmtNodesNumFromInp());

                mesh.SetBulkMeshSurfaceMeshType(GetSubElmtMeshTypeFromInp());
                mesh.SetBulkMeshLineMeshType(GetSubSubElmtMeshTypeFromInp());
            }
            mesh.SetBulkMeshMeshType(GetElmtMeshTypeFromInp());
            mesh.SetBulkMeshMeshTypeName(GetElmtMeshTypeNameFromInp());
            // allocate memory for connectivity and other array
            mesh.GetBulkMeshElmtConnPtr().resize(_nElmts,vector<int>(0));
            mesh.GetBulkMeshElmtVTKCellTypeListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtPhyIDListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtMeshTypeListPtr().resize(_nElmts,MeshType::NULLTYPE);
            mesh.GetBulkMeshElmtDimListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtVolumeListPtr().resize(_nElmts,0.0);

            int phyid,vtktype,elmtid;
            MeshType meshtype;
            phyid=0;// we use 0 for all the bulk mesh !!!
            vtktype=GetElmtVTKCellTypeFromInp();
            meshtype=GetElmtMeshTypeFromInp();
            for(int e=0;e<_nBulkElmts;e++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(static_cast<int>(numbers.size())!=1+_nNodesPerBulkElmt){
                    MessagePrinter::PrintErrorTxt("error detected, your nodes number per bulk element is not match with your mesh type, please check your inp file");
                    MessagePrinter::AsFem_Exit();
                }
                elmtid=static_cast<int>(numbers[0]);
                bulkelmtid.push_back(elmtid);

                mesh.GetBulkMeshElmtConnPtr()[elmtid-1].resize(_nNodesPerBulkElmt+1,0);
                mesh.GetBulkMeshElmtConnPtr()[elmtid-1][0]=_nNodesPerBulkElmt;
                for(int j=1;j<=_nNodesPerBulkElmt;j++){
                    mesh.GetBulkMeshElmtConnPtr()[elmtid-1][j]=static_cast<int>(numbers[j]);
                }

                // set bulk elements' other properties
                mesh.GetBulkMeshElmtMeshTypeListPtr()[elmtid-1]=meshtype;
                mesh.GetBulkMeshElmtPhyIDListPtr()[elmtid-1]=phyid;
                mesh.GetBulkMeshElmtVTKCellTypeListPtr()[elmtid-1]=vtktype;
                mesh.GetBulkMeshElmtDimListPtr()[elmtid-1]=_nMaxDim;
                mesh.GetBulkMeshElmtVolumeListPtr()[elmtid-1]=0.0;//后续待计算
            }
        }// end-of-element-block
        else if(str.find("*Surface,")!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            string SurfaceSetName;
            nLowerDimElmts=GetSurfaceElmtsNumFromInp();
            int surfaceid,i,edge;
            vector<int> ElmtIDList;
            i=str.find_last_of("=");
            SurfaceSetName=str.substr(i+1,string::npos);
            SurfaceSetName.erase(remove(SurfaceSetName.begin(), SurfaceSetName.end(), '\t'), SurfaceSetName.end());
            ElmtIDList=GetSurfaceElmtIDViaSurfaceNameFromInp(SurfaceSetName);
			SurfaceElsets.push_back(make_pair(SurfaceSetName, ElmtIDList));
            edge=GetSurfaceEdgeIDViaSurfaceNameFromInp(SurfaceSetName);

            if(edge||ElmtIDList.size()||surfaceid){}
            // TODO: here we need some strategy to split the boundary element from the surface set
        }
        else if(str.find("*Nset,")!=string::npos){
            int i,j;
            string nodesetname;
			if (str.find("\"") != string::npos) { str = StringUtils::RemoveSymbolFromStr(str, '\"'); }//*Nset, nset="Edge-Seeds-2"
            i=str.find("=");j=str.find_last_of(",");//区分在part里的set是为赋材料属性，且与后面Nset命名有重复，这里给其名字后添加“mat”，相应在读取Nset时也要将命名区分part和assembly
			if (str.find("instance") != string::npos&&str.find("generate") != string::npos) {
				str = str.substr(0, j - 1);
				j = str.find_last_of(",");
			}//instance nset用generate方式时 有3个逗号
			nodesetname = str.substr(i + 1, j - i);//只有一个逗号时，j<i，无效则会取到字符串尾
			if (nodesetname.find(",") != string::npos) { nodesetname.pop_back(); }
			if (endpartFlag) { nodesetname = nodesetname + "part"; }//
			           
            nodesetname.erase(remove(nodesetname.begin(), nodesetname.end(), '\t'), nodesetname.end());
            nodesetname.erase(remove(nodesetname.begin(), nodesetname.end(), '\n'), nodesetname.end());
            PhysicalName2NodeIDsSet[nodesetname].clear();
            if(str.find("generate")!=string::npos){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(static_cast<int>(numbers.size())!=3){
                    MessagePrinter::PrintErrorTxt("invalid node id index information, you need min-id,max-id,increament for generate node set. Please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
                for(i=static_cast<int>(numbers[0]);i<=static_cast<int>(numbers[1]);i+=static_cast<int>(numbers[2])){
                    PhysicalName2NodeIDsSet[nodesetname].push_back(i);
                }
            }
            else{
                for(i=0;i<numeric_limits<int>::max()-1;i++){
                    getline(_in,str);
					if (str.find("*") != string::npos) { 
						int iPos = _in.tellg();//得到当前位置
						_in.seekg(iPos - str.length()-2);//将文件指针移到上一行
						break; 
					}
                    numbers=StringUtils::SplitStrNum(str);
                    for(const auto &val:numbers){
                        PhysicalName2NodeIDsSet[nodesetname].push_back(static_cast<int>(val));
                    }
                }
            }
        }
		else if (str.find("*Elset,") != string::npos) {
			int i, j;
			vector<string> strvec;
			string elsetname;
			if(str.find("\"") != string::npos){str= StringUtils::RemoveSymbolFromStr(str, '\"');}//elset="Edge Seeds-1", generate
			//i = str.find("="); j = str.find_last_of(",");
			//if (str.find("internal") != string::npos) { 
			//	str= str.substr(0, j-1);
			//	j = str.find_last_of(",");
			//}//内部压力分组 有4个逗号
			//elsetname = str.substr(i + 1, j - i);//有的没有第二个逗号，会自动取到字符串尾;而有第二个逗号时索引会大1，因为最后的空格

			strvec.clear();
			strvec = StringUtils::SplitStr(str, ',');
			i = strvec[1].find("=");
			elsetname = strvec[1].substr(i + 1);

			if (elsetname.find(",") != string::npos) { elsetname.pop_back(); }
			if (endpartFlag) { elsetname = elsetname + "part"; }//区分在part里的set是为赋材料属性，且与后面Nset命名有重复，这里给其名字后添加“part”，相应在读取Nset时也要将命名区分part和assembly

			elsetname.erase(remove(elsetname.begin(), elsetname.end(), '\t'), elsetname.end());
			elsetname.erase(remove(elsetname.begin(), elsetname.end(), '\n'), elsetname.end());
			PhysicalName2ElmtIDsSet[elsetname].clear();
			if (str.find("generate") != string::npos) {
				getline(_in, str);
				numbers = StringUtils::SplitStrNum(str);
				if (static_cast<int>(numbers.size()) != 3) {
					MessagePrinter::PrintErrorTxt("invalid elment id index information, you need min-id,max-id,increament for generate elment set. Please check your input file");
					MessagePrinter::AsFem_Exit();
				}
				for (i = static_cast<int>(numbers[0]); i <= static_cast<int>(numbers[1]); i += static_cast<int>(numbers[2])) {
					PhysicalName2ElmtIDsSet[elsetname].push_back(i);
				}
			}
			else {
				for (i = 0; i < numeric_limits<int>::max() - 1; i++) {
					getline(_in, str);
					//if (str.find("*") != string::npos) break;
					if (str.find("*") != string::npos) {
						int iPos = _in.tellg();//得到当前位置
						_in.seekg(iPos - str.length()-2);//将文件指针移到上一行
						break;
					}
					numbers = StringUtils::SplitStrNum(str);
					for (const auto &val : numbers) {
						PhysicalName2ElmtIDsSet[elsetname].push_back(static_cast<int>(val));
					}
				}
			}
		}
		else if (str.find("*End Part") != string::npos) {
			endpartFlag = false;
		}
		else if (str.find("*Beam Section,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Solid Section,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Shell Section,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Amplitude,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Material,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Boundary") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Cload,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Dsload,") != string::npos) {
			//int i, j;


		}
		else if (str.find("*Dload,") != string::npos) {
			//int i, j; 


		}
	}// end-of-file-reading

	//mesh._BeamSectionESetsandMatAdList;
	//mesh._SolidSectionESetsandMatAList;
	//mesh._ShellSectionESetsandMatTList;
	t_LoadStrList = GetLoadStrListFromInp();
	t_CLoadList= GetCLoadFromInp();
	t_DsLoadList= GetDsLoadFromInp();
	t_DLoadsList= GetDloadsFromInp();
	t_MateList= GetMateFromInp();
	t_AmpsTimeList = GetAmpsTimeFromInp();

	t_BoundaryDoFList = GetBoundaryDofFromInp();
	t_BeamSectionESetsandMatAdList= GetBeamSectionESetsandMatAdFromInp();
	t_SolidSectionESetsandMatAList= GetSolidSectionESetsandMatAFromInp();
	t_ShellSectionESetsandMatTList= GetShellSectionESetsandMatTFromInp();
	//mesh.GetCLoadListPtr().resize(t_CLoadList.size());
	//mesh.GetDsLoadListPtr().resize(t_DsLoadList.size());
	//mesh.GetDLoadsListPtr().resize(t_DLoadsList.size());
	mesh.GetMateListPtr().resize(t_MateList.size());
	//mesh.GetAmpsTimeListPtr().resize(t_AmpsTimeList.size());
	//mesh.GetBoundaryDoFListPtr().resize(t_BoundaryDoFList.size());
	//mesh.GetBeamSectionESetsandMatAdListPtr().resize(t_BeamSectionESetsandMatAdList.size());
	//mesh.GetSolidSectionESetsandMatAListPtr().resize(t_SolidSectionESetsandMatAList.size());
	//mesh.GetShellSectionESetsandMatTListPtr().resize(t_ShellSectionESetsandMatTList.size());
	for (auto it : t_CLoadList) { mesh.GetCLoadListPtr().push_back(it);}
	for (auto it : t_DsLoadList) { mesh.GetDsLoadListPtr().push_back(it); }
	for (auto it : t_DLoadsList) { mesh.GetDLoadsListPtr().push_back(it); }

	for (auto it : t_LoadStrList) { mesh.GetLoadInfoListPtr().push_back(it); }

	//for (auto it : t_MateList) { mesh.GetMateListPtr().push_back(it); }
	int si = 0;
	for (auto it : t_MateList) { 
		si++;
		//mesh.GetMateListPtr()[si-1]->damp_alpha.push_back(it.damp_alpha);//表达式必须包含类类型
		mesh.GetMateListPtr()[si - 1].matname = it.matname;
		mesh.GetMateListPtr()[si - 1].damp_alpha=it.damp_alpha;
		mesh.GetMateListPtr()[si - 1].damp_beta = it.damp_beta;
		mesh.GetMateListPtr()[si - 1].density = it.density;
		mesh.GetMateListPtr()[si - 1].poisson = it.poisson;
		mesh.GetMateListPtr()[si - 1].Youngs_modulus = it.Youngs_modulus;
	}
	for (auto it : t_AmpsTimeList) { mesh.GetAmpsTimeListPtr().push_back(it); }
	for (auto it : t_BoundaryDoFList) { mesh.GetBoundaryDoFListPtr().push_back(it); }

	for (auto it : t_BeamSectionESetsandMatAdList) { mesh.GetBeamSectionESetsandMatAdListPtr().push_back(it); }
	for (auto it : t_SolidSectionESetsandMatAList) { mesh.GetSolidSectionESetsandMatAListPtr().push_back(it); }
	for (auto it : t_ShellSectionESetsandMatTList) { mesh.GetShellSectionESetsandMatTListPtr().push_back(it); }

    _nPhysicGroups=0;//单元组名字从0开始
    mesh.GetBulkMeshPhysicalGroupNameListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupIDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupDimListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2DimListPrt().clear();
    mesh.GetBulkMeshPhysicalGroupName2IDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupIDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().clear();
    mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().clear();

    //********************************************************************
    mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back("alldomain");
    mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(_nPhysicGroups);
    mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(_nMaxDim);
    mesh.GetBulkMeshPhysicalGroupName2DimListPrt().push_back(make_pair("alldomain",_nMaxDim));
    mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",_nPhysicGroups));
	mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(_nPhysicGroups, "alldomain"));
    mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().push_back(make_pair("alldomain",_nNodesPerBulkElmt));
    mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkelmtid));//第0组为所有单元组
	//_nPhysicGroups += 1;// for our bulk elements


	if(SurfaceElsets.size()>0){
		string surelsetname;
		vector<int> surIDs;
		for (auto it : SurfaceElsets) {
			surelsetname=it.first;
			surIDs= it.second;
			PhysicalName2ElmtIDsSet[surelsetname].clear();
			for (auto sid : surIDs) { PhysicalName2ElmtIDsSet[surelsetname].push_back(sid); }
		}
	}
	string phyElname;
	for (auto it : PhysicalName2ElmtIDsSet) {
		phyElname = it.first;
		phyElname.erase(remove(phyElname.begin(), phyElname.end(), '\t'), phyElname.end());
		phyElname.erase(remove(phyElname.begin(), phyElname.end(), '\n'), phyElname.end());
		//phyname.pop_back();// inp file has the last unknown character	会导致第一个set-1的“1”被删掉。。。 不用后续组又会多个,
		if(phyElname.find(",") != string::npos){ phyElname.pop_back(); }

        _nPhysicGroups+=1;
        mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back(phyElname);
        mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(_nPhysicGroups);
        mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(_nMaxDim);
        mesh.GetBulkMeshPhysicalGroupName2DimListPrt().push_back(make_pair(phyElname,_nMaxDim));
        mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair(phyElname,_nPhysicGroups));
		mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(_nPhysicGroups, phyElname));
        mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().push_back(make_pair(phyElname, _nNodesPerBulkElmt));
        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyElname,it.second));
	}

    string phyname;
    _nNodeSetPhysicalGroups=0;
    for(auto it:PhysicalName2NodeIDsSet){
        phyname=it.first;
        phyname.erase(remove(phyname.begin(), phyname.end(), '\t'),phyname.end());
        phyname.erase(remove(phyname.begin(), phyname.end(), '\n'),phyname.end());
        //phyname.pop_back();// inp file has the last unknown character	会导致第一个set-1的“1”被删掉。。。
		if (phyname.find(",") != string::npos) { phyname.pop_back(); }

//        _nPhysicGroups+=1;
//        mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back(phyname);
//        mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(_nPhysicGroups);
//        mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(0);
//        mesh.GetBulkMeshPhysicalGroupName2DimListPrt().push_back(make_pair(phyname,0));
//        mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair(phyname,_nPhysicGroups));
//        mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().push_back(make_pair(phyname,1));
//        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyname,it.second));

        _nNodeSetPhysicalGroups+=1;//节点组名字从1开始
        mesh.GetBulkMeshNodeSetPhysicalNameListPtr().push_back(phyname);
        mesh.GetBulkMeshNodeSetPhysicalIDListPtr().push_back(_nNodeSetPhysicalGroups);
        mesh.GetBulkMeshNodeSetPhysicalID2NameListPtr().push_back(make_pair(_nNodeSetPhysicalGroups,phyname));
        mesh.GetBulkMeshNodeSetPhysicalName2IDListPtr().push_back(make_pair(phyname,_nNodeSetPhysicalGroups));
        mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().push_back(make_pair(phyname,it.second));

    }
    mesh.SetBulkMeshPhysicalGroupNums(_nPhysicGroups+1);//第一组编号为0，所以总组数+1
    mesh.SetBulkMeshNodeSetPhysicalGroupNums(_nNodeSetPhysicalGroups);//第一组编号为1

	//calculate the element volume, and the accumulated node volume of (elment V averaged to connected nodes)
	//default: the whole model has only one meshtype.
	MeshType mmeshtype;
	mmeshtype=mesh.GetBulkMeshBulkElmtType();
	if(mmeshtype == MeshType::EDGE2 || mmeshtype == MeshType::BEAM2 || mmeshtype == MeshType::TRUSS2 || mmeshtype == MeshType::PIPE2){
		//单元号->节点号->节点坐标->单元长度->节点分配到的长度自加（因为单元会共用节点）	//MeshOrder高阶单元暂空
		int i, eID, nIDa, nIDb;
		vector<Vector3d> elNcoords;
		//vector<double> ellen;
		elVol.clear();
		elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
		elVol.resize(_nBulkElmts, 0.0);
		for (i = 1; i <= _nBulkElmts; ++i) {//_nBulkElmts=mesh.GetBulkMeshBulkElmtsNum()
			eID=bulkelmtid[i-1];//单元号
			nIDa= mesh.GetBulkMeshElmtConnPtr()[eID - 1][1];//第1节点号
			nIDb= mesh.GetBulkMeshElmtConnPtr()[eID - 1][2];//第2节点号
			elNcoords[0](1) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDa - 1) * 3 + 1 - 1];//edge两端坐标xyz
			elNcoords[0](2) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDa - 1) * 3 + 2 - 1];
			elNcoords[0](3) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDa - 1) * 3 + 3 - 1];
			elNcoords[1](1) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDb - 1) * 3 + 1 - 1];
			elNcoords[1](2) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDb - 1) * 3 + 2 - 1];
			elNcoords[1](3) = mesh.GetBulkMeshNodeCoordsPtr()[(nIDb - 1) * 3 + 3 - 1];
			//elVol[eID - 1]  =  sqrt((elNcoords[1](1) - elNcoords[0](1))*(elNcoords[1](1) - elNcoords[0](1)) +
			//						(elNcoords[1](2) - elNcoords[0](2))*(elNcoords[1](2) - elNcoords[0](2)) +
			//						(elNcoords[1](3) - elNcoords[0](3))*(elNcoords[1](3) - elNcoords[0](3)));
			elVol[eID - 1] =(elNcoords[1]-elNcoords[0]).norm();
			mesh.GetBulkMeshElmtVolumeListPtr()[eID - 1] = elVol[eID - 1];
			mesh.GetBulkMeshNodeeqVolumeListPtr()[nIDa - 1] += 0.5*elVol[eID - 1];
			mesh.GetBulkMeshNodeeqVolumeListPtr()[nIDb - 1] += 0.5*elVol[eID - 1];
		}

	}
	else if (mmeshtype == MeshType::TRI3) {
		int i, j, k, eID;
		vector<int>	nID;
		vector<Vector3d> elNcoords;
		elVol.clear();
		nID.resize(_nNodesPerBulkElmt, 0);
		elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
		elVol.resize(_nBulkElmts, 0.0);
		for (i = 1; i <= _nBulkElmts; i++) {//_nBulkElmts=mesh.GetBulkMeshBulkElmtsNum()
			eID = bulkelmtid[i - 1];//单元号
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {//tri3 j=1,2,3
				nID[j-1]=mesh.GetBulkMeshElmtConnPtr()[eID - 1][j];//节点号 注意j=0是单元节点数，j=n依次是节点号
				for (k = 1; k <= 3; k++) {
					elNcoords[j-1](k) = mesh.GetBulkMeshNodeCoordsPtr()[(nID[j - 1] - 1) * 3 + k - 1];//elmt各点坐标xyz
				}
			}
			elVol[eID - 1] = 0.5*((elNcoords[1]-elNcoords[0])^(elNcoords[2]-elNcoords[0])).norm();//S=|a×b|/2=0.5*|(Y1Z2-Y2Z1,Z1X2-Z2X1,X1Y2-X2Y1)|
			mesh.GetBulkMeshElmtVolumeListPtr()[eID - 1] = elVol[eID - 1];
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {
				mesh.GetBulkMeshNodeeqVolumeListPtr()[nID[j - 1] - 1] += elVol[eID - 1]/_nNodesPerBulkElmt;
			}
		}
	}
	else if (mmeshtype == MeshType::TRI6) {}
	else if (mmeshtype == MeshType::QUAD4) {
		int i, j, k, eID;
		vector<int>	nID;
		vector<Vector3d> elNcoords;
		elVol.clear();
		nID.resize(_nNodesPerBulkElmt, 0);
		elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
		elVol.resize(_nBulkElmts, 0.0);
		for (i = 1; i <= _nBulkElmts; i++) {//
			eID = bulkelmtid[i - 1];//单元号
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {//quad4 j=1,2,3,4
				nID[j - 1] = mesh.GetBulkMeshElmtConnPtr()[eID - 1][j];
				for (k = 1; k <= 3; k++) {
					elNcoords[j - 1](k) = mesh.GetBulkMeshNodeCoordsPtr()[(nID[j - 1] - 1) * 3 + k - 1];//elmt各点坐标xyz
				}
			}
			elVol[eID - 1] = 0.5*((elNcoords[1] - elNcoords[0]) ^ (elNcoords[2] - elNcoords[0])).norm()+
							 0.5*((elNcoords[2] - elNcoords[0]) ^ (elNcoords[3] - elNcoords[0])).norm();//S=Sabc+Sacd
			mesh.GetBulkMeshElmtVolumeListPtr()[eID - 1] = elVol[eID - 1];
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {
				mesh.GetBulkMeshNodeeqVolumeListPtr()[nID[j - 1] - 1] += elVol[eID - 1] / _nNodesPerBulkElmt;
			}
		}
	}
	else if (mmeshtype == MeshType::QUAD8) {}
	else if (mmeshtype == MeshType::TET4) {
		int i, j, k, eID;
		vector<int>	nID;
		vector<Vector3d> elNcoords;
		elVol.clear();
		nID.resize(_nNodesPerBulkElmt, 0);
		elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
		elVol.resize(_nBulkElmts, 0.0);
		for (i = 1; i <= _nBulkElmts; i++) {//
			eID = bulkelmtid[i - 1];//单元号
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {//tet4 j=1,2,3,4
				nID[j - 1] = mesh.GetBulkMeshElmtConnPtr()[eID - 1][j];
				for (k = 1; k <= 3; k++) {
					elNcoords[j - 1](k) = mesh.GetBulkMeshNodeCoordsPtr()[(nID[j - 1] - 1) * 3 + k - 1];//elmt各点坐标xyz
				}
			}
			elVol[eID - 1] = 1.0/6.0*abs(((elNcoords[1] - elNcoords[0]) ^ (elNcoords[2] - elNcoords[0]))*(elNcoords[3] - elNcoords[0])) ;//S=|(a×b)·c|/6
			mesh.GetBulkMeshElmtVolumeListPtr()[eID - 1] = elVol[eID - 1];
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {
				mesh.GetBulkMeshNodeeqVolumeListPtr()[nID[j - 1] - 1] += elVol[eID - 1] / _nNodesPerBulkElmt;
			}
		}
	}
	else if (mmeshtype == MeshType::HEX8) {
		int i, j, k, eID;
		vector<int>	nID;
		vector<Vector3d> elNcoords;
		elVol.clear();
		nID.resize(_nNodesPerBulkElmt, 0);
		elNcoords.resize(_nNodesPerBulkElmt, Vector3d());
		elVol.resize(_nBulkElmts, 0.0);
		for (i = 1; i <= _nBulkElmts; i++) {//
			eID = bulkelmtid[i - 1];//单元号
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {//hex8 j=1,2,3,4  5,6,7,8
				nID[j - 1] = mesh.GetBulkMeshElmtConnPtr()[eID - 1][j];
				for (k = 1; k <= 3; k++) {
					elNcoords[j - 1](k) = mesh.GetBulkMeshNodeCoordsPtr()[(nID[j - 1] - 1) * 3 + k - 1];//elmt各点坐标xyz
				}
			}
			//取内部一中点P，分解为6个棱锥=12个四面体；P123/P134	P567/P578	P156/P162	P487/P473	P148/P185	P237/P276
			Vector3d P(0.0);//取为1与7中点
			P = (elNcoords[1 - 1]+ elNcoords[7 - 1])*0.5;
			elVol[eID - 1] = 1.0 / 6.0*(abs(((elNcoords[0] - P) ^ (elNcoords[1] - P))*(elNcoords[2] - P)) +
										abs(((elNcoords[0] - P) ^ (elNcoords[2] - P))*(elNcoords[3] - P)) +
										abs(((elNcoords[4] - P) ^ (elNcoords[5] - P))*(elNcoords[6] - P)) +
										abs(((elNcoords[4] - P) ^ (elNcoords[6] - P))*(elNcoords[7] - P)) +
										abs(((elNcoords[0] - P) ^ (elNcoords[4] - P))*(elNcoords[5] - P)) + 
										abs(((elNcoords[0] - P) ^ (elNcoords[5] - P))*(elNcoords[1] - P)) +

										abs(((elNcoords[3] - P) ^ (elNcoords[7] - P))*(elNcoords[6] - P)) + 
										abs(((elNcoords[3] - P) ^ (elNcoords[6] - P))*(elNcoords[2] - P)) + 
										abs(((elNcoords[0] - P) ^ (elNcoords[3] - P))*(elNcoords[7] - P)) + 
										abs(((elNcoords[0] - P) ^ (elNcoords[7] - P))*(elNcoords[4] - P)) + 
										abs(((elNcoords[1] - P) ^ (elNcoords[2] - P))*(elNcoords[6] - P)) + 
										abs(((elNcoords[1] - P) ^ (elNcoords[6] - P))*(elNcoords[5] - P)) );//S=sigma_(1-12)[|(a×b)·c|/6]

			mesh.GetBulkMeshElmtVolumeListPtr()[eID - 1] = elVol[eID - 1];//计算体积 验算有微小的差异
			for (j = 1; j <= _nNodesPerBulkElmt; j++) {
				mesh.GetBulkMeshNodeeqVolumeListPtr()[nID[j - 1] - 1] += elVol[eID - 1] / _nNodesPerBulkElmt;
			}
		}
	}
	else if (mmeshtype == MeshType::TET10) {}
	else if (mmeshtype == MeshType::HEX20) {}
	else if (mmeshtype == MeshType::QUAD9) {}
	else {
		MessagePrinter::PrintWarningTxt("Failed to calculate the elment volume. Unsupported element type for Volume Calculation");
	}
	double modelVolume;
	modelVolume=accumulate(elVol.begin(), elVol.end(), 0.0);
	mesh.SetTotalVolume(modelVolume);//====实际上在FormBulkFE 277行有通过形函数计算ElmtVolume _BulkVolumes。。。========================


    return true;
}