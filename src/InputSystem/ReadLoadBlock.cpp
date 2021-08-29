//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.13
//+++ Purpose: This function can read the [load] block and its 
//+++          subblock from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/StringUtils.h"
#include "InputSystem/InputSystem.h"


bool InputSystem::ReadLoadBlock(ifstream &in, string str, const int &lastendlinenum, int &linenum, Mesh &mesh, LoadSystem &loadSystem)
{
	// each load subblock should looks like:
	//   [load1]  i.e. name  or [loadINabaqus]
	//     type=cload [dload, dsload,user1,user2,user3...]
	//     setsname="" [nodeids=1 2 3 4 ...]
	//	   amplitude=""
	//     direction= 0 0 1
	//		...
	//   [end]
	// important: now , str already contains [loads] !!!

	bool HasLoadBlock = false;
	bool abaFlag = false;
	LoadBlock loadblock;
	string tempstr, str0, infostr, datastr;
	//vector<double> number;
	//vector<string> loadlist;
	string msg;
	string typestr;

	//bool HaskeyName = false;
	bool HasTypeName = false;
	bool HasSetNameorIDs = false;
	bool HasLoadData = false;

	loadblock.Init();
	vector<std::pair<string, std::vector<string>> > _LoadFromInp;
	int			   _Index=0;
	string         _LoadBlockName="none";
	
	bool		   _IsOPnew;
	bool           _IsFollower;
	bool           _IsConstArea;
	
	string         _LoadSetName;
	//struct PartUnit
	//{
	//	vector<int> nodesNo;//
	//	vector<int> elesNo;
	//	vector<string> sDof;//
	//	vector<double> values;

	//}_Loadunit;
	int inodeNo;//
	int ieleNo;
	//string iDof;//
	//double ivalue;

	string         _LoadkeyName;
	string         _LoadTypeName, _LoadAmpName;
	LoadType       _LoadType;
	vector<string> _DofsName;
	vector<double> _LoadValues;
	Vector3d _direction;

	vector<string> strvec;
	string substr;
	int ii;

	// now str="[loads]"
	while (linenum <= lastendlinenum) {
		getline(in, str); linenum += 1;
		str0 = str;
		str = StringUtils::StrToLower(str0);
		str = StringUtils::RemoveStrSpace(str);
		if (StringUtils::IsCommentLine(str) || str.size() < 1) continue;

		bool HasTypeName = false;
		bool HasSetNameorIDs = false;
		bool HasLoadData = false;

		if ((str.find('[') != string::npos&&str.find(']') != string::npos) &&
			(str.find("[end]") == string::npos&&str.find("[END]") == string::npos)) {
			HasLoadBlock = false;
			tempstr = StringUtils::RemoveStrSpace(str);//[]����loads    ���￪ʼ�ĳ�ABAQUS ���������ж� Ȼ���滻��������
			if (tempstr.size() < 3) {
				MessagePrinter::PrintErrorInLineNumber(linenum);
				MessagePrinter::PrintErrorTxt("no load block name can be found in [loads] sub block");
				MessagePrinter::AsFem_Exit();
				return false;
			}
			else if (tempstr.size() == 14 && tempstr.find("loadinabaqus") != string::npos) {//loadINabaqus
				abaFlag = true;
				_LoadFromInp = mesh.GetLoadInfoListPtr();
				//����label����Ϊtrue,�Ȼ�ȡ������Ϣ���ں���ѭ������addLoadblock
				HasLoadBlock = true; HasTypeName = true; HasSetNameorIDs = true; HasLoadData = true;
			}
			else {
				loadblock.Init();
				loadblock._LoadBlockName = tempstr.substr(1, tempstr.size() - 1 - 1);
				HasloadBlock = true;
			}
			while (str.find("[end]") == string::npos&&str.find("[END]") == string::npos) {
				if ((str.find("loadinabaqus") != string::npos)) { getline(in, str); linenum += 1; break; }//
				getline(in, str); linenum += 1;
				str0 = str;
				str = StringUtils::RemoveStrSpace(str0);
				str = StringUtils::StrToLower(str);
				if (StringUtils::IsCommentLine(str) || str.size() < 1) continue;

			}
			if (HasTypeName&&HasSetNameorIDs&&HasLoadData&&HasLoadBlock) {//���һ��Load���ݶ�ȡ����ӵ�loadSystem��
				if (abaFlag) {
					for (auto it : _LoadFromInp) {//��ȡ����load block����ֵ��Ȼ��add
						loadblock.Init();
						_Index++;
						loadblock._LoadIndex = _Index;
						loadblock._LoadBlockName = "load_"+to_string(_Index);
						infostr = it.first;
						datastr = it.second;
						strvec = StringUtils::SplitStr(infostr, ',');//keyword line
						if (infostr.find("*Cload") != string::npos) {
							loadblock._LoadkeyName = "Cload";
							for (auto sk : strvec) {
								if (sk.find("op=NEW") != string::npos) {
								//if (sk.find("op=") != string::npos) {
									//ii = sk.find_first_of('=');
									//substr = sk.substr(ii + 1, sk.length() - 1);
									//substr = StringUtils::RemoveStrSpace(substr);
									////if (substr.find("NEW") != string::npos) {
									////	loadblock._IsOPnew = "NEW";//MOD
									////}
									//loadblock._IsOPnew = substr;// bool not string
									loadblock._IsOPnew = true;
								}
								else if (sk.find("amplitude=") != string::npos) {
									ii = sk.find_first_of('=');
									substr = sk.substr(ii + 1, sk.length() - 1);
									substr = StringUtils::RemoveStrSpace(substr);
									loadblock._LoadAmpName = substr;
								}
								else if (sk.find("follower") != string::npos) {
									loadblock._IsFollower = true;
								}
							}

						}
						else if (infostr.find("*Dload") != string::npos) {
							loadblock._LoadkeyName = "Dload";
							for (auto sk : strvec) {
								if (sk.find("op=NEW") != string::npos) {
									loadblock._IsOPnew = true;
								}
								else if (sk.find("amplitude=") != string::npos) {
									ii = sk.find_first_of('=');
									substr = sk.substr(ii + 1, sk.length() - 1);
									substr = StringUtils::RemoveStrSpace(substr);
									loadblock._LoadAmpName = substr;
								}
								//else if (sk.find("follower") != string::npos) {// body force no follower
								//	loadblock._IsFollower = true;
								//}
							}
						}
						else if (infostr.find("*Dsload") != string::npos) {
							loadblock._LoadkeyName = "Dsload";
							loadblock._IsFollower = true;//������Ĭ����follower
							loadblock._IsConstArea = false;
							for (auto sk : strvec) {
								if (sk.find("op=NEW") != string::npos) {
									loadblock._IsOPnew = true;
								}
								else if (sk.find("amplitude=") != string::npos) {
									ii = sk.find_first_of('=');
									substr = sk.substr(ii + 1, sk.length() - 1);
									substr = StringUtils::RemoveStrSpace(substr);
									loadblock._LoadAmpName = substr;
								}
								else if (sk.find("follower=NO") != string::npos) {
									loadblock._IsFollower = false;
								}
								else if (sk.find("constant resultant=YES") != string::npos) {
									loadblock._IsConstArea = true;
								}
							}
						}
						else {
							MessagePrinter::PrintErrorInLineNumber(linenum);
							msg = "not supported Load keytype ["infostr"] in [" + loadblock._LoadIndex + "] sub block, please check your input !!";
							MessagePrinter::PrintErrorTxt(msg);
							MessagePrinter::AsFem_Exit();
						}

						_DofsName.clear();
						_LoadValues.clear();
						_direction.setZero();
						strvec.clear();
						for (auto di: datastr) {
							strvec = StringUtils::SplitStr(di, ',');//data line
							if (di.find("PART-1") != string::npos) {//PART-1-1. �������/��Ԫ��ֵ��������set
								ii = strvec[0].find_first_of('.');
								substr = di.substr(ii + 1, strvec[0].length() - 1);
								substr = StringUtils::RemoveStrSpace(substr);
								if (loadblock._LoadkeyName == "Cload") {//��� �ǵ�Ԫ
									inodeNo=atoi(substr.c_str());
									loadblock._Partunit.nodesNo.push_back(inodeNo);
								}
								else {
									ieleNo = atoi(substr.c_str());
									loadblock._Partunit.elesNo.push_back(ieleNo);
								}
								//substr = StringUtils::RemoveStrSpace(strvec[1]);
								//_DofsName.push_back(substr);
								////_LoadValues.push_back(atof(strvec[2].c_str()));
								//substr = StringUtils::RemoveStrSpace(strvec[2]);
								//_LoadValues.push_back(atof(substr.c_str()));
							}
							else {//set
								loadblock._LoadSetName= strvec[0];//ÿ��dataline��ͬ���ظ���ֵ
							}
							substr = StringUtils::RemoveStrSpace(strvec[1]);//������1-6�У��������ظ�
							_DofsName.push_back(substr);
							substr = StringUtils::RemoveStrSpace(strvec[2]);
							_LoadValues.push_back(atof(substr.c_str()));
							if (strvec.size > 3 && (_direction.norm()==0.0)) {
							//if (strvec.size > 3&&(strvec[1].find("GRAV") != string::npos|| strvec[1].find("TR") != string::npos)) {//strvec.size==6, have directions
								for (int k = 1; k < 4;k++) {
									substr = StringUtils::RemoveStrSpace(strvec[k+2]);
									_direction(k) = atof(substr.c_str());
								}
							}
						}

						loadblock._direction = _direction;
						if (loadblock._LoadSetName.size()>0) {//setģʽ
							loadblock._DofsName = _DofsName;
							loadblock._LoadValues= _LoadValues;
						}
						else {
							loadblock._Partunit.sDof = _DofsName;
							loadblock._Partunit.values = _LoadValues;
						}
						if (_DofsName[0]<='3') {
							loadblock._LoadTypeName = "CONCENTRATEDLC";
							loadblock._LoadType = LoadType::CONCENTRATEDLC;
						}
						else if (_DofsName[0] <= '6'&&_DofsName[0] > '3') {
							loadblock._LoadTypeName = "MOMENTLC";
							loadblock._LoadType = LoadType::MOMENTLC;
						}
						else if (_DofsName[0].find("TRSHR") != string::npos) {
							loadblock._LoadTypeName = "SURFACETRACTIONLC";
							loadblock._LoadType = LoadType::SURFACETRACTIONLC;
						}
						else if (_DofsName[0].find("TRVEC") != string::npos) {
							loadblock._LoadTypeName = "SURFACETRACTIONLC";
							loadblock._LoadType = LoadType::SURFACETRACTIONLC;
						}
						else if (_DofsName[0].find("GRAV") != string::npos) {
							loadblock._LoadTypeName = "GRAVITYLC";
							loadblock._LoadType = LoadType::GRAVITYLC;
						}
						else if (_DofsName[0].find("B") != string::npos) {
							loadblock._LoadTypeName = "BODYLC";
							loadblock._LoadType = LoadType::BODYLC;
						}

						loadSystem.AddLoadBlock2List(loadblock);
					}
				}
				else {// inner input; need to add
					HasLoadBlock = true;
					//......
					loadSystem.AddLoadBlock2List(loadblock);
				}
			}
			else {
				msg = "information is not complete in [loads] sub block, some information is missing in [" + loadblock._LoadBlockName + "]";
				MessagePrinter::PrintErrorTxt(msg, false);
				if (!HasTypeName) {
					msg = "no type found in [" + loadblock._LoadBlockName + "] sub block, 'dofname=' must be given in each [loads] sub block";
					MessagePrinter::PrintErrorTxt(msg);
				}
				if (!HasSetNameorIDs) {
					msg = "no SetNameorIDs found in [" + loadblock._LoadBlockName + "] sub block, 'setname=setname or IDs=ids' must be given in each [loads] sub block";
					MessagePrinter::PrintErrorTxt(msg);
				}
				if (!HasLoadData) {
					msg = "no load data found in [" + loadblock._LoadBlockName + "] sub block, 'loadvalue= 1.0 2.0 3.0 ...' must be given in each [loads] sub block";
					MessagePrinter::PrintErrorTxt(msg);
				}
				return false;
			}
		}

	}
	return true;
}
