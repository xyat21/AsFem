//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Xiaoyuan
//+++ Date   : 2021.07.13
//+++ Purpose: This function can read the [amplitude] block and its 
//+++          subblock from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"


bool InputSystem::ReadAmpBlock(ifstream &in, string str, const int &lastendlinenum, int &linenum, Mesh &mesh, AmpSystem &ampSystem){
	// each amplitude subblock should looks like:
	//   [amplitude1]  i.e. name  or [ampINabaqus]
	//     type=tabular [periodic,user1,user2,user3...]
	//     smooth=0.0 [default is 0.0, so it is not necessary to be given!!!]
	//	   //fixinterval=0.0
	//     params=0 0  0.1  0.5  0.3  0.4  0.5  0.8 ...
	//   [end]
	// important: now , str already contains [amps] !!!

	bool HasAmpBlock = false;
	bool abaFlag = false;
	AmpBlock ampblock;
	string tempstr, str0;
	//vector<double> number;
	//vector<string> amplist;
	string msg;
	string typestr;

	bool HasTypeName = false;
	bool HasParams = false;
	bool HasSmooth = false;

	ampblock.Init();
	vector<std::pair<std::vector<std::string>, std::vector<double>> > _AmpFromInp;
	vector<double> _AmpParams;
	int			   _Index;
	double		   _Smoothdata;
	double		   _FixInterval;
	string		   _AmpBlockName;
	string		   _AmpTypeName;

	// now str="[amps]"
	while (linenum <= lastendlinenum) {
		getline(in, str); linenum += 1;
		str0 = str;
		str = StringUtils::StrToLower(str0);
		str = StringUtils::RemoveStrSpace(str);
		if (StringUtils::IsCommentLine(str) || str.size() < 1) continue;

		bool HasTypeName = false;
		bool HasParams = false;
		bool HasSmooth = false;

		if ((str.find('[') != string::npos&&str.find(']') != string::npos) &&
			(str.find("[end]") == string::npos&&str.find("[END]") == string::npos)) {
			HasAmpBlock = false;
			tempstr = StringUtils::RemoveStrSpace(str);//[]内是amps    这里开始改成ABAQUS 下面增加判断 然后替换读入数据
			if (tempstr.size() < 3) {
				MessagePrinter::PrintErrorInLineNumber(linenum);
				MessagePrinter::PrintErrorTxt("no amplitude block name can be found in [amps] sub block");
				MessagePrinter::AsFem_Exit();
				return false;
			}
			else if (tempstr.size() == 13 && tempstr.find("ampinabaqus") != string::npos) {//ampINabaqus
				abaFlag = true;
				_AmpFromInp = mesh.GetAmpsTimeListPtr();
				//这里label都设为true,先获取基本信息，在后面循环建立addAmpblock
				HasAmpBlock = true; HasTypeName = true; HasParams = true; 
			}
			else {
				ampblock.Init();
				ampblock._AmpBlockName = tempstr.substr(1, tempstr.size() - 1 - 1);
				HasAmpBlock = true;
			}
			while (str.find("[end]") == string::npos&&str.find("[END]") == string::npos) {
				if ((str.find("ampinabaqus") != string::npos)) { getline(in, str); linenum += 1; break; }//
				getline(in, str); linenum += 1;
				str0 = str;
				str = StringUtils::RemoveStrSpace(str0);
				str = StringUtils::StrToLower(str);
				if (StringUtils::IsCommentLine(str) || str.size() < 1) continue;
				//if (str.find("type=") != string::npos) {
				//	string substr = str.substr(str.find_first_of('=') + 1);
				//	if (substr.find("dirichlet") != string::npos && substr.length() == 9) {
				//		ampblock._BCTypeName = "dirichlet";
				//		ampblock._BCType = BCType::DIRICHLETBC;
				//		HasElmt = true;
				//	}
				//	else if (substr.find("neumann") != string::npos && substr.length() == 7) {
				//		ampblock._BCTypeName = "neumann";
				//		ampblock._BCType = BCType::NEUMANNBC;
				//		HasElmt = true;
				//	}
				//	else if (substr.find("nodaldirichlet") != string::npos && substr.length() == 14) {
				//		ampblock._BCTypeName = "nodaldirichlet";
				//		ampblock._BCType = BCType::NODALDIRICHLETBC;
				//		HasElmt = true;
				//	}
				//	else if (substr.find("nodalneumann") != string::npos && substr.length() == 12) {
				//		ampblock._BCTypeName = "nodalneumann";
				//		ampblock._BCType = BCType::NODALNEUMANNBC;
				//		HasElmt = true;
				//	}
				//	else if (substr.find("pressure") != string::npos && substr.length() == 8) {
				//		ampblock._BCTypeName = "pressure";
				//		ampblock._BCType = BCType::PRESSUREBC;
				//		HasElmt = true;
				//	}
				//	else if (substr.find("user") != string::npos) {
				//		number = StringUtils::SplitStrNum(str);
				//		if (number.size() < 1) {
				//			MessagePrinter::PrintErrorInLineNumber(linenum);
				//			msg = "no user number found umat in [" + ampblock._BCBlockName + "] sub block, 'type=user1,user2,...,user10' is expected";
				//			MessagePrinter::PrintErrorTxt(msg);
				//			MessagePrinter::AsFem_Exit();
				//			return false;
				//		}

				//		if (int(number[0]) < 1 || int(number[0]) > 10) {
				//			MessagePrinter::PrintErrorInLineNumber(linenum);
				//			msg = "user number is invalid in [" + ampblock._BCBlockName + "] sub block, 'type=user1,user2,...,user10' is expected";
				//			MessagePrinter::PrintErrorTxt(msg);
				//			MessagePrinter::AsFem_Exit();

				//			return false;
				//		}
				//		else {
				//			ampblock._BCTypeName = str.substr(str.find_first_of('=') + 1);
				//			switch (int(number[0])) {
				//			case 1:
				//				ampblock._BCType = BCType::USER1BC;
				//				break;
				//			case 2:
				//				ampblock._BCType = BCType::USER2BC;
				//				break;
				//			case 3:
				//				ampblock._BCType = BCType::USER3BC;
				//				break;
				//			case 4:
				//				ampblock._BCType = BCType::USER4BC;
				//				break;
				//			case 5:
				//				ampblock._BCType = BCType::USER5BC;
				//				break;
				//			case 6:
				//				ampblock._BCType = BCType::USER6BC;
				//				break;
				//			case 7:
				//				ampblock._BCType = BCType::USER7BC;
				//				break;
				//			case 8:
				//				ampblock._BCType = BCType::USER8BC;
				//				break;
				//			case 9:
				//				ampblock._BCType = BCType::USER9BC;
				//				break;
				//			case 10:
				//				ampblock._BCType = BCType::USER10BC;
				//				break;
				//			default:
				//				ampblock._BCType = BCType::NULLBC;
				//				break;
				//			}
				//			HasElmt = true;
				//		}
				//	}
				//	else {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "unsupported type in [" + ampblock._BCBlockName + "] sub block, 'type=" + substr + "' is invalid";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//	}
				//}
				//else if (str.find("dof=") != string::npos) {
				//	if (!HasElmt) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "type= is not given yet in [" + ampblock._BCBlockName + "] sub block, 'dof=' should be given after 'type=' in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//	}
				//	if (str.size() < 5) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "dof name not found in [" + ampblock._BCBlockName + "] sub block, 'dof=' should be given after 'type=' in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//	}
				//	else {
				//		ampblock._DofName = str.substr(4);
				//		if (!dofHandler.IsValidDofName(ampblock._DofName)) {
				//			MessagePrinter::PrintErrorInLineNumber(linenum);
				//			msg = "invalid dof name found in [" + ampblock._BCBlockName + "] sub block, dof must be one of the names list in [dofs] block";
				//			MessagePrinter::PrintErrorTxt(msg);
				//			MessagePrinter::AsFem_Exit();
				//			return false;
				//		}
				//		ampblock._DofID = dofHandler.GetDofIDviaDofName(ampblock._DofName);

				//		HasDof = true;
				//	}
				//}
				//else if (str.find("boundary=") != string::npos) {
				//	if (!HasElmt) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "'type=' is not given yet in [" + ampblock._BCBlockName + "] sub block, 'boundary=' should be given after 'type=' in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//	}
				//	if (str.size() < 10) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "boundary name not found in [" + ampblock._BCBlockName + "] sub block, 'boundary=bcname' should be given in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//		return false;
				//	}
				//	else {
				//		int i = str0.find_first_of('=');
				//		string substr = str0.substr(i + 1, str0.length());
				//		amplist = StringUtils::SplitStr(substr, ' ');
				//		if (amplist.size() > 0)
				//		{
				//			if (StringUtils::IsUniqueStrVec(amplist)) {
				//				ampblock._BoundaryNameList = amplist;
				//				HasBoundary = true;
				//			}
				//			else {
				//				MessagePrinter::PrintErrorInLineNumber(linenum);
				//				msg = "duplicated boundary name in [" + ampblock._BCBlockName + "] sub block";
				//				MessagePrinter::PrintErrorTxt(msg);
				//				MessagePrinter::AsFem_Exit();
				//			}
				//		}
				//		else {
				//			HasBoundary = false;
				//			MessagePrinter::PrintErrorInLineNumber(linenum);
				//			msg = "no boundary name found in [" + ampblock._BCBlockName + "] sub block, 'boundary=bc1 bc2 ...' is expected";
				//			MessagePrinter::PrintErrorTxt(msg);
				//			MessagePrinter::AsFem_Exit();
				//		}
				//	}
				//}
				//else if (str.find("value=") != string::npos) {
				//	if (!HasElmt) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "no 'type=' found in [" + ampblock._BCBlockName + "] sub block, 'value=' must be given after 'type=' in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//	}
				//	number = StringUtils::SplitStrNum(str);
				//	if (number.size() < 1 && str.find("t") == string::npos) {
				//		MessagePrinter::PrintErrorInLineNumber(linenum);
				//		msg = "no bc value found in [" + ampblock._BCBlockName + "] sub block, 'value=real value' should be given in [bcs] sub block";
				//		MessagePrinter::PrintErrorTxt(msg);
				//		MessagePrinter::AsFem_Exit();
				//		return false;
				//	}
				//	else {
				//		HasValue = true;
				//		int i;
				//		i = str.find_first_of("=");
				//		string substr = str.substr(i + 1);
				//		// cout<<"ampblock="<<ampblock._BCBlockName<<", str="<<str<<", substr="<<substr<<endl;
				//		if (substr.find("t") != string::npos&&substr.find("*") == string::npos) {
				//			ampblock._IsTimeDependent = true;
				//			ampblock._BCValue = 1.0;
				//		}
				//		else if (substr.find("-t") != string::npos&&substr.find("*") == string::npos) {
				//			ampblock._IsTimeDependent = true;
				//			ampblock._BCValue = -1.0;
				//		}
				//		else {
				//			if (StringUtils::IsValidExpression(substr)) {
				//				ampblock._BCValue = number[0];
				//				ampblock._IsTimeDependent = true;
				//			}
				//			else {
				//				ampblock._BCValue = number[0];
				//				ampblock._IsTimeDependent = false;
				//			}
				//		}
				//	}
				//}
				//HasAmpBlock = true; HasTypeName = true; HasParams = true; 
			}
			if (HasTypeName&&HasParams&&HasAmpBlock) {//完成一个Amp数据读取后，添加到ampSystem中
				if (abaFlag) {
					for (auto it : _AmpFromInp) {//
						_Index++;
						typestr.clear();
						ampblock.Init();
						ampblock._AmpIndex = _Index;
						ampblock._AmpBlockName = it.first[0];//name/type/smooth/fixInterval
						//HasAmpBlock = true;
						ampblock._AmpTypeName = it.first[1];//
						typestr = it.first[1];
						//HasTypeName = true;
						if (typestr.find("TABULAR") != string::npos) {
							ampblock._AmpType = AmpType::TABULAR;
						}
						else if (typestr.find("EQUALLY") != string::npos) {//EQUALLY SPACED, fixed interval=0.3
							ampblock._AmpType = AmpType::EQUALSPACED;
							ampblock._Smoothdata  = stod(it.first[2]);
							ampblock._FixInterval = stod(it.first[3]);
						}
						else if (typestr.find("PERIODIC") != string::npos) {
							ampblock._AmpType = AmpType::PERIODIC;
						}
						else if (typestr.find("MODULATED") != string::npos) {
							ampblock._AmpType = AmpType::MODULATED;
						}
						else if (typestr.find("DECAY") != string::npos) {
							ampblock._AmpType = AmpType::DECAY;
						}
						else if (typestr.find("SOLUTION") != string::npos) {//SOLUTION DEPENDENT
							ampblock._AmpType = AmpType::SOLUTIONBASED;
						}
						else if (typestr.find("SMOOTH") != string::npos) {//SMOOTH STEP
							ampblock._AmpType = AmpType::SMOOTHSTEP;
						}
						else if (typestr.find("SPECTRUM") != string::npos) {
							ampblock._AmpType = AmpType::SPECTRUM;
						}
						else if (typestr.find("PSD") != string::npos) {
							ampblock._AmpType = AmpType::PSD;
						}
						else if (typestr.find("USER") != string::npos) {//currently not defined !!
							ampblock._AmpType = AmpType::USER1A;
							msg = "Currently not supported Amp type ["typestr"] in [" + ampblock._AmpBlockName + "] sub block, please check your input !!";
							MessagePrinter::PrintErrorTxt(msg);
						}
						else {
							MessagePrinter::PrintErrorInLineNumber(linenum);
							msg = "not supported Amp type ["typestr"] in [" + ampblock._AmpBlockName + "] sub block, please check your input !!";
							MessagePrinter::PrintErrorTxt(msg);
							MessagePrinter::AsFem_Exit();
						}

						if (it.first.size()>=3) {// concerning smooth or fixIntercal  value
							//name type smooth fixInterval
							HasSmooth = true;
							ampblock._Smoothdata = stod(it.first[2]);
						}

						ampblock._AmpParams = it.second;
						//HasParams = true;
						//if (!HasValue) ampblock._Smoothdata = 0.0;
						ampSystem.AddAmpBlock2List(ampblock);
					}
				}
				else {// inner input; need to add
					HasAmpBlock = true;
					if (!HasSmooth) ampblock._Smoothdata = 0.0;
					ampSystem.AddAmpBlock2List(ampblock);
				}
			}
			else {
				msg = "information is not complete in [amps] sub block, some information is missing in [" + ampblock._AmpBlockName + "]";
				MessagePrinter::PrintErrorTxt(msg, false);
				if (!HasTypeName) {
					msg = "no type found in [" + ampblock._AmpBlockName + "] sub block, 'type=' must be given in each [amps] sub block";
					MessagePrinter::PrintErrorTxt(msg);
				}
				if (!HasParams) {
					msg = "no paramters found in [" + ampblock._AmpBlockName + "] sub block, 'params=param_data' must be given in each [amps] sub block";
					MessagePrinter::PrintErrorTxt(msg);
				}
				return false;
			}
		}
	}

	return true;
}
