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
//+++ Purpose: implement the lagrange shape function in 2d case, in
//+++          this code, we support:
//+++            quad4, quad8 and quad9 mesh
//+++            tri3 and tri6 mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.23  add some shp Ni for some future-may-need elements

#include "FE/LagrangeShapeFun.h"

void LagrangeShapeFun::Compute2DLagrangeShapeFun(const double &xi,const double &eta,const Nodes &nodes,bool flag){
    _DetJac=0.0;
    switch (GetMeshType()){
        case MeshType::QUAD4:{//N1-N4 start from left bottom corner��count as per anticlockwise, right-hand rule as normal direction
            (*this)(1,0)=(1.0-xi)*(1.0-eta)/4.0;
            (*this)(2,0)=(1.0+xi)*(1.0-eta)/4.0;
            (*this)(3,0)=(1.0+xi)*(1.0+eta)/4.0;
            (*this)(4,0)=(1.0-xi)*(1.0+eta)/4.0;

            (*this)(1,1)= (eta-1.0)/4.0;
            (*this)(1,2)= (xi -1.0)/4.0;

            (*this)(2,1)= (1.0-eta)/4.0;
            (*this)(2,2)=-(1.0+xi )/4.0;

            (*this)(3,1)= (1.0+eta)/4.0;
            (*this)(3,2)= (1.0+xi )/4.0;

            (*this)(4,1)=-(1.0+eta)/4.0;
            (*this)(4,2)= (1.0-xi )/4.0;
            break;
        }
        case MeshType::QUAD8:{
            // 2D-8Nodes rectangle element
            (*this)(1,0)=(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)/4.0;
            (*this)(2,0)=(1.0+xi)*(1.0-eta)*( xi-eta-1.0)/4.0;
            (*this)(3,0)=(1.0+xi)*(1.0+eta)*( xi+eta-1.0)/4.0;
            (*this)(4,0)=(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)/4.0;
            (*this)(5,0)=(1.0-xi*xi)*(1.0-eta    )/2.0;
            (*this)(6,0)=(1.0+xi   )*(1.0-eta*eta)/2.0;
            (*this)(7,0)=(1.0-xi*xi)*(1.0+eta    )/2.0;
            (*this)(8,0)=(1.0-xi   )*(1.0-eta*eta)/2.0;

            // derivatives over xi and eta
            (*this)(1,1)=(1.0-eta)*(2.0*xi+eta)/4.0;
            (*this)(1,2)=(1.0-xi )*(xi+2.0*eta)/4.0;

            (*this)(2,1)=(1.0-eta)*(2.0*xi-eta)/4.0;
            (*this)(2,2)=(1.0+xi )*(2.0*eta-xi)/4.0;

            (*this)(3,1)=(1.0+eta)*(2.0*xi+eta)/4.0;
            (*this)(3,2)=(1.0+xi )*(xi+2.0*eta)/4.0;

            (*this)(4,1)=(1.0+eta)*(2.0*xi-eta)/4.0;
            (*this)(4,2)=(1.0-xi )*(2.0*eta-xi)/4.0;

            (*this)(5,1)=xi*(eta-1.0);
            (*this)(5,2)=(xi*xi-1.0)/2.0;

            (*this)(6,1)=(1.0-eta*eta)/2.0;
            (*this)(6,2)=-(1.0+xi)*eta;

            (*this)(7,1)=-xi*(1.0+eta);
            (*this)(7,2)=(1.0-xi*xi)/2.0;

            (*this)(8,1)=(eta*eta-1.0)/2.0;
            (*this)(8,2)=(xi-1.0)*eta;
            break;
        }
        case MeshType::QUAD9:{
            // 2D-9Nodes rectangle element
            (*this)(1,0)=(xi*xi-xi )*(eta*eta-eta)/4.0;
            (*this)(2,0)=(xi*xi+xi )*(eta*eta-eta)/4.0;
            (*this)(3,0)=(xi*xi+xi )*(eta*eta+eta)/4.0;
            (*this)(4,0)=(xi*xi-xi )*(eta*eta+eta)/4.0;
            (*this)(5,0)=(1.0-xi*xi)*(eta*eta-eta)/2.0;
            (*this)(6,0)=(xi*xi+xi )*(1.0-eta*eta)/2.0;
            (*this)(7,0)=(1.0-xi*xi)*(eta*eta+eta)/2.0;
            (*this)(8,0)=(xi*xi-xi )*(1.0-eta*eta)/2.0;
            (*this)(9,0)=(1.0-xi*xi)*(1.0-eta*eta);

            (*this)(1,1)=(2.0*xi-1.0)*(eta*eta-eta)/4.0;
            (*this)(1,2)=(xi*xi-xi  )*(2.0*eta-1.0)/4.0;

            (*this)(2,1)=(2.0*xi+1.0)*(eta*eta-eta)/4.0;
            (*this)(2,2)=(xi*xi+xi  )*(2.0*eta-1.0)/4.0;

            (*this)(3,1)=(2.0*xi+1.0)*(eta*eta+eta)/4.0;
            (*this)(3,2)=(xi*xi+xi  )*(2.0*eta+1.0)/4.0;

            (*this)(4,1)=(2.0*xi-1.0)*(eta*eta+eta)/4.0;
            (*this)(4,2)=(xi*xi-xi  )*(2.0*eta+1.0)/4.0;

            (*this)(5,1)=-xi*(eta*eta-eta);
            (*this)(5,2)=(1.0-xi*xi )*(2.0*eta-1.0)/2.0;

            (*this)(6,1)=(2.0*xi+1.0)*(1.0-eta*eta)/2.0;
            (*this)(6,2)=-(xi*xi+xi )*eta;

            (*this)(7,1)=-xi*(eta*eta+eta);
            (*this)(7,2)=(1.0-xi*xi )*(2.0*eta+1.0)/2.0;

            (*this)(8,1)=(2.0*xi-1.0)*(1.0-eta*eta)/2.0;
            (*this)(8,2)=-(xi*xi-xi )*eta;

            (*this)(9,1)=-2.0*xi*(1.0-eta*eta);
            (*this)(9,2)=-2.0*eta*(1.0-xi*xi);
            break;
        }
        case MeshType::TRI3:{
            (*this)(1,0)=1.0-xi-eta;
            (*this)(1,1)=-1.0;
            (*this)(1,2)=-1.0;

            (*this)(2,0)=xi;
            (*this)(2,1)=1.0;
            (*this)(2,2)=0.0;

            (*this)(3,0)=eta;
            (*this)(3,1)=0.0;
            (*this)(3,2)=1.0;
            break;
        }
        case MeshType::TRI6:{
            // taken from: http://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
            (*this)(1,0)= (1.0-xi-eta)*(1.0-2*xi-2*eta);
            (*this)(1,1)=-3.0+4.0*eta+4.0*xi;
            (*this)(1,2)=-3.0+4.0*eta+4.0*xi;

            (*this)(2,0)= xi*(2.0*xi-1.0);
            (*this)(2,1)= 4.0*xi-1.0;
            (*this)(2,2)= 0.0;

            (*this)(3,0)= eta*(2.0*eta-1.0);
            (*this)(3,1)= 0.0;
            (*this)(3,2)= 4.0*eta-1.0;

            (*this)(4,0)= 4.0*xi*(1.0-xi-eta);
            (*this)(4,1)= 4.0*(1.0-eta-2*xi);
            (*this)(4,2)=-4.0*xi;

            (*this)(5,0)= 4.0*xi*eta;
            (*this)(5,1)= 4.0*eta;
            (*this)(5,2)= 4.0*xi;

            (*this)(6,0)= 4.0*eta*(1.0-xi-eta);
            (*this)(6,1)=-4.0*eta;
            (*this)(6,2)= 4.0*(1-2*eta-xi);
            
            break;
        }

		case MeshType::S3: {//	||Active degrees of freedom 1, 2, 3, 4, 5, 6

			break;
		}
		case MeshType::S3R: {

			break;
		}
		case MeshType::S4: {

			break;
		}
		case MeshType::S4R: {

			break;
		}
		case MeshType::M3D3: {//VFIFE	||Active degrees of freedom 1, 2, 3
			(*this)(1, 0) = 1.0 - xi - eta;
			(*this)(1, 1) = -1.0;
			(*this)(1, 2) = -1.0;

			(*this)(2, 0) = xi;
			(*this)(2, 1) = 1.0;
			(*this)(2, 2) = 0.0;

			(*this)(3, 0) = eta;
			(*this)(3, 1) = 0.0;
			(*this)(3, 2) = 1.0;
			break;
		}
		case MeshType::M3D4: {//	||Active degrees of freedom 1, 2, 3
			(*this)(1, 0) = (1.0 - xi)*(1.0 - eta) / 4.0;
			(*this)(2, 0) = (1.0 + xi)*(1.0 - eta) / 4.0;
			(*this)(3, 0) = (1.0 + xi)*(1.0 + eta) / 4.0;
			(*this)(4, 0) = (1.0 - xi)*(1.0 + eta) / 4.0;

			(*this)(1, 1) = (eta - 1.0) / 4.0;
			(*this)(1, 2) = (xi - 1.0) / 4.0;

			(*this)(2, 1) = (1.0 - eta) / 4.0;
			(*this)(2, 2) = -(1.0 + xi) / 4.0;

			(*this)(3, 1) = (1.0 + eta) / 4.0;
			(*this)(3, 2) = (1.0 + xi) / 4.0;

			(*this)(4, 1) = -(1.0 + eta) / 4.0;
			(*this)(4, 2) = (1.0 - xi) / 4.0;

			break;
		}
		case MeshType::SC6R: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::SFM3D3: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::SFM3D4: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::SFM3D4R: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::CPE3: {//	||Active degrees of freedom 1, 2

			break;
		}
		case MeshType::CPE4: {

			break;
		}
		case MeshType::CPE4R: {

			break;
		}
		case MeshType::CPE6M: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::CPS3: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
		case MeshType::CPS4: {

			break;
		}
		case MeshType::CPS4R: {

			break;
		}
		case MeshType::CPS6M: {//	||Active degrees of freedom 1, 2, 3

			break;
		}
        default:
            break;
    }

    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    for(int i=1;i<=GetShapeFunNums();i++){
        _dxdxi+=(*this)(i,1)*nodes(i,1);
        _dydxi+=(*this)(i,1)*nodes(i,2);
        _dzdxi+=(*this)(i,1)*nodes(i,3);

        _dxdeta+=(*this)(i,2)*nodes(i,1);
        _dydeta+=(*this)(i,2)*nodes(i,2);
        _dzdeta+=(*this)(i,2)*nodes(i,3);
    }
   
    _DetJac=(_dydxi*_dzdeta-_dydeta*_dzdxi)*(_dydxi*_dzdeta-_dydeta*_dzdxi)
           +(_dzdxi*_dxdeta-_dzdeta*_dxdxi)*(_dzdxi*_dxdeta-_dzdeta*_dxdxi)
           +(_dxdxi*_dydeta-_dxdeta*_dydxi)*(_dxdxi*_dydeta-_dxdeta*_dydxi);
    _DetJac=sqrt(_DetJac);

    if(abs(_DetJac)<1.0e-15){
        MessagePrinter::PrintErrorTxt("singular element in 2D case, this error occurs in your 2D shape function calculation");
        MessagePrinter::AsFem_Exit();
    }

    

    _Jac[0][0]= _dxdxi;_Jac[0][1]= _dydxi;
    _Jac[1][0]=_dxdeta;_Jac[1][1]=_dydeta;//�ſɱȾ���

    _XJac[0][0]= _Jac[1][1]/_DetJac;
    _XJac[0][1]=-_Jac[0][1]/_DetJac;
    _XJac[1][0]=-_Jac[1][0]/_DetJac;
    _XJac[1][1]= _Jac[0][0]/_DetJac;//	J����ʽ/|J|  -> ���� [J]^-1

    double temp;
    for(int i=1;i<=GetShapeFunNums();i++){
        if(flag){
            temp=(*this)(i,1)*_XJac[0][0]+(*this)(i,2)*_XJac[0][1];
            (*this)(i,2)=(*this)(i,1)*_XJac[1][0]+(*this)(i,2)*_XJac[1][1];
            (*this)(i,1)=temp;
        }
       _shape_value[i-1]=(*this)(i,0);//Ni
       _shape_grad[i-1].setZero();
       _shape_grad[i-1](1)=(*this)(i,1);// sum dNi/dx
       _shape_grad[i-1](2)=(*this)(i,2);// sum dNi/dy
	   //3D element ��Ӧ����dNi/dz
    }
}