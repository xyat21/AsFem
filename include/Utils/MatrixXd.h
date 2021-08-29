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
//+++ Date   : 2020.10.18
//+++ Reviewer : Xiaoyuan @ 2021.08.20
//+++ Purpose: Define the general Matrix  in AsFem
//+++          we mainly use this for the calculation of jacobian
//+++          If one wants to use Eigen's MatrixXd, please use
//+++          Eigen::MatrixXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Date   : 2021.08.28  add some new operation funs

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

#include "Eigen/Eigen"

#include "Utils/MessagePrinter.h"
#include "Utils/VectorXd.h"

using namespace std;


class MatrixXd{
public:
    MatrixXd();
    MatrixXd(const MatrixXd &a);
    MatrixXd(const int &m,const int &n);
    MatrixXd(const int &m,const int &n,const double &val);
    void Resize(const int &m,const int &n){
        _vals.resize(m*n,0.0);_M=m;_N=n;_MN=m*n;
    }
    void Resize(const int &m,const int &n,const double &val){
        _vals.resize(m*n,val);_M=m;_N=n;_MN=m*n;
    }
    double* GetDataPtr(){
        return _vals.data();
    }
    inline int GetM()const{return _M;}
    inline int GetN()const{return _N;}
    void Clean(){_vals.clear();}
    //*****************************************
    //*** Operator overload
    //*****************************************
    inline double& operator()(const int &i,const int &j){
        return _vals[(i-1)*_N+j-1];
    }
    inline double operator()(const int &i,const int &j)const{
        return _vals[(i-1)*_N+j-1];
    }
    inline double& operator[](const int &i){
        return _vals[i-1];
    }
    inline double operator[](const int &i)const{
        return _vals[i-1];
    }
    //*****************************************
    //*** For basic mathematic operator
    //*****************************************
    //*** for =
    inline MatrixXd& operator=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=val;
        return *this;
    }
    inline MatrixXd& operator=(const MatrixXd &a){
        if(_M==0&&_N==0){
            _M=a.GetM();_N=a.GetN();
            _MN=_M*_N;_vals.resize(_MN,0.0);
            for(int i=0;i<_MN;++i) _vals[i]=a._vals[i];
            return *this;
        }
        else{
            if(_M!=a.GetM()&&_N!=a.GetN()){
                MessagePrinter::PrintErrorTxt("a=b cant be applied for two matrix with different size");
                MessagePrinter::AsFem_Exit();
            }
            else{
                for(int i=0;i<_MN;++i) _vals[i]=a._vals[i];
                return *this;
            }
        }
        return *this;
    }
    //****************************
    //*** for +
    inline MatrixXd operator+(const double &val)const{
        MatrixXd temp(_M,_N);
        for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]+val;
        return temp;
    }
    inline MatrixXd operator+(const MatrixXd &a)const{
        MatrixXd temp(_M,_N);
        if(_M!=a.GetM()&&_N!=a.GetN()){
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two matrix with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]+a._vals[i];
            return temp;
        }
        return temp;
    }
    //*** for +=
    inline MatrixXd& operator+=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=_vals[i]+val;
        return *this;
    }
    inline MatrixXd& operator+=(const MatrixXd &a){
        MatrixXd temp(_M,_N);
        if(_M!=a.GetM()&&_N!=a.GetN()){
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two matrix with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_MN;++i) _vals[i]=_vals[i]+a._vals[i];
            return *this;
        }
        return *this;
    }
    //****************************
    //*** for -
    inline MatrixXd operator-(const double &val)const{
        MatrixXd temp(_M,_N);
        for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]-val;
        return temp;
    }
    inline MatrixXd operator-(const MatrixXd &a)const{
        MatrixXd temp(_M,_N);
        if(_M!=a.GetM()&&_N!=a.GetN()){
            MessagePrinter::PrintErrorTxt("a-b cant be applied for two matrix with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]-a._vals[i];
            return temp;
        }
        return temp;
    }
    //*** for -=
    inline MatrixXd& operator-=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=_vals[i]-val;
        return *this;
    }
    inline MatrixXd& operator-=(const MatrixXd &a){
        MatrixXd temp(_M,_N);
        if(_M!=a.GetM()&&_N!=a.GetN()){
            MessagePrinter::PrintErrorTxt("a-b cant be applied for two matrix with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_MN;++i) _vals[i]=_vals[i]-a._vals[i];
            return *this;
        }
        return *this;
    }
    //****************************
    //*** for *
    inline MatrixXd operator*(const double &val)const{
        MatrixXd temp(_M,_N);
        for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]*val;
        return temp;
    }
    inline VectorXd operator*(const VectorXd &a)const{
        VectorXd temp(_M,0.0);
        if(_N!=a.GetM()){
            MessagePrinter::PrintErrorTxt("A*b should be applied for A matrix has the same cols as b vector!");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=1;i<=_M;i++){
                temp(i)=0.0;
                for(int j=1;j<=_N;j++){
                    temp(i)+=(*this)(i,j)*a(j);
                }
            }
            return temp;
        }
        return temp;
    }

	inline Vector3d operator*(const Vector3d &a)const {
		Vector3d temp(0.0);
		if (_N != 3 || _M != 3) {
			MessagePrinter::PrintErrorTxt("A*b should be applied for A matrix 3x3 has the same cols as b vector 1x3 !");
			MessagePrinter::AsFem_Exit();
		}
		else {
			for (int i = 1; i <= 3; i++) {
				temp(i) = 0.0;
				for (int j = 1; j <= _M; j++) {
					temp(i) += (*this)(i, j)*a(j);
				}
			}
			return temp;
		}
		return temp;
	}

	inline MatrixXd operator*(const MatrixXd &a)const{
        MatrixXd temp(_M,a.GetN());
        if(_N!=a.GetM()){
            MessagePrinter::PrintErrorTxt("A*B should be applied for A matrix has the same cols as the rows of B matrix!");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=1;i<=_M;i++){
                for(int j=1;j<=a.GetN();j++){
                    temp(i,j)=0.0;
                    for(int k=1;k<=a.GetM();k++){
                        temp(i,j)+=(*this)(i,k)*a(k,j);
                    }
                }
            }
            return temp;
        }
        return temp;
    }
    //*** for *=
    inline MatrixXd& operator*=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=_vals[i]*val;
        return *this;
    }
    //****************************
    //*** for /
    inline MatrixXd operator/(const double &val)const{
        MatrixXd temp(_M,_N);
        for(int i=0;i<_MN;++i) temp._vals[i]=_vals[i]/val;
        return temp;
    }
    //*** for /=
    inline MatrixXd& operator/=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=_vals[i]/val;
        return *this;
    }

    void setZero(){
        for(int i=0;i<_MN;++i) _vals[i]=0.0;
    }
    void setRandom(){
        srand(time(0));
        for(int i=0;i<_MN;++i) _vals[i]=static_cast<double>(1.0*rand()/RAND_MAX);
    }
    inline MatrixXd Inverse()const{
        Eigen::MatrixXd Mat(_M,_N),MatInv(_M,_N);
        MatrixXd temp(_M,_N);

        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        MatInv=Mat.inverse();
        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                temp(i,j)=MatInv.coeff(i-1,j-1);
            }
        }
        return temp;
    }
    inline double Det()const{
        Eigen::MatrixXd Mat(_M,_N);
        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        return Mat.determinant();
    }

	inline MatrixXd Transpose()const {
		MatrixXd temp(_N, _M);
		for (int i = 1; i <= _M; i++) {
			for (int j = 1; j <= _N; j++) {
				temp(j, i) = (*this)(i, j);
			}
		}
		return temp;
	}

	//*** for column and row based operator
	inline Vector3d IthRow(const int &i)const {
		if (_N != 3) {
			MessagePrinter::PrintErrorTxt("A IthRow should be applied for A matrix _Mx3 has 3 cols !");
			MessagePrinter::AsFem_Exit();
		}
		Vector3d temp(0.0);
		temp(1) = (*this)(i, 1);
		temp(2) = (*this)(i, 2);
		temp(3) = (*this)(i, 3);
		return temp;
	}
	inline Vector3d IthCol(const int &i)const {
		if (_M != 3) {
			MessagePrinter::PrintErrorTxt("A IthCol should be applied for A matrix 3x_N has 3 rows !");
			MessagePrinter::AsFem_Exit();
		}
		Vector3d temp(0.0);
		temp(1) = (*this)(1, i);
		temp(2) = (*this)(2, i);
		temp(3) = (*this)(3, i);
		return temp;
	}

private:
    vector<double> _vals;
    int _M,_N,_MN;
};