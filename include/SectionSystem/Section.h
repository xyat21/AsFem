//****************************************************************
// Purpose: This file contains the class definition for Sections
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#ifndef SECTION_H_
#define SECTION_H_
#include <iostream>
#include <string>
#include <cmath>

// obviously it is ratio of the circumference of a circle to its diameter, ^_-
const double PI = 3.14159265358979323846;


class BaseSection
{
	/* Abstract Base class for Section Module; all class inherited from
	* BaseSection should implement a series method to calculation section
	* characteritics, see all the pure virtual functions.
	*
	*/
protected:
	std::string type_;          // section type name
	int id_;                    // section id
	double A_;                  // area
	double Sy_, Sz_;            // moment of area
	double CGy_, CGz_;          // shape center
	double SHy_, SHz_;          // shear center
	double Iyy_, Izz_, Iyz_;    // inertial moments and product

protected:
	virtual void calcArea() = 0;
	virtual void calcSy() = 0;
	virtual void calcSz() = 0;
	virtual void calcCGy() = 0;
	virtual void calcCGz() = 0;
	virtual void calcSHy() = 0;
	virtual void calcSHz() = 0;
	virtual void calcIyy() = 0;
	virtual void calcIzz() = 0;
	virtual void calcIyz() = 0;
	void initiate();

public:
	BaseSection(int id);
	virtual ~BaseSection();
	virtual void info() const;
	std::string type() const;
	int id() const;
	double A() const;
	double Sy() const;
	double Sz() const;
	double CGy() const;
	double CGz() const;
	double SHy() const;
	double SHz() const;
	double Iyy() const;
	double Izz() const;
	double Iyz() const;
};


struct CustomSectionParas
{
	/* struct for store section characteristics, see CustomSection class */
	double A;
	double Sy, Sz;
	double SHy, SHz;
	double CGy, CGz;
	double Iyy, Izz, Iyz;
};


class CustomSection : public BaseSection
{
	/* class for user defined section, all the section characteristics should
	* be inputed by users by using CustomSectionParas struct.
	* Parameter lists:
	*   A - area
	*   Sy, Sz - Moment of area of y and z axis, respectively.
	*   SHy, SHz - Shear center
	*   CGy, CGz - shape center
	*   Iyy, Izz - moment of inertia
	*   Iyz - inertial product
	*/
protected:
	virtual void calcArea();
	virtual void calcSy();
	virtual void calcSz();
	virtual void calcCGy();
	virtual void calcCGz();
	virtual void calcSHy();
	virtual void calcSHz();
	virtual void calcIyy();
	virtual void calcIzz();
	virtual void calcIyz();

public:
	CustomSection(int id, CustomSectionParas &paras);
	virtual ~CustomSection();
	virtual void info() const;
};


class Ctube : public BaseSection
{
	/* circular tube section
	*  ri - internal diameter
	*  ro - outer diameter
	*
	*           ^ y
	*           |
	*        ...|...
	*       . . | . .
	*      . .  |  . .
	* ----.-.---|---.-.---------> z
	*      . .  |  . .
	*       . ..|.. .
	*        ...|...
	*           |
	*/
private:
	double ri_;
	double ro_;

protected:
	virtual void calcArea();
	virtual void calcSy();
	virtual void calcSz();
	virtual void calcCGy();
	virtual void calcCGz();
	virtual void calcSHy();
	virtual void calcSHz();
	virtual void calcIyy();
	virtual void calcIzz();
	virtual void calcIyz();

public:
	Ctube(int id, double ri, double ro);
	virtual ~Ctube();
	double ri() const;
	double ro() const;
	virtual void info() const;
};


class AngleSteel : public BaseSection
{
	/* Section Parameters
	* Sy, Sz - calculated in xOy coordinate system
	* CGy, CGZ - presented in xOy coordiante system
	* Iyy, Izz - calculated in xcOyc coordinate system (principle central axis)
	*
	*     y ^ t1
	*       |____
	*       |   |
	*       |   |
	*       |   |    ^ yc
	*  w1   |   |    |
	*       |A1 |    |
	*       |   |    |
	*       |   | ---|-----> zc
	*       |   |    |shape center
	*       |   |______________
	*       |   |     A2      | t2
	*   ____|___|_____________|________> z
	*       |
	*       |        w2
	*
	*/
private:
	double w1_;
	double w2_;
	double t1_;
	double t2_;

protected:
	virtual void calcArea();
	virtual void calcSy();
	virtual void calcSz();
	virtual void calcCGy();
	virtual void calcCGz();
	virtual void calcSHy();
	virtual void calcSHz();
	virtual void calcIyy();
	virtual void calcIzz();
	virtual void calcIyz();

public:
	AngleSteel(int id, double w1, double w2, double t1, double t2);
	virtual ~AngleSteel();
	double w1() const;
	double w2() const;
	double t1() const;
	double t2() const;
	virtual void info() const;
};


class Rectangle : public BaseSection
{
	/* Seciton Parameters
	* xOy coordinate system is same as xcOyc;
	*
	*
	*              ^ y
	*        ______|______
	*        |     |     |
	*        |     |     |
	*        |     |     |  height
	*    ----|-----|-----|-------> z
	*        |     |     |
	*        |     |     |
	*        |_____|_____|
	*              |
	*           width
	*
	*/
protected:
	double width_;
	double height_;

protected:
	virtual void calcArea();
	virtual void calcSy();
	virtual void calcSz();
	virtual void calcCGy();
	virtual void calcCGz();
	virtual void calcSHy();
	virtual void calcSHz();
	virtual void calcIyy();
	virtual void calcIzz();
	virtual void calcIyz();

public:
	Rectangle(int id, double width, double height);
	virtual ~Rectangle();
	double width() const;
	double height() const;
	virtual void info() const;
};


#endif
