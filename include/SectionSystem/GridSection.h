/* ************************************************************************* **
**    OpenVFIFE - Open System for Vector Form Instrinsic                     **
**                Finite Element Method (VFIFE)                              **
**                GinkGo(Tan Biao)                                           **
**                                                                           **
**                                                                           **
** (C) Copyright 2021, The GinkGo(Tan Biao). All Rights Reserved.            **
**                                                                           **
** Commercial use of this program without express permission of              **
** GinkGo(Tan Biao), is strictly prohibited.  See                            **
** file 'COPYRIGHT'  in main directory for information on usage and          **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.                  **
**                                                                           **
** Developed by:                                                             **
**      Tan Biao (ginkgoltd@outlook.com)                                     **
**                                                                           **
** ************************************************************************* */

// $Date: 2021-06-07$
// Written: Tan Biao
// Revised:
//
// Purpose: This file contains the class definition for SectionGrid, including
// SectionCell, SectionResponse, RectangleGrid, AngleSteelGrid

// The interface:
//

#ifndef _GRIDSECTION_H_
#define _GRIDSECTION_H_
#include <vector>
#include "section.h"
#include "material.h"
#include "integrator.h"// should modify to replace original headers.


// define a struct to restore point coordinate
struct SectionPoint { double y, z; };



struct CellNodes
{
	/* store cell corner nodes, WARNING: Please follow the instructions of node
	*  sequence. Again, FOLLOW THE SEQUENCE AS ILLUSTRAED IN THE FIGURE.
	*  ^
	*  |
	*  |    p1 _______ p2
	*  |       |     |
	*  |       | * * |
	*  |       |     |
	*  |       | * * |
	*  |    p3 |_____| p4
	*  |
	*  |------------> z
	*/
	SectionPoint p1, p2, p3, p4;
};


struct SectionGridParameter
{
	/* Parameters for section grid generation
	* sections - number of sections
	* cells - number of cells
	* points - number of integration points of each cell(along one direction)
	* x1, x2, ..., x8 - paramers for grid generation for different section
	*
	* Retangle Section: only x1, x2 are used
	*   x1 - number of cells along section width
	*   x2 - number of cells along section height
	*
	* AngleSteel Section: only x1, x2 are used
	*   x1 - number of cells along vertical leg (including corner)
	*   x2 - number of cells along horizontal leg
	*
	* You can use x1, ..., x8 for other types of Section which is regular.
	* However, Please add a explanation here for other users.
	*/
	int sections;       // number of sections
	int cells;          // number of cells
	int points;         // number of integration points of each cell
	int x1, x2, x3, x4, x5, x6, x7, x8; // control paramers
};


class SectionCell
{
	/*
	* only quadrilateral cell has been implemented here. If you want to use
	* other type of SectionCell such as triangle, please add new type of
	* SectionCell
	*
	*  ^ y
	*  |
	*  |    p1 _______ p2             * integration points
	*  |       |     |               -------
	*  |       | * * |               |     |  section cell
	*  |       |     |               |     |
	*  |       | * * |               |_____|
	*  |    p3 |_____| p4
	*  |
	*  |------------> z
	*
	* TODO:
	* Notes: The SectionCell class only support same material, however it is
	*        quite plain to assign a different type of material to each
	*        integration point. Just modified material_ and add a function for
	*        assign materials. and change SectionRespose::setCellMaterial()
	*/
private:
	int number_integral_points_;// number of integral points along one axis
	int fractured_fiber;
	const BaseMaterial* material_;
	const GaussLegendre* legendre_;
	inline void setIntegratePointCoordinate();

public:
	// record coordinates of integral points
	std::vector<std::vector<SectionPoint>> integral_points;
	SectionPoint cp1, cp2, cp3, cp4;   // cell node coordinates
	double area;            // cell area

	// record stress and strain at each integral points
	std::vector<std::vector<double*>> mstress; // stress - current stress
	std::vector<std::vector<double*>> mstrain; // strain - current strain
	std::vector<std::vector<double*>> mps; // ps - current plastic strain
	std::vector<std::vector<double*>> malpha;  // alpha - back stress
	std::vector<std::vector<double*>> mq;      // q - hardening funciton

	// restore degree of damage of cell
	std::vector<std::vector<bool>> damage;     // if fiber fractured or not
	double damage_degree;

	// record total internal force, [N, miy, miz, mjy, mjz]
	std::vector<double> force;

public:
	// points - integral points
	SectionCell(const BaseMaterial* mat, const CellNodes &cn, int p = 2);
	SectionCell(const CellNodes &cn, int p = 2);
	virtual ~SectionCell();

	inline void setMaterial(const BaseMaterial* mat);
	virtual void response(const Eigen::RowVectorXd &transfer,
		const Eigen::VectorXd &deform);
};


class SectionResponse
{
	/* Section is composed by cells, see RectangleGrid. SectionResponse class
	* is responsible for gathering cells results.
	*
	* damageDegree() is trying to calculate the damage degree of the section,
	* which is defined as follow:
	*   $ \lambda = \frac{N_{collapse}}{N} $
	* $ N_{collapse} $ - collapsed integration points
	* $ N $ - all integration points
	*
	*/
public:
	std::vector<SectionCell*> cells;
	std::vector<double> force;

public:
	virtual ~SectionResponse();
	virtual void setCellMaterial(const BaseMaterial* mat);
	virtual double damageDegree() const;
	virtual void response(const Eigen::RowVectorXd &transfer,
		const Eigen::VectorXd &deform);
};


class SectionGrid
{
	// interface of GridSection
public:
	virtual void initiateResponseContainers(const BaseSection* sect,
		SectionResponse* s) const = 0;
};


class RectangleGrid : public SectionGrid
{
	// only support same material type
	/* Grid Parameters
	* x1 - number of cells along vertical direction (y)
	* x2 - number of cells along horizontal direction (z)
	*
	*              ^ y
	*        ______|______
	*        |     |     |
	*        |  1  |  2  |
	*        |     |     |  height
	*    ----|-----------|-------> z
	*        |     |     |
	*        |  3  |  4  |
	*        |_____|_____|
	*              |
	*           width
	*
	*/

private:
	// rows - number of section cells along width
	// cols - number of section cells along height
	int rows_, cols_;
	int integrals_points_;

	inline void checkSectionType(const BaseSection* sect) const;

public:
	RectangleGrid(const SectionGridParameter &p);
	virtual ~RectangleGrid();

	virtual void initiateResponseContainers(const BaseSection* sect,
		SectionResponse* s) const;
};


class AngleSteelGrid : public SectionGrid
{
	/* Grid Parameters
	* x1 - number of cells along vertical direction (y), exclude the corner
	* x2 - number of cells along horizontal direction (z), exclude the corner
	*
	*         t1    ^ y
	*       _____   |
	*       |   |   |
	*       | 1 |   |    1______2
	*       |___|   |     |    |
	*  w1   |   |   |     |    |
	*       | 2 |   |    3|____|4
	*       |___|   |
	*       |   |   |
	* ------|-3-|---|---------------------> z
	*       |___|___|______________
	*       |   |   1  |  2  |  3  |
	*       |___|___|__|_____|_____|  t2
	*               |
	*               | w2
	*
	*/
private:
	int x1_, x2_;
	int integrals_points_;

	inline void checkSectionType(const BaseSection* sect) const;

public:
	AngleSteelGrid(const SectionGridParameter &para);
	virtual ~AngleSteelGrid();

	virtual void initiateResponseContainers(const BaseSection* sect,
		SectionResponse* s) const;
};


// creator of GridSection
class CreateGrid
{
	/* C++ has no reflection mechanism, which means the program can not create
	* the GridSection by the section type. CreateGridSection gives a very
	* simple implementation for it. If you create a new GridSetion class,
	* Please add your
	*/
public:
	static SectionGrid* getInstance(const BaseSection* sect,
		const SectionGridParameter &para)
	{
		if (sect->type() == "Rectangle")
		{
			// return new GridRectangle(x1);
			return new RectangleGrid(para);
		}
		else if (sect->type() == "AngleSteel")
		{
			return new AngleSteelGrid(para);
		}
		else
		{
			std::cerr << "Error: Unrecognized Section Type ";
			std::cerr << sect->type() << "." << std::endl;
			return nullptr;
		}
	}
};


#endif