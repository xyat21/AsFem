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

#include "SectionSystem/GridSection.h"


using namespace std;


/*****************************************************************************/
/*                                  Helpers                                  */
/*****************************************************************************/
inline void initVector2D(vector<vector<double*>> &v, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		vector<double*> tmp;
		for (int j = 0; j < cols; j++)
		{
			double* val = new double;
			*val = 0;
			tmp.push_back(val);
		}
		v.push_back(tmp);
	}
}

inline void initCellDamage(vector<vector<bool>> &v, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		vector<bool> tmp;
		for (int j = 0; j < cols; j++)
		{
			tmp.push_back(false);
		}
		v.push_back(tmp);
	}
}

// for test
void printVector2D(const vector<vector<double*>> &v)
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < v[0].size(); j++)
		{
			cout << *(v[i][j]) << ", ";
		}
		cout << endl;
	}
	cout << endl;
};

void printVector2D(const vector<vector<double>> &v)
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < v[0].size(); j++)
		{
			cout << v[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;
};


/*****************************************************************************/
/*                                 SectionCell                               */
/*****************************************************************************/
inline void SectionCell::setIntegratePointCoordinate()
{
	// initiate coordiantes of integral points
	// [[(y1, z1), (y1, z2),... (y1, zn)],
	//  [(y2, z1), (y2, z2),... (y2, zn)],
	//  ...
	//  [(yn, z1), (yn, z2),... (yn, zn)]]
	for (int i = 0; i < number_integral_points_; i++)
	{
		vector<SectionPoint> tmp;
		for (int j = 0; j < number_integral_points_; j++)
		{
			double y = (cp1.y - cp3.y) / 2.0 * legendre_->x(i) +
				(cp1.y + cp3.y) / 2.0;
			double z = (cp2.z - cp1.z) / 2.0 * legendre_->x(j) +
				(cp1.z + cp2.z) / 2.0;
			SectionPoint p = { .y = y,.z = z };
			tmp.push_back(p);
		}
		integral_points.push_back(tmp);
	}
}

SectionCell::SectionCell(const BaseMaterial* mat, const CellNodes &cn, int p)
{
	// set material
	material_ = mat;

	// initiate integrator
	if (p <= 0)
	{
		cerr << "Warning: integral points must be positive." << endl;
	}
	number_integral_points_ = p;
	legendre_ = new GaussLegendre(p);

	// initiate cell nodes
	cp1 = cn.p1;
	cp2 = cn.p2;
	cp3 = cn.p3;
	cp4 = cn.p4;
	area = fabs((cn.p1.y - cn.p3.y) * (cn.p1.z - cn.p2.z));

	// initiate integral_points, after setting number_integral_points_
	setIntegratePointCoordinate();

	// initiate stress and strain at each integral points
	initVector2D(mstress, p, p);
	initVector2D(mstrain, p, p);
	initVector2D(mps, p, p);
	initVector2D(malpha, p, p);
	initVector2D(mq, p, p);
	initCellDamage(damage, p, p);
	fractured_fiber = 0;

	// initiate force container
	for (int i = 0; i < 5; i++)
	{
		force.push_back(0);
	}
}

SectionCell::SectionCell(const CellNodes &cn, int p)
{
	// initiate integrator
	if (p <= 0)
	{
		cerr << "Warning: integral points must be positive." << endl;
	}
	number_integral_points_ = p;
	legendre_ = new GaussLegendre(p);

	// initiate cell nodes
	cp1 = cn.p1;
	cp2 = cn.p2;
	cp3 = cn.p3;
	cp4 = cn.p4;
	area = fabs((cn.p1.y - cn.p3.y) * (cn.p1.z - cn.p2.z));

	// initiate integral_points, after setting number_integral_points_
	setIntegratePointCoordinate();

	// initiate stress and strain at each integral points
	initVector2D(mstress, p, p);
	initVector2D(mstrain, p, p);
	initVector2D(mps, p, p);
	initVector2D(malpha, p, p);
	initVector2D(mq, p, p);
	initCellDamage(damage, p, p);
	fractured_fiber = 0;

	// initiate force container
	for (int i = 0; i < 5; i++)
	{
		force.push_back(0);
	}
}

SectionCell::~SectionCell()
{
	// release materials_ container in case of memory leakages
	for (int i = 0; i < mstress.size(); i++)
	{
		for (int j = 0; j < mstress[0].size(); j++)
		{
			// release mstress_
			delete mstress[i][j];
			mstress[i][j] = NULL;

			// release mstrain_
			delete mstrain[i][j];
			mstrain[i][j] = NULL;

			// release mps_
			delete mps[i][j];
			mps[i][j] = NULL;

			// release malpha_
			delete malpha[i][j];
			malpha[i][j] = NULL;

			// release mq_
			delete mq[i][j];
			mq[i][j] = NULL;
		}
	}

	// clear mstress_
	mstress.clear();
	mstress.shrink_to_fit();

	// clear mstrain_
	mstrain.clear();
	mstrain.shrink_to_fit();

	// clear mps_
	mps.clear();
	mps.shrink_to_fit();

	// clear malpha_
	malpha.clear();
	malpha.shrink_to_fit();

	// clear mq_
	mq.clear();
	mq.shrink_to_fit();

	// release GaussLegendre integrator
	delete legendre_;
}

inline void SectionCell::setMaterial(const BaseMaterial* mat)
{
	material_ = mat;
}

void SectionCell::response(const Eigen::RowVectorXd &transfer,
	const Eigen::VectorXd &deform)
{
	// step 1: calculate strain increment and response of each integral points
	for (int i = 0; i < number_integral_points_; i++)
	{
		for (int j = 0; j < number_integral_points_; j++)
		{
			// compute the strain at each integral points
			double ds = transfer(0) * deform(0)
				+ transfer(1) * deform(1) * integral_points[i][j].z
				+ transfer(2) * deform(2) * integral_points[i][j].y
				+ transfer(3) * deform(3) * integral_points[i][j].z
				+ transfer(4) * deform(4) * integral_points[i][j].y;

			if (damage[i][j])
			{
				// fiber already fractured
				*(mstress[i][j]) = 0;
			}
			else
			{
				// fiber still working
				if (material_->isFractured(*(mstrain[i][j]) + ds))
				{
					// fiber fracture
					damage[i][j] = true;
					*(mstress[i][j]) = 0;
					fractured_fiber += 1;
				}
				else
				{
					// get material response of each integral points
					double eq = material_->Eq(ds, mstress[i][j], mstrain[i][j],
						mps[i][j], malpha[i][j], mq[i][j]);
				}
			}
		}
	}
	damage_degree = fractured_fiber / number_integral_points_
		/ number_integral_points_;

	// step 2: integration of cell stress
	// reset SectionResponse.force
	for (int i = 0; i < 5; i++)
	{
		force[i] = 0;   // fjx, miy, miz, mjy, mjz
	}

	double factor = (cp1.y - cp3.y) * (cp2.z - cp1.z) / 4.0;
	for (int i = 0; i < number_integral_points_; i++)
	{
		for (int j = 0; j < number_integral_points_; j++)
		{
			force[0] += legendre_->w(i) * legendre_->w(j) * transfer(0)
				* factor * (*(mstress[i][j]));
			force[1] += legendre_->w(i) * legendre_->w(j) * transfer(1)
				* integral_points[i][j].z * factor * (*(mstress[i][j]));
			force[2] += legendre_->w(i) * legendre_->w(j) * transfer(2)
				* integral_points[i][j].y * factor * (*(mstress[i][j]));
			force[3] += legendre_->w(i) * legendre_->w(j) * transfer(3)
				* integral_points[i][j].z * factor * (*(mstress[i][j]));
			force[4] += legendre_->w(i) * legendre_->w(j) * transfer(4)
				* integral_points[i][j].y * factor * (*(mstress[i][j]));
		}
	}
}


/*****************************************************************************/
/*                               SectionResponse                             */
/*****************************************************************************/
SectionResponse::~SectionResponse()
{
	// release cells container in case of memory leakages
	for (int i = 0; i < cells.size(); i++)
	{
		delete cells[i];
		cells[i] = nullptr;
	}

	// clear cells
	cells.clear();
	cells.shrink_to_fit();
}

void SectionResponse::setCellMaterial(const BaseMaterial* mat)
{
	for (int i = 0; i < cells.size(); i++)
	{
		cells[i]->setMaterial(mat);
	}
}

double SectionResponse::damageDegree() const
{
	double damage = 0;
	for (int i = 0; i < cells.size(); i++)
	{
		damage += cells[i]->damage_degree;
	}
	return damage / cells.size();
}

void SectionResponse::response(const Eigen::RowVectorXd &transfer,
	const Eigen::VectorXd &deform)
{
	// step 1: calculate cell response
	for (int i = 0; i < cells.size(); i++)
	{
		cells[i]->response(transfer, deform);
	}

	// clear section force
	for (int i = 0; i < force.size(); i++)
	{
		force[i] = 0;
	}

	for (int i = 0; i < force.size(); i++)
	{
		for (int j = 0; j < cells.size(); j++)
		{
			force[i] += cells[j]->force[i];
		}
	}
}


/*****************************************************************************/
/*                                 GridSection                               */
/*****************************************************************************/
// interface


/*****************************************************************************/
/*                               GridRectangle                               */
/*****************************************************************************/
RectangleGrid::RectangleGrid(const SectionGridParameter &p) : SectionGrid()
{
	// check rows and cols
	if ((p.x1 < 0) || (p.x2 < 0))
	{
		cerr << "Warning: cells must be positive." << endl;
	}
	if (p.x1 * p.x2 != p.cells)
	{
		cerr << "Warning: number of cells does not satisfy the grid, plase ";
		cerr << "check the SectionGridParameter. " << endl;
	}

	rows_ = p.x1;
	cols_ = p.x2;
	integrals_points_ = p.points;
}

RectangleGrid::~RectangleGrid()
{

}

inline void RectangleGrid::checkSectionType(const BaseSection* sect) const
{
	// check section type
	if (sect->type() != "Rectangle")
	{
		cerr << "Error: Please Enter an Rectangle section object." << endl;
		exit(-1);
	}
}

void RectangleGrid::initiateResponseContainers(const BaseSection* sect,
	SectionResponse* s) const
{
	// check section type first
	checkSectionType(sect);

	// coordinate of cell nodes
	double w = ((Rectangle*)sect)->width();
	double h = ((Rectangle*)sect)->height();

	vector<vector<SectionPoint>> sp;
	for (int i = 0; i <= rows_; i++)
	{
		vector<SectionPoint> tmp;
		for (int j = 0; j <= cols_; j++)
		{
			double y = h / 2.0 - h / cols_ * i;
			double z = -w / 2.0 + w / rows_ * j;
			SectionPoint p = { .y = y,.z = z };
			tmp.push_back(p);
		}
		sp.push_back(tmp);
	}

	// create cells
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			CellNodes cn = { .p1 = sp[i][j],.p2 = sp[i][j + 1],.p3 = sp[i + 1][j],
							.p4 = sp[i + 1][j + 1] };
			SectionCell* c = new SectionCell(cn, integrals_points_);
			cout << c->cp1.y << endl;
			s->cells.push_back(c);
		}
	}

	// initiate force container
	for (int i = 0; i < 5; i++)
	{
		s->force.push_back(0);
	}
}


/*****************************************************************************/
/*                              GridAngleSteel                               */
/*****************************************************************************/
AngleSteelGrid::AngleSteelGrid(const SectionGridParameter &para) :
	SectionGrid()
{
	// check rows and cols
	if ((para.x1 < 0) || (para.x2 < 0))
	{
		cerr << "Warning: cells must be positive." << endl;
	}
	if ((para.x1 + para.x2) != para.cells)
	{
		cerr << "Warning: number of cells does not satisfy the grid, plase ";
		cerr << "check the SectionGridParameter. " << endl;
	}

	x1_ = para.x1;
	x2_ = para.x2;
	integrals_points_ = para.points;
}

AngleSteelGrid::~AngleSteelGrid()
{

}

inline void AngleSteelGrid::checkSectionType(const BaseSection* sect) const
{
	// check section type
	if (sect->type() != "AngleSteel")
	{
		cerr << "Error: Please Enter an AngleSteel section object." << endl;
		exit(-1);
	}
}

void AngleSteelGrid::initiateResponseContainers(const BaseSection* sect,
	SectionResponse* s) const
{
	// check section type first
	checkSectionType(sect);

	// coordinate of cell nodes, multiple 10e6 for eliminate truncation error
	double w1 = ((AngleSteel*)sect)->w1();
	double w2 = ((AngleSteel*)sect)->w2();
	double t1 = ((AngleSteel*)sect)->t1();
	double t2 = ((AngleSteel*)sect)->t2();
	double yc = ((AngleSteel*)sect)->CGy();
	double zc = ((AngleSteel*)sect)->CGz();

	// vertical part
	for (int i = 0; i < x1_; i++)
	{
		// point 1
		double y1 = w1 - yc - (w1 - t2) / x1_ * i;
		double z1 = -zc;
		SectionPoint p1 = { .y = y1,.z = z1 };

		// point 2
		double y2 = y1;
		double z2 = z1 + t1;
		SectionPoint p2 = { .y = y2,.z = z2 };

		// point 3
		double y3 = y1 - (w1 - t2) / x1_;
		double z3 = z1;
		SectionPoint p3 = { .y = y3,.z = z3 };

		// point 4
		double y4 = y3;
		double z4 = z2;
		SectionPoint p4 = { .y = y4,.z = z4 };

		// create cell
		CellNodes cn = { .p1 = p1,.p2 = p2,.p3 = p3,.p4 = p4 };
		SectionCell* c = new SectionCell(cn, integrals_points_);
		s->cells.push_back(c);
	}

	// corner part
	{
		// point 1
		double y1 = t2 - yc;
		double z1 = -zc;
		SectionPoint p1 = { .y = y1,.z = z1 };

		// point 2
		double y2 = y1;
		double z2 = z1 + t1;
		SectionPoint p2 = { .y = y2,.z = z2 };

		// point 3
		double y3 = y1 - t2;
		double z3 = z1;
		SectionPoint p3 = { .y = y3,.z = z3 };

		// point 4
		double y4 = y3;
		double z4 = z2;
		SectionPoint p4 = { .y = y4,.z = z4 };

		// create cell
		CellNodes cn = { .p1 = p1,.p2 = p2,.p3 = p3,.p4 = p4 };
		SectionCell* c = new SectionCell(cn, integrals_points_);
		s->cells.push_back(c);
	}

	// horizontal part
	for (int i = 0; i < x2_; i++)
	{
		// point 1
		double y1 = t2 - yc;
		double z1 = t1 - zc + (w2 - t1) / x2_ * i;
		SectionPoint p1 = { .y = y1,.z = z1 };

		// point 2
		double y2 = y1;
		double z2 = z1 + (w2 - t1) / x2_;
		SectionPoint p2 = { .y = y2,.z = z2 };

		// point 3
		double y3 = y1 - t2;
		double z3 = z1;
		SectionPoint p3 = { .y = y3,.z = z3 };

		// point 4
		double y4 = y3;
		double z4 = z2;
		SectionPoint p4 = { .y = y4,.z = z4 };

		// create cell
		CellNodes cn = { .p1 = p1,.p2 = p2,.p3 = p3,.p4 = p4 };
		SectionCell* c = new SectionCell(cn, integrals_points_);
		s->cells.push_back(c);
	}

	// initiate force container
	for (int i = 0; i < 5; i++)
	{
		s->force.push_back(0);
	}
}
