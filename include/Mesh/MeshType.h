//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//+++ Date   : 2021.08.21  add abaqus meshtype names

#ifndef ASFEM_MESHTYPE_H
#define ASFEM_MESHTYPE_H

enum class MeshType{
    NULLTYPE,
    POINT, //MASS	Point mass
    EDGE2,
    EDGE3,
    EDGE4,
    EDGE5,
    TRI3,
    TRI6,
    QUAD4,
    QUAD8,
    QUAD9,
    TET4,
    TET10,
    HEX8,
    HEX20,
    HEX27,

	TRUSS2,//LINE2 ? EDGE
	TRUSS3,
	BEAM2,
	BEAM3,
	PIPE2,
	PIPE3,

	//element index reference from abaqus
	B21,	//  2-node linear beam in a plane		||Active degrees of freedom 1, 2, 6
	B22,	//  3-node quadratic beam in a plane	||Active degrees of freedom 1, 2, 6
	B31,	//	2-node linear beam in space			||Active degrees of freedom 1, 2, 3, 4, 5, 6
	B32,	//	3-node quadratic beam in space		||Active degrees of freedom 1, 2, 3, 4, 5, 6
	C3D4,	//	4-node linear tetrahedron			||Active degrees of freedom 1, 2, 3
	C3D8,	//	8-node linear brick					||Active degrees of freedom 1, 2, 3
	C3D10M,	//	10-node modified second - order tetrahedron		||Active degrees of freedom 1, 2, 3
	CPE3,	//	3-node linear plane strain triangle				||Active degrees of freedom 1, 2
	CPE4,	//	4-node bilinear plane strain quadrilateral												||Active degrees of freedom 1, 2
	CPE4R,	//	4-node bilinear plane strain quadrilateral, reduced integration, hourglass control		||Active degrees of freedom 1, 2
	CPE6M,	//	6-node modified second - order plane strain triangle		||Active degrees of freedom 1, 2
	CPS3,	//	3-node linear plane stress triangle							||Active degrees of freedom 1, 2
	CPS4,	//	4-node bilinear plane stress quadrilateral					||Active degrees of freedom 1, 2
	CPS4R,	//	4-node bilinear plane stress quadrilateral, reduced integration, hourglass control		||Active degrees of freedom 1, 2
	CPS6M,	//	6-node modified second - order plane stress triangle		||Active degrees of freedom 1, 2
	M3D3,	//	3-node triangular membrane						||Active degrees of freedom 1, 2, 3
	M3D4,	//	4-node quadrilateral membrane					||Active degrees of freedom 1, 2, 3
	//MASS	Point mass		//	||Active degrees of freedom 1, 2, 3
	PIPE21,	//	2-node linear pipe in a plane		||Active degrees of freedom 1, 2, 6
	PIPE31,	//	2-node linear pipe in space			||Active degrees of freedom 1, 2, 3, 4, 5, 6
	S3,		//  3-node triangular shell, finite membrane strains			||Active degrees of freedom 1, 2, 3, 4, 5, 6
	S3R,	//  3-node triangular shell, finite membrane strains			||Active degrees of freedom 1, 2, 3, 4, 5, 6
	S4,		//	4-node general - purpose shell, finite membrane strains		||Active degrees of freedom 1, 2, 3, 4, 5, 6
	S4R,	//	4-node general - purpose shell, finite membrane strains		||Active degrees of freedom 1, 2, 3, 4, 5, 6
	SC6R,	//	6-node triangular in - plane continuum shell wedge, general - purpose continuum shell, finite membrane strains.	||Active degrees of freedom 1, 2, 3
	SFM3D3,	//	3-node triangular surface element							||Active degrees of freedom 1, 2, 3
	SFM3D4,		//4-node quadrilateral surface element						||Active degrees of freedom 1, 2, 3
	SFM3D4R,	//4-node quadrilateral surface element, reduced integration	||Active degrees of freedom 1, 2, 3
	//SPRINGA	Axial spring between two nodes ||Active degrees of freedom 1, 2, 3  //Axial spring between two nodes, whose line of action is the line joining the two nodes. This line of action may rotate in large-displacement
	T2D2,	//	2-node linear 2D truss		||Active degrees of freedom 1, 2
	T3D2	//	2-node linear 3D truss		||Active degrees of freedom 1, 2, 3

};

#endif // ASFEM_MESHTYPE_H