// Standalone C++ version of the monoparametric tiling
// Author: Guillaume Iooss

#ifndef __LINALG_H
#define __LINALG_H

#include <iostream>
#include <list>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define int64 long long int

// To disable the "assert"
// #define NDEBUG


using namespace std;


/* ~~~~~~ Matrix operations ~~~~~~ */

// Free a matrix
void freeMatrix(int64** mat, int nRow);

// Copy a matrix
int64** copyMatrix(int64** mat, int nRow, int nCol);

// Pretty-printer
void printMatrix(int64** mat, int nRow, int nCol);

// Get the submatrix of a matrix composed of the columns going from j1 to j2 (inclusive on j1, exclusive on j2)
int64** submatrixColumn(int64** A, int nRowA, int nColA, int j1, int j2);

// Swap rows i1 and i2 of a matrix
void swapRows(int64** mat, int nRow, int nCol, int i1, int i2);


// Perform the row operation on mat: L_i1 <- L_i1 + scalar*L_i2
void rowAddition(int64** mat, int nRow, int nCol, int i1, int i2, int scalar);

// Matrix multiplication mat1*mat2
int64** matrixMultiplication(int64** mat1, int nRow1, int nCol1, int64** mat2, int nRow2, int nCol2);

// Matrix-vector multiplication mat1*vect2
int64* matrixVectMultiplication(int64** mat1, int nRow1, int nCol1, int64* vect2, int nRow2);

// Vector-matrix multiplication vect1*mat2
int64* vectMatrixMultiplication(int64* vect1, int nCol1, int64** mat2, int nRow2, int nCol2);

// Inversion of a unimodular matrix
int64** inverseMatUnimod(int64** unimodMatinv, int nRow, int nCol);


/* ~~~~~~ Polyhedron ~~~~~~ */

// Datastructure for polyhedron - similar to the Polylib format
// Ineq? |  Ind  |  Param  | Const 
// --------------------------------
//  0/1  |  ...  |   ...   |  ...  
// 0 = Equality, 1 = Inequality
// Number of rows in the matrix = nConstr
// Number of cols in the matrix = 2 + nParam + nInd
struct polyhedronMPP {
	int64** polyScalar;
	
	int nConstr;        /* Number of constraints */
	
	int nInd;           /* Number of indices */
	int nParam;         /* Number of parameters */
	
	// Modulo constraints
	// Format: like Polylib, except that the first column is the modulo factor
	// Example: (i -N -1) % 2 = 0 ===> 2 | 1 | -1 | -1
	int64** modConstr;
	int nConstrMod;
};

// Constructor
polyhedronMPP* buildPolyhedron(int64** mat, int nConstr, int nInd, int nParam);


// Constructor from different matrices
polyhedronMPP* reformPoly(int64* ineqEqPart, int64** paramPart, int64** linPart, int64* constPart, int nConstr, int nParam, int nInd);


// Destructor
void freePolyhedron(polyhedronMPP* poly);

// Getters (all in one). The arrays must be preallocated outside of this function
void extractPoly(polyhedronMPP* poly, int64* ineqEqPart, int64** linPart, int64** paramPart, int64* constPart);

// Transform the equalities of the polyhedron (A=0) into 2 inequalities (A>=0 and -A>=0)
polyhedronMPP* eliminateEqualities(polyhedronMPP* poly);

// Cartesian product of two polyhedra
polyhedronMPP* cartesianProduct(polyhedronMPP* poly1, polyhedronMPP* poly2);

// Pretty-printer
void printPolyhedronMPP(polyhedronMPP* poly);


// Datastructure for affine functions
//  Ind  | Param | Const 
// ----------------------
//  ...  |  ...  |  ...  
struct affFuncMPP {
	int64** affFuncScalar;
	
	int dimOut;			/* Dimension of the output */
	
	int nInd;           /* Number of indices */
	int nParam;         /* Number of parameters */	
};

// Constructor
affFuncMPP* buildAffineFunction(int64** mat, int dimOut, int nInd, int nParam);

// Destructor
void freeAffineFunction(affFuncMPP* affFunc);

// Pretty-printer
void printAffFuncMPP(affFuncMPP* func);

// Getters (all in one)
void extractAffFunc(affFuncMPP* affFunc, int64** linPart, int64** paramPart, int64* constPart);

// Given a polyhedron and a matrix, compute the resulting polyhedron after a change of basis
polyhedronMPP* changeOfBasis(polyhedronMPP* poly, affFuncMPP* affFunc);

#endif

