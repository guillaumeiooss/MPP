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

void printVector(int64* vect, int nElem);

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

// Dot-product vect1*vect2
int64 dotProduct(int64* vect1, int nElem1, int64* vect2, int nElem2);

// Inversion of a unimodular matrix
int64** inverseMatUnimod(int64** unimodMatinv, int nRow, int nCol);



/* ~~~~~~ Rational matrix ~~~~~~ */
struct rational64 {
	int64 num;
	int64 den;				// Must always be non-negative (den > 0)
};

int64 gcd(int64 a, int64 b);
int64 ppcm(int64 a, int64 b);
int64 gcd_array(int64* elems, int nelem);
int64 ppcm_array(int64* elems, int nelem);


// Basic operations on the rational64
rational64 toRational(int64 x);
rational64 simplify(rational64 rat);
rational64 multiplyRational(rational64 a, rational64 b);
rational64 addRational(rational64 a, rational64 b);
rational64 subRational(rational64 a, rational64 b);
rational64 oppositeRational(rational64 a);
rational64 invertRational(rational64 rat);


// Transform an integral vector into a rational vector
rational64* toRationalVector(int64* vect, int nElem);

// Pretty-printer
void printMatrix(rational64** mat, int nRow, int nCol);

void printVector(rational64* vect, int nElem);

// Free a matrix
void freeMatrix(rational64** mat, int nRow);

// Transform an integral matrix into a rational matrix
rational64** toRationalMatrix(int64** mat, int nRow, int nCol);

// Transform a rational matrix into an integral matrix
int64** toIntegralMatrix(rational64** ratmat, int nRow, int nCol);

// Dot product between two rational vectors
rational64 dotProduct(rational64* vect1, rational64* vect2, int nElem);

// Matrix-scalar product
rational64** matScalProduct(rational64** mat, int nRow, int nCol, rational64 elem);

// Matrix-vector product
rational64* matVectProduct(rational64** mat, int nRow, int nCol, rational64* vect, int nElem);

// Matrix multiplication for rational matrices
rational64** matrixMultiplication(rational64** mat1, int nRow1, int nCol1, rational64** mat2, int nRow2, int nCol2);

// Inplace swapping of matrix rows
void swapRowsMatrix(rational64** DG, int nRow, int nCol, int r1, int r2);

// L_r1 <- L_r1 + coeff*L_r2
void rowAddition(rational64** DG, int nRow, int nCol, int r1, int r2, rational64 coeff);

// Transpose a matrix
rational64** transpose(rational64** mat, int nRow, int nCol);


// Get the inverse and the determinant of a matrix by the Gauss pivot method
//		A is a square matrix, invA is a preallocated square matrix
//		The determinant is returned. The inverse is stored in invA
int64 inverseDet(int64** A, rational64** invA, int nRow);



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

// Build a rectangle [|\vec{0}; b*\vec{sizes} [|
polyhedronMPP* rectangularShape(int64* sizes, int nDim);

// Build the lattice of tile origins, corresponding to the tiling of tile [|\vec{0}; b*\vec{sizes} [|
int64** rectangularOriginLattice(int64* sizes, int nDim);

// Build a parallelogram, using hyperplanes (which are the columns of the square matrix "hyperplanes")
// Contraints: { \vec{i} | 0 <= \vec{v_k}.\vec{i} <= sizes[k].b -1 }  sizes[k].b where \vec{v_k} is a column vector
polyhedronMPP* parallelogramShape(int64** hyperplanes, int64* sizes, int nDim);

// Build the lattice of tile origins, using hyperplanes (which are the columns of the square matrix "hyperplanes")
// Compute the (rational) inverse of "hyperplanes", and multiply it by the sizes
int64** parallelogramOriginLattice(int64** hyperplanes, int64* sizes, int nDim);


/* ~~~~~~ Affine function ~~~~~~ */

// Datastructure for affine functions
//  Ind  | Param | Const 
// ----------------------
//  ...  |  ...  |  ...  
struct affFuncMPP {
	int64** affFuncScalar;  /* Coefficients of the affine function */
	int64* divs;            /* Divisor of the whole row. If integral, equal to "1". Should always be >=1 */
	
	int dimOut;			/* Dimension of the output */
	
	int nInd;           /* Number of indices */
	int nParam;         /* Number of parameters */	
};

// Constructor
affFuncMPP* buildAffineFunction(int64** mat, int dimOut, int nInd, int nParam);

affFuncMPP* buildRationalAffineFunction(int64** mat, int64* div, int dimOut, int nInd, int nParam);

// Destructor
void freeAffineFunction(affFuncMPP* affFunc);

// Pretty-printer
void printAffFuncMPP(affFuncMPP* func);

// Simplify the divs with the coefficients, if possible
void simplifyAffFuncMPP(affFuncMPP* func);

// Getters (all in one)
void extractAffFunc(affFuncMPP* affFunc, int64** linPart, int64** paramPart, int64* constPart);

void extractRationalAffFunc(affFuncMPP* affFunc, int64** linPart, int64** paramPart, int64* constPart, int64* div);

// Given a polyhedron and a matrix, compute the resulting polyhedron after a change of basis
polyhedronMPP* changeOfBasis(polyhedronMPP* poly, affFuncMPP* affFunc);


// Get the isl relation representation of an affine function
int64** getRelationRepresentation(affFuncMPP* func);

#endif

