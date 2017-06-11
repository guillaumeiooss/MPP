#include "linAlg.h"

// Helper function
void freeMatrix(int64** mat, int nRow) {
	for (int i=0; i<nRow; i++)
		free(mat[i]);
	free(mat);
	return;
}



// Copy matrix
int64** copyMatrix(int64 **mat, int nRow, int nCol) {
	assert(nRow>=0);
	assert(nCol>=0);
	
	int64** copyMat = (int64**) malloc(nRow * sizeof(int64));
	for (int i=0; i<nRow; i++) {
		copyMat[i] = (int64*) malloc (nCol * sizeof(int64));
		for (int j=0; j<nCol; j++)
			copyMat[i][j] = mat[i][j];
	}
	return copyMat;
}


// Pretty-printer
void printMatrix(int64** mat, int nRow, int nCol) {
	for (int i=0; i<nRow; i++) {
		cout << "[ ";
		for (int j=0; j<nCol; j++) {
			cout << mat[i][j] << " ";
		}
		cout << "]" << endl;
	}
	cout << endl;
	
	return;
}

void printVector(int64* vect, int nElem) {
	cout << "[ ";
	for (int j=0; j<nElem; j++) {
		cout << vect[j] << " ";
	}
	cout << "]" << endl;
	return;
}



// Extract a submatrix made of columns
int64** submatrixColumn(int64 **A, int nRowA, int nColA, int j1, int j2) {
    assert(0<=j1 && j1<j2 && j2<=nColA);
    
	int64 **B = (int64**) malloc(nRowA * sizeof(int64*));
	for (int i=0; i<nRowA; i++) {
		B[i] = (int64*) malloc((j2-j1) * sizeof(int64));
		for (int j=j1; j<j2; j++)
			B[i][j-j1] = A[i][j];
	}
	return B;
}



// Swap 2 rows in a matrix
void swapRows(int64 **mat, int nRow, int nCol, int i1, int i2) {
	assert(0<=i1 && i1<nRow);
	assert(0<=i2 && i2<nRow);
	
	int64 *rowTemp = (int64*) malloc(nCol * sizeof(int64));
	
	// rowTemp <- Row i1
	for (int j=0; j<nCol; j++)
		rowTemp[j] = mat[i1][j];
	
	// Row i1 <- Row i2
	for (int j=0; j<nCol; j++)
		mat[i1][j] = mat[i2][j];
	
	// Row i2 <- rowTemp
	for (int j=0; j<nCol; j++)
		mat[i2][j] = rowTemp[j];
	
	free(rowTemp);
	
	return;
}



// Row addition: L_i1 <- L_i1 + scalar*L_i2
void rowAddition(int64 **mat, int nRow, int nCol, int i1, int i2, int scalar) {
	assert(0<=i1 && i1<nRow);
	assert(0<=i2 && i2<nRow);
	
	for (int j=0; j<nCol; j++)
		mat[i1][j] += scalar * mat[i2][j];
	return;
}



// Matrix multiplication
int64** matrixMultiplication(int64 **mat1, int nRow1, int nCol1, int64 **mat2, int nRow2, int nCol2) {
	assert(nCol1==nRow2);
	
	int64 **retMat = (int64**) malloc(nRow1 * sizeof(int64*));
	
	for (int i=0; i<nRow1; i++) {
		retMat[i] = (int64*) malloc(nCol2 * sizeof(int64));
		
		for (int j=0; j<nCol2; j++) {
			retMat[i][j] = 0;
			for (int k=0; k<nCol1; k++)
				retMat[i][j] += mat1[i][k] * mat2[k][j];
		}
	}
	return retMat;
}


// Matrix-vector multiplication mat1*vect2
int64* matrixVectMultiplication(int64** mat1, int nRow1, int nCol1, int64* vect2, int nRow2) {
	assert(nCol1==nRow2);
	
	int64* retVect = (int64*) malloc(nRow1 * sizeof(int64));
	
	for (int i=0; i<nRow1; i++) {
		retVect[i] = 0;
		for (int j=0; j<nCol1; j++)
			retVect[i] += mat1[i][j] * vect2[j];
	}
	return retVect;
}


// Vector-matrix multiplication vect1*mat2
int64* vectMatrixMultiplication(int64* vect1, int nCol1, int64** mat2, int nRow2, int nCol2) {
	assert(nCol1==nRow2);
	
	int64* retVect = (int64*) malloc(nCol2 * sizeof(int64));
	
	for (int i=0; i<nCol2; i++) {
		retVect[i] = 0;
		for (int j=0; j<nCol1; j++)
			retVect[i] += vect1[j] * mat2[j][i];
	}
	return retVect;
}



// Inversion of a unimodular matrix
int64** inverseMatUnimod(int64 **unimodMatinv, int nRow, int nCol) {
	assert(nRow==nCol);
	
	// Copy the matrix into A
	int64** A = copyMatrix(unimodMatinv, nRow, nCol);
	
	if (nRow==1)
		return A;
	
	// First, deal with the left upper corner
	
	// Left unimodular inverse of the matrix
	int64** left = (int64**) malloc(nRow * sizeof(int64*));
	for (int i=0; i<nRow; i++) {
		left[i] = (int64*) malloc(nRow * sizeof(int64));
		for (int j=0; j<nRow; j++)
			left[i][j] = 0;
	}
	
	for (int i=0; i<nRow; i++) {
		long sc = 1;
		
		if (A[i][0]!=0) {
			sc = A[i][0] / abs(A[i][0]);
		}
		for (int j=0; j<nRow; j++) {
			A[i][j] *= sc;
			if (i==j)
				left[i][j] = sc;
		}
	}
	
	// Massage the column i
	for (int i=0; i<nRow; i++) {
		if (A[i][0] != 0) {
			swapRows(A, nRow, nRow, 0, i);
			swapRows(left, nRow, nRow, 0, i);
			break;
		}
	}
	
	for (int j=1; j<nRow; j++) {
		while (A[j][0] % A[0][0] !=0) {
			long scalar = (long) -A[j][0] / A[0][0];
			
			rowAddition(A, nRow, nRow, j, 0, scalar);
			rowAddition(left, nRow, nRow, j, 0, scalar);
			swapRows(A, nRow, nRow, 0, j);
			swapRows(left, nRow, nRow, 0, j);
		}
		long scalar = (long) -A[j][0]/A[0][0];
		rowAddition(A, nRow, nRow, j, 0, scalar);
		rowAddition(left, nRow, nRow, j, 0, scalar);
	}
	
	// Massage the row i
	int64** smaller = (int64**) malloc((nRow-1)*sizeof(int64*));
	for (int i=0; i<nRow-1; i++) {
		smaller[i] = (int64*) malloc((nRow-1)*sizeof(int64));
		for (int j=0; j<nRow-1; j++)
			smaller[i][j] = A[i+1][j+1];
	}
	
	// vect: 1*(nRow-1)
	int64** vect = (int64**) malloc(sizeof(int64*));
	vect[0] = (int64*) malloc((nRow-1)*sizeof(int64));
	for (int i=0; i<nRow-1; i++)
		vect[0][i] = A[0][i+1];
	
	// smaller: (nRow-1) * (nRow-1)
	smaller = inverseMatUnimod(smaller, nRow-1, nRow-1);
	vect = matrixMultiplication(vect, 1, nRow-1, smaller, nRow-1, nRow-1);
	
	// left2 : nRow * nRow
	int64** left2 = (int64**) malloc(nRow * sizeof(int64*));
	for (int i=0; i<nRow; i++) {
		left2[i] = (int64*) malloc(nRow * sizeof(int64));
		
		for (int j=0; j<nRow; j++) {
			if (i==0 && i!=j)
					left2[i][j] = -vect[0][j-1];
			else if (i==j)
				left2[i][i] = 1;
			else
				left2[i][j] = 0;
		}
	}
	
	// left3 : nRow * nRow
	int64** left3 = (int64**) malloc(nRow * sizeof(int64*));
	for (int i=0; i<nRow; i++)
		left3[i] = (int64*) malloc(nRow * sizeof(int64));
	
	left3[0][0] = 1;
	for (int i=1; i<nRow; i++) {
		left3[0][i] = 0;
		left3[i][0] = 0;
	}
	for (int i=1; i<nRow; i++)
		for (int j=1; j<nRow; j++)
			left3[i][j] = smaller[i-1][j-1];
	
	int64** temp = matrixMultiplication(left2, nRow, nRow, left, nRow, nRow);
	int64** result = matrixMultiplication(left3, nRow, nRow, temp, nRow, nRow);
	
	// Free temporary arrays
	freeMatrix(A,nRow);
	freeMatrix(left,nRow);
	freeMatrix(smaller,nRow-1);
	freeMatrix(vect,1);
	freeMatrix(left2,nRow);
	freeMatrix(left3,nRow);
	freeMatrix(temp,nRow);
	
	return result;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// Constructor
polyhedronMPP* buildPolyhedron(int64** mat, int nConstr, int nInd, int nParam) {
	polyhedronMPP *res = (polyhedronMPP*) malloc(sizeof(polyhedronMPP));
	
	res->polyScalar = mat;
	res->nConstr = nConstr;
	res->nInd = nInd;
	res->nParam = nParam;
	
	res->nConstrMod = 0;
	res->modConstr = (int64**) malloc(0*sizeof(int64*));
	return res;
}



void freePolyhedron(polyhedronMPP* poly) {
	freeMatrix(poly->polyScalar, poly->nConstr);
	freeMatrix(poly->modConstr, poly->nConstrMod);
	free(poly);
	return;
}



void printPolyhedronMPP(polyhedronMPP* poly) {
	for (int i=0; i<poly->nConstr; i++) {
		cout << "		(";
		for (int j=0; j<poly->nInd+poly->nParam+2; j++)
			cout << " " << poly->polyScalar[i][j];
		cout << " )" << endl;
	}
	cout << endl;
}



affFuncMPP* buildAffineFunction(int64** mat, int dimOut, int nInd, int nParam) {
	affFuncMPP *res = (affFuncMPP*) malloc(sizeof(affFuncMPP));
	
	res->affFuncScalar = mat;
	res->dimOut = dimOut;
	res->nInd = nInd;
	res->nParam = nParam;
	
	return res;
}



void freeAffineFunction(affFuncMPP* affFunc) {
	freeMatrix(affFunc->affFuncScalar, affFunc->dimOut);
	free(affFunc);
	return;
}



void printAffFuncMPP(affFuncMPP* func) {
	for (int i=0; i<func->dimOut; i++) {
		cout << "		(";
		for (int j=0; j<func->nInd+func->nParam+1; j++)
			cout << " " << func->affFuncScalar[i][j];
		cout << " )" << endl;
	}
	cout << endl;
}



void extractPoly(polyhedronMPP* poly, int64* ineqEqPart, int64** linPart, int64** paramPart, int64* constPart) {
	int nParam = poly->nParam;
	int nInd = poly->nInd;
	int nConstr = poly->nConstr;
	int64** matScalar = poly->polyScalar;
	
	for (int i=0; i<nConstr; i++) {
		ineqEqPart[i] = matScalar[i][0];
		for (int j=0; j<nInd; j++)
			linPart[i][j] = matScalar[i][1+j];
		for (int j=0; j<nParam; j++)
			paramPart[i][j] = matScalar[i][1+nInd+j];
		constPart[i] = matScalar[i][1+nInd+nParam];
	}
	return;
}



// Transform the equalities of the polyhedron (A=0) into 2 inequalities (A>=0 and -A>=0)
polyhedronMPP* eliminateEqualities(polyhedronMPP* poly) {
	int nColumns = 2+poly->nInd+poly->nParam;
	
	// Get the equality constraints
	list<int64*> lequalityConstr;
	
	for (int i=0; i<poly->nConstr; i++)
	if (poly->polyScalar[i][0]==0) {
		// If we have an equality: copy the constraint in equalityConst
		int64* eqConstr = (int64*) malloc(nColumns * sizeof(int64));
		for (int j=0; j<nColumns; j++)
			eqConstr[j] = poly->polyScalar[i][j];
		lequalityConstr.push_back(eqConstr);
	}
	
	int n_eq_constr = lequalityConstr.size();
	
	// Create the new matrix of constraint
	int nConstr_nEq = poly->nConstr + n_eq_constr;
	int64** nMatConstr = (int64**) malloc(nConstr_nEq * sizeof(int64*)) ;
	for (int i=0; i<poly->nConstr; i++) {
		nMatConstr[i] = (int64*) malloc(nColumns * sizeof(int64));
		
		// Modify the constraint A==0 into A>=0 (by changing the first column)
		nMatConstr[i][0] = 1;
		for (int j=1; j<nColumns; j++)
			nMatConstr[i][j] = poly->polyScalar[i][j];
	}
	for (int i=0; i<n_eq_constr; i++) {
		int shiftedI = i+poly->nConstr;
		nMatConstr[shiftedI] = (int64*) malloc(nColumns * sizeof(int64));
		
		int64* eqConstr = lequalityConstr.front(); // lequalityConstr.get(i)
		
		nMatConstr[shiftedI][0] = 1;		// Is an inequality
		for (int j=1; j<nColumns; j++)
			nMatConstr[shiftedI][j] = -eqConstr[j];
		
		// Free temporary structure
		lequalityConstr.pop_front();
		free(eqConstr);
	}
	
	return buildPolyhedron(nMatConstr, nConstr_nEq, poly->nInd, poly->nParam);
}


polyhedronMPP* cartesianProduct(polyhedronMPP* poly1, polyhedronMPP* poly2) {
	assert(poly1->nParam==poly2->nParam);
	
	int nInd1 = poly1->nInd;
	int nInd2 = poly2->nInd;
	int nParam = poly1->nParam;
	int nConstr1 = poly1->nConstr;
	int nConstr2 = poly2->nConstr;
	
	// Extracting the information from both polyhedra
	int64* ineqEqPart1 = (int64*) malloc(nConstr1 * sizeof(int64));
	int64** linPart1 = (int64**) malloc(nConstr1 * sizeof(int64*));
	for (int i=0; i<nConstr1; i++)
		linPart1[i] = (int64*) malloc(nInd1 * sizeof(int64));
	int64** paramPart1 = (int64**) malloc(nConstr1 * sizeof(int64*));
	for (int i=0; i<nConstr1; i++)
		paramPart1[i] = (int64*) malloc(nParam * sizeof(int64));
	int64* constPart1 = (int64*) malloc(nConstr1 * sizeof(int64));
	extractPoly(poly1, ineqEqPart1, linPart1, paramPart1, constPart1);
	
	int64* ineqEqPart2 = (int64*) malloc(nConstr2 * sizeof(int64));
	int64** linPart2 = (int64**) malloc(nConstr2 * sizeof(int64*));
	for (int i=0; i<nConstr2; i++)
		linPart2[i] = (int64*) malloc(nInd2 * sizeof(int64));
	int64** paramPart2 = (int64**) malloc(nConstr2 * sizeof(int64*));
	for (int i=0; i<nConstr2; i++)
		paramPart2[i] = (int64*) malloc(nParam * sizeof(int64));
	int64* constPart2 = (int64*) malloc(nConstr2 * sizeof(int64));
	extractPoly(poly2, ineqEqPart2, linPart2, paramPart2, constPart2);
	
	
	// Building the cartesian product of both polyhedra
	int nConstrProd = nConstr1 + nConstr2;
	int nIndProd = nInd1 + nInd2;
	int nColumnProd = 2+nParam+nIndProd;
	
	int64** matProd = (int64**) malloc(nConstrProd * sizeof(int64*));
	for (int i=0; i<nConstrProd; i++)
		matProd[i] = (int64*) malloc(nColumnProd * sizeof(int64));
	
	for (int i=0; i<nConstr1; i++) {
		matProd[i][0] = ineqEqPart1[i];
		for (int j=0; j<nInd1; j++)
			matProd[i][1+j] = linPart1[i][j];
		for (int j=0; j<nParam; j++)
			matProd[i][1+nIndProd+j] = paramPart1[i][j];
		matProd[i][nColumnProd-1] = constPart1[i];
	}
	for (int i=0; i<nConstr2; i++) {
		matProd[nConstr1+i][0] = ineqEqPart2[i];
		for (int j=0; j<nInd2; j++)
			matProd[nConstr1+i][1+nInd1+j] = linPart2[i][j];
		for (int j=0; j<nParam; j++)
			matProd[nConstr1+i][1+nIndProd+j] = paramPart2[i][j];
		matProd[nConstr1+i][nColumnProd-1] = constPart2[i];
	}
	polyhedronMPP* polyProd = buildPolyhedron(matProd, nConstrProd, nIndProd, nParam);
	
	return polyProd;
}


void extractAffFunc(affFuncMPP* affFunc, int64** linPart, int64** paramPart, int64* constPart) {
	int nParam = affFunc->nParam;
	int nInd = affFunc->nInd;
	int dimOut = affFunc->dimOut;
	int64** matScalar = affFunc->affFuncScalar; 
	
	for (int i=0; i<dimOut; i++) {
		for (int j=0; j<nInd; j++)
			linPart[i][j] = matScalar[i][j];
		for (int j=0; j<nParam; j++)
			paramPart[i][j] = matScalar[i][nInd+j];
		constPart[i] = matScalar[i][nInd+nParam];
	}
	
	return;
}



polyhedronMPP* reformPoly(int64* ineqEqPart, int64** paramPart, int64** linPart, int64* constPart, int nConstr, int nParam, int nInd) {
	int64** matRet = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++) {
		matRet[i] = (int64*) malloc((2+nInd+nParam)*sizeof(int64));
		
		matRet[i][0] = ineqEqPart[i];
		for (int j=0; j<nInd; j++)
			matRet[i][j+1] = linPart[i][j];
		for (int j=0; j<nParam; j++)
			matRet[i][1+nInd+j] = paramPart[i][j];
		matRet[i][nInd+nParam+1] = constPart[i];
	}
	polyhedronMPP* polyRet = buildPolyhedron(matRet, nConstr, nInd, nParam);
	
	return polyRet;
}



// Given a polyhedron and a matrix, compute the resulting polyhedron after a change of basis
polyhedronMPP* changeOfBasis(polyhedronMPP* poly, affFuncMPP* affFunc) {
	assert(affFunc->nInd == affFunc->dimOut);
	
	assert(affFunc->nParam == poly->nParam);
	assert(affFunc->nInd == poly->nInd);
	
	// Extraction of constraints
	int64* ineqEqPart = (int64*) malloc(poly->nConstr * sizeof(int64));
	int64** linPart = (int64**) malloc(poly->nConstr * sizeof(int64));
	for (int i=0; i<poly->nConstr; i++)
		linPart[i] = (int64*) malloc(poly->nInd * sizeof(int64));
	int64** paramPart = (int64**) malloc(poly->nConstr * sizeof(int64));
	for (int i=0; i<poly->nConstr; i++)
		paramPart[i] = (int64*) malloc(poly->nParam * sizeof(int64));
	int64* constPart = (int64*) malloc(poly->nConstr * sizeof(int64));
	extractPoly(poly, ineqEqPart, linPart, paramPart, constPart);
	
	// Extraction of mat
	int64** linPartMatAff = (int64**) malloc(affFunc->dimOut * sizeof(int64*));
	for (int i=0; i<affFunc->dimOut; i++)
		linPartMatAff[i] = (int64*) malloc(affFunc->nInd * sizeof(int64));
	int64** paramPartMatAff = (int64**) malloc(affFunc->dimOut * sizeof(int64*));
	for (int i=0; i<affFunc->dimOut; i++)
		paramPartMatAff[i] = (int64*) malloc(affFunc->nParam * sizeof(int64));
	int64* constMatAff = (int64*) malloc(affFunc->dimOut * sizeof(int64));
	extractAffFunc(affFunc, linPartMatAff, paramPartMatAff, constMatAff);
	
	int64** linPartMatInv = inverseMatUnimod(linPartMatAff, affFunc->dimOut, affFunc->nInd);
	
	// Computes nLinPartMat = linPart*linPartMatInv
	int64** nLinPartMat = matrixMultiplication(linPart, poly->nConstr, poly->nInd, linPartMatInv, affFunc->nInd, affFunc->nInd);
	
	// Computes nParamPartMat = paramPart - linPart * linPartMatInv * paramPartMatAff
	int64** tempMatLin = matrixMultiplication(linPartMatInv, affFunc->nInd, affFunc->nInd, paramPartMatAff, affFunc->dimOut, affFunc->nParam);
	int64** nParamPartMat = matrixMultiplication(linPart, poly->nConstr, poly->nInd, tempMatLin, affFunc->nInd, affFunc->nParam);
	for (int i=0; i<poly->nConstr; i++)
		for (int j=0; j<affFunc->nParam; j++)
			nParamPartMat[i][j] = paramPart[i][j] - nParamPartMat[i][j];
	
	// Computes nConstPartMat = constPart - linPart * linPartMatInv * constMatAff
	int64* nConstPartMat = (int64*) malloc(poly->nConstr*sizeof(int64));
	
	int64** tempMatConst = matrixMultiplication(linPart, poly->nConstr, poly->nInd, linPartMatInv, affFunc->nInd, affFunc->dimOut);
	for (int i=0; i<poly->nConstr; i++) {
		nConstPartMat[i] = constPart[i];
		for (int k=0; k<affFunc->dimOut; k++)
			nConstPartMat[i] -= tempMatConst[i][k]*constMatAff[k];
	}
	
	// Reform the polyhedron by putting back the matrices in place
	polyhedronMPP* polyRet = reformPoly(ineqEqPart, nParamPartMat, nLinPartMat, nConstPartMat,
		poly->nConstr, poly->nParam, poly->nInd);
	
	
	// Free temporary arrays
	free(ineqEqPart);
	freeMatrix(linPart,poly->nConstr);
	freeMatrix(paramPart,poly->nConstr);
	free(constPart);
	
	freeMatrix(linPartMatAff,affFunc->dimOut);
	freeMatrix(paramPartMatAff,affFunc->dimOut);
	free(constMatAff);
	freeMatrix(linPartMatInv,affFunc->dimOut);
	
	freeMatrix(nLinPartMat,poly->nConstr);
	
	freeMatrix(tempMatLin,affFunc->nInd);
	freeMatrix(nParamPartMat,poly->nConstr);
	
	free(nConstPartMat);
	freeMatrix(tempMatConst,poly->nConstr);
	
	return polyRet;
}


polyhedronMPP* rectangularShape(int64* sizes, int nDim) {
	int nConstr = 2*nDim;
	int nCol_mat = 2+nDim+1;
	
	int64** matPolyRet = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		matPolyRet[i] = (int64*) malloc(nCol_mat * sizeof(int64));
	for (int i=0; i<nConstr; i++)
		for (int j=0; j<nCol_mat; j++)
			matPolyRet[i][j] = 0;
	
	// (First lines)
	for (int i=0; i<nDim; i++) {
		matPolyRet[i][0] = 1;					// Ineq
		matPolyRet[i][1+i] = 1;					// Lin
	}
	
	// (Second lines)
	for (int i=0; i<nDim; i++) {
		matPolyRet[nDim+i][0] = 1;				// Ineq
		matPolyRet[nDim+i][1+i] = -1;			// Lin
		matPolyRet[nDim+i][1+nDim] = sizes[i];	// Param "b"
		matPolyRet[nDim+i][1+nDim+1] = -1;		// Const
	}
	
	polyhedronMPP* polyRet = buildPolyhedron(matPolyRet, nConstr, nDim, 1);
	return polyRet;
}


polyhedronMPP* parallelogramShape(int64** hyperplanes, int64* sizes, int nDim) {
	int nConstr = 2*nDim;
	int nCol_mat = 2+nDim+1;
	
	int64** matPolyRet = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		matPolyRet[i] = (int64*) malloc(nCol_mat * sizeof(int64));
	for (int i=0; i<nConstr; i++)
		for (int j=0; j<nCol_mat; j++)
			matPolyRet[i][j] = 0;
	
	// (First lines)
	for (int i=0; i<nDim; i++) {
		matPolyRet[i][0] = 1;									// Ineq
		for (int j=0; j<nDim; j++)
			matPolyRet[i][1+j] = hyperplanes[j][i];				// Lin
	}
	
	// (Second lines)
	for (int i=0; i<nDim; i++) {
		matPolyRet[nDim+i][0] = 1;								// Ineq
		for (int j=0; j<nDim; j++)
			matPolyRet[nDim+i][1+j] = -hyperplanes[j][i];		// Lin
		matPolyRet[nDim+i][1+nDim] = sizes[i];					// Param "b"
		matPolyRet[nDim+i][1+nDim+1] = -1;						// Const
	}
	
	polyhedronMPP* polyRet = buildPolyhedron(matPolyRet, nConstr, nDim, 1);
	return polyRet;
}




