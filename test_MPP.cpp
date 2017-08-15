// Test of the CART functions
// Copyright Guillaume Iooss, 2014, All right reserved.

#include "MPP_rect.h"
#include "lexmin.h"
#include "MPP_gen.h"

void printoutDomain(list<list<polyhedronMPP*> > llDom) {
	cout << "Intersection between the " << llDom.size() << " following union(s) of polyhedra:" << endl;
	for (list<list<polyhedronMPP*> >::iterator lDom=llDom.begin(); lDom!=llDom.end(); lDom++) {
		cout << "	[[[" << endl;
		for (list<polyhedronMPP*>::iterator poly=lDom->begin(); poly!=lDom->end(); poly++)
			printPolyhedronMPP(*poly);
		cout << "]]]" << endl;
	}
	return;
}


// Example 1: { i,j | N-i-j-1>=0 } with tiles of size b*b
void test_MPP_Poly_Ex1() {
	int nInd = 2;
	int nParam = 1;
	int nConstr = 1;
	
	int64** mat = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		mat[i] = (int64*) malloc((2+nInd+nParam)* sizeof(int64));
	mat[0][0] = 1;
	mat[0][1] = -1;
	mat[0][2] = -1;
	mat[0][3] = 1;
	mat[0][4] = -1;
	polyhedronMPP *polyScalar = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 1;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	
	list<list<polyhedronMPP*> > resultDom = getRectangularTiledDomain(polyScalar, scale, opt);
	printoutDomain(resultDom);
	// RESULT: union of 3 polyhedra:
	//(0/1 | alpha beta | ii jj | Nbl Nloc | b | cnst)    // 0 = equality / 1 = inequality
	//( 0 | -1 -1 |  0  0 |  1 0 | 0 | -1 )
	//( 1 |  0  0 | -1 -1 |  0 1 | 1 | -1 )
	//( 1 |  0  0 |  1  0 |  0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 |  0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 |  0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 |  0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 |  0 0 | 1 | -3 )
	//		 => {alpha,beta, ii,jj | alpha+beta=Nbl-1 && Nloc+b-1>=ii+jj && 0<=(ii,jj)<b && b>=3}
	// and
	//( 0 | -1 -1 |  0  0 | 1 0 | 0 |  0 )
	//( 1 |  0  0 | -1 -1 | 0 1 | 0 | -1 )
	//( 1 |  0  0 |  1  0 | 0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 | 0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 | 0 0 | 1 | -3 )
	//		=> {alpha,beta, ii,jj | alpha+beta=Nbl && Nloc-1>=ii+jj && 0<=(ii,jj)<b && b>=3}
	// and
	//( 1 | -1 -1 |  0  0 | 1 0 | 0 | -2 )
	//( 1 |  0  0 |  1  0 | 0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 | 0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 | 0 0 | 1 | -3 )
	//		=> {alpha,beta, ii,jj | alpha+beta<=Nbl-2 && 0<=(ii,jj)<b && b>=3}

}


// Example 2: {i,j | i-j>=0 && 2j-i>=0 && N-1-j>=0 } with tiles of size b*b
void test_MPP_Poly_Ex2() {
	int nInd = 2;
	int nParam = 1;
	int nConstr = 3;
	
	int64** mat = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		mat[i] = (int64*) malloc((2+nInd+nParam)* sizeof(int64));
	mat[0][0] = 1;
	mat[0][1] = 1;
	mat[0][2] = -1;
	mat[0][3] = 0;
	mat[0][4] = 0;
	
	mat[1][0] = 1;
	mat[1][1] = -1;
	mat[1][2] = 2;
	mat[1][3] = 0;
	mat[1][4] = 0;
	
	mat[2][0] = 1;
	mat[2][1] = 0;
	mat[2][2] = -1;
	mat[2][3] = 1;
	mat[2][4] = -1;
	
	polyhedronMPP *polyScalar = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 1;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	
	list<list<polyhedronMPP*> > resultDom = getRectangularTiledDomain(polyScalar, scale, opt);
	printoutDomain(resultDom);
	// RESULT: intersection between 3 unions of polyhedra (each union corresponding to a given constraint)
	// (we omit systematically the last 5 lines of constaints, which are 0<=(ii,jj)<b && b>=3)
	//
	// * First union:
	//(0/1 | alpha beta | ii jj | Nbl Nloc | b | cnst)    // 0 = equality / 1 = inequality
	//( 0 | 1 -1 | 0  0 | 0 0 |  0 |  0 )
	//( 1 | 0  0 | 1 -1 | 0 0 |  0 |  0 )
	// and
	//( 1 | 1 -1 | 0  0 | 0 0 |  0 | -1 )
	//=> {alpha,beta,ii,jj | alpha=beta && ii>=jj && 0<=ii,jj<b && b>=3 }
	//				U {i,j,alpha,beta | alpha>=beta+1 && 0<=ii,jj<b && b>=3 }
	//
	// * Second union:
	//(0/1 | alpha beta | ii jj | Nbl Nloc | b | cnst)    // 0 = equality / 1 = inequality
	//( 0 | -1 2 |  0 0 | 0 0 |  0 |  0 )
	//( 1 |  0 0 | -1 2 | 0 0 |  0 |  0 )
	// and
	//( 0 | -1 2 |  0 0 | 0 0 |  0 |  1 )
	//( 1 |  0 0 | -1 2 | 0 0 | -1 |  0 )
	// and
	//( 1 | -1 2 |  0 0 | 0 0 |  0 | -1 )
	//=> {alpha,beta,ii,jj | 2.beta-alpha=0 && 2.jj-ii>=0 && 0<=ii,jj<b }
	//				U {alpha,beta,ii,jj | 2.beta-alpha=-1 && 2.jj-ii>=b && 0<=ii,jj<b }
	//				U {alpha,beta,ii,jj | 2.beta-alpha>=1 && 0<=ii,jj<b }
	//
	// * Third union:
	//( 0 | 0 -1 | 0  0 | 1 0 |  0 |  0 )
	//( 1 | 0  0 | 0 -1 | 0 1 |  0 | -1 )
	// and
	//( 1 | 0 -1 | 0  0 | 1 0 |  0 | -1 )
	//=> {alpha,beta,ii,jj | Nbl=beta && Nloc>=jj+1 && 0<=ii,jj<b }
	//				U {i,j,alpha,beta | Nbl-1>=beta && 0<=ii,jj<b }
	return;
}


/* ------------------------------------------ */

void printoutFunction(map<polyhedronMPP*, affFuncMPP*> mFunc) {
	cout << mFunc.size() << " branches:" << endl;
	for (map<polyhedronMPP*, affFuncMPP*>::const_iterator itMap=mFunc.begin(); itMap!=mFunc.end(); itMap++) {
		cout << "IF" << endl;
		printPolyhedronMPP(itMap->first);
		cout << "	THEN" << endl;
		printAffFuncMPP(itMap->second);
	}
}


// Example 1: (i,j -> i-1, j-1)  (Jacobi1D dependence)
void test_MPP_Func_Ex1() {
	int nInd = 2;
	int nParam = 0;
	int dimOut = 2;
	
	int64** mat = (int64**) malloc(dimOut * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		mat[i] = (int64*) malloc((nInd+nParam+1)* sizeof(int64));
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = -1;
	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = -1;
	affFuncMPP *affScalar = buildAffineFunction(mat, dimOut, nInd, nParam);
	
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 1;
	
	int64* scaleIm = (int64*) malloc(dimOut*sizeof(int64));
	scaleIm[0] = 1; scaleIm[1] = 1;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	
	
	map<polyhedronMPP*, affFuncMPP*> resultFunc = getRectangularTiledFunction(affScalar, scale, scaleIm, opt);
	printoutFunction(resultFunc);
	// RESULT: 4 branches
	// We ignore the last 5 constraints (0<=(ii,jj)<b && b>=3) of the conditions
	// 1) IF
	//(0/1 | alpha beta | ii jj | b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  0 | 1 | -1 )
	//( 1 | 0 0 |  0  1 | 1 | -1 )
	//( 1 | 0 0 | -1  0 | 0 |  0 )
	//( 1 | 0 0 |  0 -1 | 0 |  0 )
	// THEN
	//(alpha beta | ii jj | b | cnst)
	//( 1 0 | 0 0 | 0 | -1 )
	//( 0 1 | 0 0 | 0 | -1 )
	//( 0 0 | 1 0 | 1 | -1 )
	//( 0 0 | 0 1 | 1 | -1 )
	//	=> IF (ii==0 && jj==0 && b>=3) THEN (alpha-1, beta-1, ii+b-1, jj+b-1)
	//
	// 2) IF
	//(0/1 | alpha beta | ii jj | b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  0 | 0 | -1 )
	//( 1 | 0 0 |  0  1 | 1 | -1 )
	//( 1 | 0 0 | -1  0 | 1 |  0 )
	//( 1 | 0 0 |  0 -1 | 0 |  0 )
	// THEN
	//(alpha beta | ii jj | b | cnst)
	//( 1 0 | 0 0 | 0 |  0 )
	//( 0 1 | 0 0 | 0 | -1 )
	//( 0 0 | 1 0 | 0 | -1 )
	//( 0 0 | 0 1 | 1 | -1 )
	//	=> IF (ii>=1 && jj==0 && b>=3) THEN (alpha,beta-1, ii-1, jj+b-1)
	//
	// 3) IF
	//(0/1 | alpha beta | ii jj | b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  0 | 1 | -1 )
	//( 1 | 0 0 |  0  1 | 0 | -1 )
	//( 1 | 0 0 | -1  0 | 0 |  0 )
	//( 1 | 0 0 |  0 -1 | 1 |  0 )
	// THEN
	//(alpha beta | ii jj | b | cnst)
	//( 1 0 | 0 0 | 0 | -1 )
	//( 0 1 | 0 0 | 0 |  0 )
	//( 0 0 | 1 0 | 1 | -1 )
	//( 0 0 | 0 1 | 0 | -1 )
	//	=> IF (ii==0 && jj>=1 && b>=3) THEN (alpha-1,beta, ii+b-1, jj-1)
	//
	// 4) IF
	//(0/1 | alpha beta | ii jj | b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  0 | 0 | -1 )
	//( 1 | 0 0 |  0  1 | 0 | -1 )
	//( 1 | 0 0 | -1  0 | 1 |  0 )
	//( 1 | 0 0 |  0 -1 | 1 |  0 )
	// THEN
	//(alpha beta | ii jj | b | cnst)
	//( 1 0 | 0 0 | 0 |  0 )
	//( 0 1 | 0 0 | 0 |  0 )
	//( 0 0 | 1 0 | 0 | -1 )
	//( 0 0 | 0 1 | 0 | -1 )
	//	=> IF (ii>=1 && jj>=1 && b>=3) THEN (alpha,beta, ii-1, jj-1)
}

// Example 2: (i,j -> i+j, i)
void test_MPP_Func_Ex2() {
	int nInd = 2;
	int nParam = 1;
	int dimOut = 2;
	
	int64** mat = (int64**) malloc(dimOut * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		mat[i] = (int64*) malloc((nInd+nParam+1)* sizeof(int64));
	mat[0][0] = 1; mat[0][1] = 1; mat[0][2] = 0; mat[0][3] = 0;
	mat[1][0] = 1; mat[1][1] = 0; mat[1][2] = 0; mat[1][3] = 0;
	affFuncMPP *affScalar = buildAffineFunction(mat, dimOut, nInd, nParam);
	
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 2;
	
	int64* scaleIm = (int64*) malloc(dimOut*sizeof(int64));
	scaleIm[0] = 1; scaleIm[1] = 2;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	opt->errorIfModulo = false;
	
	
	map<polyhedronMPP*, affFuncMPP*> resultFunc = getRectangularTiledFunction(affScalar, scale, scaleIm, opt);
	printoutFunction(resultFunc);
	
	// RESULT: 3 branches
	// 1) IF
	//(0/1 | alpha beta | ii jj | Nb Nl b | cnst)    // 0 = equality / 1 = inequality
	//( 1  | 0 0 |  1  1 | 0 0 0 |  0 )
	//( 1  | 0 0 |  1  0 | 0 0 0 |  0 )
	//( 1  | 0 0 | -1 -1 | 0 0 1 | -1 )
	//( 1  | 0 0 | -1  0 | 0 0 2 | -2 )
	//( 1  | 0 0 |  1  0 | 0 0 0 |  0 )
	//( 1  | 0 0 |  0  1 | 0 0 0 |  0 )
	//( 1  | 0 0 | -1  0 | 0 0 1 | -1 )
	//( 1  | 0 0 |  0 -1 | 0 0 2 | -1 )
	//( 1  | 0 0 |  0  0 | 0 0 1 | -3 )
	// THEN
	//(alpha beta | ii jj | Nb Nl b | cnst)
	//( 1 2 | 0 0 | 0 0 0 | 0 )
	//( 1 0 | 0 0 | 0 0 0 | 0 ) / 2
	//( 0 0 | 1 1 | 0 0 0 | 0 )
	//( 0 0 | 1 0 | 0 0 0 | 0 )
	// => IF (0<=ii+jj<=b-1 && 0<=ii<=2b-2 && 0<=ii<b && 0<=jj<2b && b>=3)
	//			THEN (alpha+2.beta, alpha/2, ii+jj, ii)
	//
	// 2) IF
	//(0/1 | alpha beta | ii jj | Nb Nl b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  1 | 0 0 -1 |  0 )
	//( 1 | 0 0 |  1  0 | 0 0  0 |  0 )
	//( 1 | 0 0 | -1 -1 | 0 0  2 | -1 )
	//( 1 | 0 0 | -1  0 | 0 0  2 | -2 )
	//( 1 | 0 0 |  1  0 | 0 0  0 |  0 )
	//( 1 | 0 0 |  0  1 | 0 0  0 |  0 )
	//( 1 | 0 0 | -1  0 | 0 0  1 | -1 )
	//( 1 | 0 0 |  0 -1 | 0 0  2 | -1 )
	//( 1 | 0 0 |  0  0 | 0 0  1 | -3 )
	// THEN
	//(alpha beta | ii jj | Nb Nl b | cnst)
	//( 1 2 | 0 0 | 0 0  0 | 1 )
	//( 1 0 | 0 0 | 0 0  0 | 0 ) / 2
	//( 0 0 | 1 1 | 0 0 -1 | 0 )
	//( 0 0 | 1 0 | 0 0  0 | 0 )
	// => IF (b<=ii+jj<=2b-1 && 0<=ii<=2b-2 && 0<=ii<b && 0<=jj<2b && b>=3)
	//			THEN (alpha+2.beta+1, alpha/2, ii+jj-b, ii)
	//
	// 3) IF
	//(0/1 | alpha beta | ii jj | Nb Nl b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | 0 0 |  1  1 | 0 0 -2 |  0 )
	//( 1 | 0 0 |  1  0 | 0 0  0 |  0 )
	//( 1 | 0 0 | -1 -1 | 0 0  3 | -1 )
	//( 1 | 0 0 | -1  0 | 0 0  2 | -2 )
	//( 1 | 0 0 |  1  0 | 0 0  0 |  0 )
	//( 1 | 0 0 |  0  1 | 0 0  0 |  0 )
	//( 1 | 0 0 | -1  0 | 0 0  1 | -1 )
	//( 1 | 0 0 |  0 -1 | 0 0  2 | -1 )
	//( 1 | 0 0 |  0  0 | 0 0  1 | -3 )
	// THEN
	//(alpha beta | ii jj | Nb Nl b | cnst)
	//( 1 2 | 0 0 | 0 0  0 | 2 )
	//( 1 0 | 0 0 | 0 0  0 | 0 ) / 2
	//( 0 0 | 1 1 | 0 0 -2 | 0 )
	//( 0 0 | 1 0 | 0 0  0 | 0 )
	// => IF (2b<=ii+jj<=3b-1 && 0<=ii<=2b-2 && 0<=ii<b && 0<=jj<2b && b>=3)
	//			THEN (alpha+2.beta+2, alpha/2, ii+jj-2b, ii)
}




/* ------------------------------------------ */

void test_LinAlg_1() {
	int64** A = (int64**) malloc(3*sizeof(int64*));
	for (int i=0; i<3; i++)
		A[i] = (int64*) malloc(3*sizeof(int64));
	A[0][0] = 1; A[0][1] = 0; A[0][2] = 1;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 1;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 1;
	
	int64** Ainv = inverseMatUnimod(A, 3, 3);
	
	printMatrix(Ainv, 3, 3);
}

void test_LinAlg_2() {
	int64** A = (int64**) malloc(2*sizeof(int64*));
	for (int i=0; i<2; i++)
		A[i] = (int64*) malloc(2*sizeof(int64));
	A[0][0] = 0; A[0][1] = 1;
	A[1][0] = 1; A[1][1] = 1;
	
	int64** Ainv = inverseMatUnimod(A, 2, 2);
	
	printMatrix(Ainv, 2, 2);
}

void test_LinAlg_3() {
	int nDim = 3;
	int64** A = (int64**) malloc(nDim*sizeof(int64*));
	for (int i=0; i<nDim; i++)
		A[i] = (int64*) malloc(nDim*sizeof(int64));
	A[0][0] = 1; A[0][1] = 1; A[0][2] = 1;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 2;
	
	rational64** invA = (rational64**) malloc(nDim*sizeof(rational64*));
	for (int i=0; i<nDim; i++)
		invA[i] = (rational64*) malloc(nDim*sizeof(rational64));
	for (int i=0; i<nDim; i++)
		for (int j=0; j<nDim; j++) {
			invA[i][j].num = 0;
			invA[i][j].den = 0;
		}
	
	int64 det = inverseDet(A, invA, nDim);
	
	cout << "Determinant = " << det << endl;
	printMatrix(invA, nDim, nDim);
}



/* ------------------------------------------ */

void test_lexminmax_1() {
	// Triangle {i,j | 0<=i && 0<=j && 0<=N-i-j}
	int nInd = 2;
	int nParam = 1;
	int nConstr = 3;
	
	int64** mat = (int64**) malloc(3*sizeof(int64*));
	for (int i=0; i<3; i++)
		mat[i] = (int64*) malloc(5*sizeof(int64));
	mat[0][0] = 1; mat[0][1] = 1; mat[0][2] = 0; mat[0][3] = 0; mat[0][4] = 0;
	mat[1][0] = 1; mat[1][1] = 0; mat[1][2] = 1; mat[1][3] = 0; mat[1][4] = 0;
	mat[2][0] = 1; mat[2][1] =-1; mat[2][2] =-1; mat[2][3] = 1; mat[2][4] = 0;
	
	polyhedronMPP* poly = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	rational64** matLexMax = lexmax(poly);
	printMatrix(matLexMax, poly->nInd, poly->nParam+1);
	// [ 1/1 0/1 ]
	// [ 0/1 0/1 ]
	// => (N,0) is the lexmax
	
	cout << endl;
	
	rational64** matLexMin = lexmin(poly);
	printMatrix(matLexMin, poly->nInd, poly->nParam+1);
	// [ 0/1 0/1 ]
	// [ 0/1 0/1 ]
	// => (0,0) is the lexmin
}

void test_lexminmax_2() {
	// Triangle {i,j | 0<=j && 0<=i-j && 0<=N-2i-j}
	int nInd = 2;
	int nParam = 1;
	int nConstr = 3;
	
	int64** mat = (int64**) malloc(3*sizeof(int64*));
	for (int i=0; i<3; i++)
		mat[i] = (int64*) malloc(5*sizeof(int64));
	mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 1; mat[0][3] = 0; mat[0][4] = 0;
	mat[1][0] = 1; mat[1][1] = 1; mat[1][2] = -1; mat[1][3] = 0; mat[1][4] = 0;
	mat[2][0] = 1; mat[2][1] =-2; mat[2][2] =-1; mat[2][3] = 1; mat[2][4] = 0;
	
	polyhedronMPP* poly = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	rational64** matLexMin = lexmin(poly);
	printMatrix(matLexMin, poly->nInd, poly->nParam+1);
	// [ 0/1 0/1 ]
	// [ 0/1 0/1 ]
	// => (0,0) is the lexmin
	
	rational64** matLexMax = lexmax(poly);
	printMatrix(matLexMax, poly->nInd, poly->nParam+1);
	// [ 1/2 0/1 ]
	// [ 0/1 0/1 ]
	// => (N/2, 0) is the rational lexmax
	
	cout << endl;
}

void test_minmax_1() {
	// Triangle {i,j | 0<=i && 0<=j && 0<=N-i-j}
	// Objective function: max (2i+j)
	int nInd = 2;
	int nParam = 1;
	int nConstr = 3;
	
	// poly = {i,j | 0<=i && 0<=j && 0<=N-i-j}
	int nConstr_mat = 3;
	int nCol_mat = nInd+nParam+2;
	int64** mat = (int64**) malloc(nConstr_mat*sizeof(int64*));
	for (int i=0; i<nConstr_mat; i++)
		mat[i] = (int64*) malloc(nCol_mat*sizeof(int64));
	mat[0][0] = 1; mat[0][1] = 1; mat[0][2] = 0; mat[0][3] = 0; mat[0][4] = 0;
	mat[1][0] = 1; mat[1][1] = 0; mat[1][2] = 1; mat[1][3] = 0; mat[1][4] = 0;
	mat[2][0] = 1; mat[2][1] =-1; mat[2][2] =-1; mat[2][3] = 1; mat[2][4] = 0;
	
	polyhedronMPP* poly = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	// obj = (i,j-> 2i+j)
	int nRow_matObj = 1;
	int nCol_matObj = nInd+nParam+1;
	int64** matObj = (int64**) malloc(nRow_matObj*sizeof(int64*));
	for (int i=0; i<nRow_matObj; i++)
		matObj[i] = (int64*) malloc(nCol_matObj*sizeof(int64));
	matObj[0][0] = 2; matObj[0][1] = 1; matObj[0][2] = 0; matObj[0][3] = 0;
	
	affFuncMPP* obj = buildAffineFunction(matObj, 1, nInd, nParam);
	
	rational64** matArgMax = argmax(obj, poly);
	printMatrix(matArgMax, poly->nInd, poly->nParam+1);
	// [ 1/1 0/1 ]
	// [ 0/1 0/1 ]
	// => maximum reached at (N,0)
	
	rational64** matMax = max(obj, poly);
	printMatrix(matMax, obj->dimOut, poly->nParam+1);
	// [ 2/1 0/1 ]
}

void test_shape_rect_1() {
	int nDim = 3;
	int64* sizes = (int64*) malloc(nDim * sizeof(int64));
	sizes[0] = 4; sizes[1] = 2; sizes[2] = 16;
	
	polyhedronMPP* shape = rectangularShape(sizes, nDim);
	printPolyhedronMPP(shape);
	//( 1  1  0  0  0  0 )
	//( 1  0  1  0  0  0 )
	//( 1  0  0  1  0  0 )
	//( 1 -1  0  0  4 -1 )
	//( 1  0 -1  0  2 -1 )
	//( 1  0  0 -1 16 -1 )
}

void test_shape_para_1() {
	int nDim = 3;
	int64** hyperplanes = (int64**) malloc(nDim * sizeof(int64*));
	for (int i=0; i<nDim; i++)
		hyperplanes[i] = (int64*) malloc(nDim * sizeof(int64));
	hyperplanes[0][0] =  1; hyperplanes[0][1] =  0; hyperplanes[0][2] =  0;
	hyperplanes[1][0] = -1; hyperplanes[1][1] =  1; hyperplanes[1][2] =  0;
	hyperplanes[2][0] =  0; hyperplanes[2][1] =  0; hyperplanes[2][2] =  1;
	
	int64* sizes = (int64*) malloc(nDim * sizeof(int64));
	sizes[0] = 4; sizes[1] = 2; sizes[2] = 16;
	
	
	polyhedronMPP* shape = parallelogramShape(hyperplanes, sizes, nDim);
	printPolyhedronMPP(shape);
	//( 1  1 -1  0  0  0 )
	//( 1  0  1  0  0  0 )
	//( 1  0  0  1  0  0 )
	//( 1 -1  1  0  4 -1 )
	//( 1  0 -1  0  2 -1 )
	//( 1  0  0 -1 16 -1 )
}


void test_MPP_Gen_Poly_Ex1() {
	// Example 1 (from "test_MPP_Poly_Ex1") : { i,j | N-i-j-1>=0 } with square tiles of size b*b
	int nInd = 2;
	int nParam = 1;
	int nConstr = 1;
	
	int64** mat = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		mat[i] = (int64*) malloc((2+nInd+nParam)* sizeof(int64));
	mat[0][0] = 1;
	mat[0][1] = -1;
	mat[0][2] = -1;
	mat[0][3] = 1;
	mat[0][4] = -1;
	polyhedronMPP *polyScalar = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 1;
	polyhedronMPP* shape = rectangularShape(scale, nInd);
	
	int64** lattice = (int64**) malloc( (nInd+1) * sizeof(int64*));
	for (int i=0; i<nInd+1; i++)
		lattice[i] = (int64*) malloc(nInd * sizeof(int64));
	lattice[0][0] = 1; lattice[0][1] = 0;
	lattice[1][0] = 0; lattice[1][1] = 1;
	lattice[2][0] = 1; lattice[2][1] = 1;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	
	list<list<polyhedronMPP*> > resultDom = getTiledDomain(polyScalar, shape, lattice, opt);
	printoutDomain(resultDom);
	// RESULT: Intersection of 1 union of 3 polyhedra:
	//(0/1 | alpha beta | ii jj | Nbl Nloc | b | cnst)    // 0 = equality / 1 = inequality
	//( 1 | -1 -1 |  0  0 | 1 0 | 0 | -2 )
	//( 1 |  0  0 |  1  0 | 0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 | 0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 | 0 0 | 1 | -3 )
	//		=> {alpha,beta, ii,jj | alpha+beta<=Nbl-2 && 0<=(ii,jj)<b && b>=3}
	// and
	//( 0 | -1 -1 |  0  0 |  1 0 | 0 | -1 )
	//( 1 |  0  0 | -1 -1 |  0 1 | 1 | -1 )
	//( 1 |  0  0 |  1  0 |  0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 |  0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 |  0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 |  0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 |  0 0 | 1 | -3 )
	//		 => {alpha,beta, ii,jj | alpha+beta=Nbl-1 && Nloc+b-1>=ii+jj && 0<=(ii,jj)<b && b>=3}
	// and
	//( 0 | -1 -1 |  0  0 | 1 0 | 0 |  0 )
	//( 1 |  0  0 | -1 -1 | 0 1 | 0 | -1 )
	//( 1 |  0  0 |  1  0 | 0 0 | 0 |  0 )
	//( 1 |  0  0 |  0  1 | 0 0 | 0 |  0 )
	//( 1 |  0  0 | -1  0 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0 -1 | 0 0 | 1 | -1 )
	//( 1 |  0  0 |  0  0 | 0 0 | 1 | -3 )
	//		=> {alpha,beta, ii,jj | alpha+beta=Nbl && Nloc-1>=ii+jj && 0<=(ii,jj)<b && b>=3}
}


void test_MPP_Gen_Func_Ex1() {
	// Example 1: (i,j -> i, j) with:
	//   - square 3b*2b tiles on the input space
	//   - hexagonal 45° "4b*2b" tiles on the output space
	
	int nInd = 2;
	int nParam = 0;
	int dimOut = 2;
	
	int64** mat = (int64**) malloc(dimOut * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		mat[i] = (int64*) malloc((nInd+nParam+1)* sizeof(int64));
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	affFuncMPP *affScalar = buildAffineFunction(mat, dimOut, nInd, nParam);
	
	// Input space partitioning
	int64* scale = (int64*) malloc(nInd*sizeof(int64));
	scale[0] = 1; scale[1] = 1;
	polyhedronMPP *shapeIn = rectangularShape(scale, 2);
	int64** latticeIn = rectangularOriginLattice(scale, 2);
	
	// Output space partitioning
	// Tile shape: {ii,jj | 0<=ii-jj<4b && -b<=jj<b && 0<=ii+jj<4b }
	int nConstrHex = 6;
	int nCol_matHex = 2 + dimOut + 1;
	int64** matHex = (int64**) malloc(nConstrHex * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		matHex[i] = (int64*) malloc(nCol_matHex * sizeof(int64));
	matHex[0][0] = 1; matHex[0][1] = -1; matHex[0][2] =  1; matHex[0][3] =  4; matHex[0][4] = -1;	// (0<=-ii+jj+4b-1)
	matHex[1][0] = 1; matHex[1][1] =  1; matHex[1][2] = -1; matHex[1][3] =  0; matHex[1][4] =  0;	// (0<=ii-jj)
	matHex[2][0] = 1; matHex[2][1] =  0; matHex[2][2] =  1; matHex[2][3] = -1; matHex[2][4] =  0;	// (0<=jj-b)
	matHex[3][0] = 1; matHex[3][1] =  0; matHex[3][2] = -1; matHex[3][3] =  1; matHex[3][4] = -1;	// (0<=b-jj-1)
	matHex[4][0] = 1; matHex[4][1] =  1; matHex[4][2] =  1; matHex[4][3] =  0; matHex[4][4] =  0;	// (0<=ii+jj)
	matHex[5][0] = 1; matHex[5][1] = -1; matHex[5][2] = -1; matHex[5][3] =  4; matHex[5][4] = -1;	// (0<=4b-ii-jj-1)
	polyhedronMPP* shapeOut = buildPolyhedron(matHex, nConstrHex, dimOut, 1);
	
	// Lattice of tile origin - basis is: (3b, b), (3b, -b)
	int64** latticeOut = (int64**) malloc((dimOut+1) * sizeof(int64*));
	for (int i=0; i<dimOut+1; i++)
		latticeOut[i] = (int64*) malloc(dimOut * sizeof(int64));
	latticeOut[0][0] = 3; latticeOut[0][1] =  3;
	latticeOut[1][0] = 1; latticeOut[1][1] = -1;
	latticeOut[2][0] = 1; latticeOut[2][1] =  1;
	
	optionMPP* opt = (optionMPP*) malloc(sizeof(optionMPP));
	opt->kMinMaxOption = 0;
	opt->areParamDiv = false;
	opt->minBlSizeParam = 3;
	
	map<polyhedronMPP*, affFuncMPP*> resultFunc = getTiledFunction(affScalar,
			shapeIn, latticeIn, shapeOut, latticeOut, opt);
	printoutFunction(resultFunc);
	
	// TODO: test that
	
	
}


/* ------------------------------------------ */

int main() {
	//test_MPP_Poly_Ex1();
	//test_MPP_Poly_Ex2();
	//test_MPP_Func_Ex1();
	//test_MPP_Func_Ex2();
	
	//test_LinAlg_1();
	//test_LinAlg_2();
	//test_LinAlg_3();
	
	
	//test_lexminmax_1();
	//test_lexminmax_2();
	//test_minmax_1();
	
	//test_shape_rect_1();
	//test_shape_para_1();
	
	test_MPP_Gen_Poly_Ex1();
	//test_MPP_Gen_Func_Ex1();
	
	// TODO: other test for the polyhedron case
	// TODO: test for the affine function case
	
	return 0;
}

