// Test of the CART functions
// Copyright Guillaume Iooss, 2014, All right reserved.

#include "MPP_rect.h"
#include "lexmin.h"


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
	mat[1][0] = 1; mat[1][1] =-1; mat[1][2] =-1; mat[1][3] = 1; mat[1][4] = 0;
	
	polyhedronMPP* poly = buildPolyhedron(mat, nConstr, nInd, nParam);
	
	rational64** matLexMax = lexmax(poly);
	printMatrix(matLexMax, poly->nInd, poly->nParam+1);
	// [ 0/1 0/1 ]
	// [ 1/1 0/1 ]
	// => (0,N) is the lexmax
	
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
	mat[1][0] = 1; mat[1][1] =-2; mat[1][2] =-1; mat[1][3] = 1; mat[1][4] = 0;
	
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


/* ------------------------------------------ */

int main() {
	//test_MPP_Poly_Ex1();
	//test_MPP_Poly_Ex2();
	//test_MPP_Func_Ex1();
	//test_LinAlg_1();
	//test_LinAlg_2();
	
	//test_lexminmax_1();
	test_lexminmax_2();
	
	return 0;
}
