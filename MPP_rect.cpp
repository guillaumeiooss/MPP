// Standalone C++ version of Monoparametric tiling
// Author: Guillaume Iooss

#include "MPP_rect.h"

// Rectangular case
list<list<polyhedronMPP*> > getRectangularTiledDomain(polyhedronMPP *polyScalar, int64 *scale, optionMPP *option) {
	
	// Assertions
	assert(option->minBlSizeParam>=1);
	
	int nParam = polyScalar->nParam;
	int nInd = polyScalar->nInd;
	int nConstr = polyScalar->nConstr;
	
	// DEBUG
	//cout << "nConstr = " << nConstr << " | nInd = " << nInd << " | nParam = " << nParam << endl;
	
	// Returned structure
	list<list<polyhedronMPP*> > resDom;
	
	// Trivial case: the polyhedron do not have any index
	//	=> We return a similar polyhedron (with the right amount of parameters)
	if (nInd==0) {
		int64** matConst = (int64**) malloc(1 * sizeof(int64*));
		matConst[0] = (int64*) malloc( (3+2*polyScalar->nParam) * sizeof(int64));
		polyhedronMPP* nullPoly = buildPolyhedron(matConst, 1, 0, 1+2*nParam);
		
		list<polyhedronMPP*> templist;
		templist.push_back(nullPoly);
		resDom.push_back(templist);
		
		return resDom;
	}
	
	// Trivial case: the polyhedron do not have any constraint
	//	=> We return the universe
	if (nConstr==0) {
		int64** matConst = (int64**) malloc(0 * sizeof(int64*));
		polyhedronMPP* nullPoly = buildPolyhedron(matConst, 0, nInd, 1+2*nParam);
		
		list<polyhedronMPP*> templist;
		templist.push_back(nullPoly);
		resDom.push_back(templist);
		
		return resDom;
	}
	
	// Extraction of the information
	polyhedronMPP* polyScalar_noEq = eliminateEqualities(polyScalar);
	nConstr = polyScalar_noEq->nConstr;
	
	int64* ineqEqPart = (int64*) malloc(nConstr * sizeof(int64));
	int64** linPart = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		linPart[i] = (int64*) malloc(nInd * sizeof(int64));
	int64** paramPart = (int64**) malloc(nConstr * sizeof(int64*));
	for (int i=0; i<nConstr; i++)
		paramPart[i] = (int64*) malloc(nParam * sizeof(int64));
	int64* constPart = (int64*) malloc(nConstr * sizeof(int64));
	extractPoly(polyScalar_noEq, ineqEqPart, linPart, paramPart, constPart);
	
	// Computation of kmin/kmax
	vector<long> kmax(nConstr);
	vector<long> kmin(nConstr);
	
	for (int i=0; i<nConstr; i++) {
		long tempkmin = 0;		// ( \sum_i QD-_i + \sum_i Q^p-_i )
		long tempkmax = 0;		// ( \sum_i QD+_i + \sum_i Q^p+_i )
		
		for (int j=0; j<nInd; j++) {
			if (linPart[i][j]>0)
				tempkmax += scale[j] * linPart[i][j];
			else if (linPart[i][j]<0)
				tempkmin += scale[j] * linPart[i][j];
		}
		
		for (int j=0; j<nParam; j++) {
			if (paramPart[i][j]>0)
				tempkmax += paramPart[i][j];
			else if (paramPart[i][j]<0)
				tempkmin += paramPart[i][j];
		}
		
		// kmax
		if (constPart[i]-tempkmax>=0) {
			if (option->kMinMaxOption==0)
				kmax[i] = tempkmax + (long) floor( ((double) (constPart[i] - tempkmax)) / ((double) option->minBlSizeParam) );
			else
				kmax[i] = tempkmax;
		} else
			kmax[i] = tempkmax-1;
		
		// kmin
		if (constPart[i]-tempkmin>=0)
			kmin[i] = tempkmin;
		else {
			if (option->kMinMaxOption==0)
				kmin[i] = tempkmin + (long) floor( ((double) (constPart[i] - tempkmin)) / ((double) option->minBlSizeParam) );
			else
				kmin[i] = tempkmin-1;
		}
	}
	
	
	/* DEBUG - kmin / kmax
	cout << "DEBUG: kmin - kmax" << endl;
	for (int i=0; i<nConstr; i++)
		cout << "kmin[" << i << "] = " << kmin[i] << "  |  kmax[" << i << "] = " << kmax[i] << endl;
	cout << endl;
	//*/
	
	// We consider the ith constraint (outer intersection)
	for (int i=0; i<nConstr; i++) {
		list<polyhedronMPP*> domIthList;
		
		if (option->areParamDiv) {
			// We consider all the different cases arising with the constraint (inner union)
			for (long k=kmin[i]+1; k<=kmax[i]; k++) {
				// * Case k > kmin
				
				// In Polylib format, the matrix of constraint is:
				//
				// ( eq |  \alpha | ii | \rho   |  pp  |  b  | const)    <= Row to know which column corresponds to what...
				//
				// [ 0    Q_i,*.D   0    Qp_i,*    0      0     k   ]
				// [ 1      0       Q      0     Qp_i,*  -k     q   ]
				// [ 1      0      Id      0       0      0     0   ]
				// [ 1      0     -Id      0       0     D.1   -1   ]
				// [ 0      0       0      0      Id      0     0   ]
				// [ 1      0       0      0       0      1   -minb ]
				//
				// First column: 0 = equality / 1 = inequality
				// \alpha = blocked index parameter / ii = local index parameter
				// \rho = blocked parameters / pp = local parameters / b = block size parameter
				
				int nRow_blConstr = 3+2*nInd+nParam;
				int nColumn_blConstr = 3+2*nParam+2*nInd;
				int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
				for (int l=0; l<nRow_blConstr; l++)
				 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
				for (int i=0; i<nRow_blConstr; i++)
					for (int j=0; j<nColumn_blConstr; j++)
						blockedConstr[i][j] = 0;
				
				// (First line)
				for (int j=0; j<nInd; j++)
					blockedConstr[0][1+j] = linPart[i][j] * ((long) scale[j]);
				for (int j=0; j<nParam; j++) {
					blockedConstr[0][1+2*nInd+j] = paramPart[i][j];
				}
				blockedConstr[0][nColumn_blConstr-1] = k;
				
				// (Second line)
				blockedConstr[1][0] = 1;
				for (int j=0; j<nInd; j++)
					blockedConstr[1][1+nInd+j] = linPart[i][j];
				for (int j=0; j<nParam; j++)
					blockedConstr[1][1+2*nInd+nParam+j] = paramPart[i][j];
				blockedConstr[1][1+2*nInd+2*nParam] = -k;
				blockedConstr[1][nColumn_blConstr-1] = constPart[i];
				
				// (Third lines)
				for (int j=0; j<nInd; j++) {
					blockedConstr[2+j][0] = 1;
					blockedConstr[2+j][1+nInd+j] = 1;
				}
				
				// (Fourth lines)
				for (int j=0; j<nInd; j++) {
					blockedConstr[2+nInd+j][0] = 1;
					blockedConstr[2+nInd+j][1+nInd+j] = -1;
					blockedConstr[2+nInd+j][1+2*nInd+2*nParam] = ((long) scale[j]);
					blockedConstr[2+nInd+j][nColumn_blConstr-1] = -1;
				}
				
				// (Fifth lines)
				for (int j=0; j<nParam; j++) {
					blockedConstr[2+2*nInd+j][1+2*nInd+nParam+j] = 1;
				}
				
				// (Sixth line)
				blockedConstr[2+2*nInd+nParam][0] = 1;
				blockedConstr[2+2*nInd+nParam][1+2*nInd+2*nParam] = 1;
				blockedConstr[2+2*nInd+nParam][nColumn_blConstr-1] = -option->minBlSizeParam;
				
				polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
				domIthList.push_back(polyRet);
			}
			
			
			
			// * Case k = kmin
			
			// In Polylib format, the matrix of constraint is:
			//
			// ( eq | \alpha  | ii | \rho   |  pp  |  b  | const)    <= Row to know which column corresponds to what...
			//
			// [ 1    Q_i,*.D   0    Qp_i,*     0     0    kmin ]
			// [ 1      0       Id     0        0     0     0   ]
			// [ 1      0      -Id     0        0    D.1   -1   ]
			// [ 1      0       0      0        0     1   -minb ]
			// [ 0      0      0       0       Id     0     0   ]
			//
			// First column: 0=equality / 1 = inequality
			// \alpha = blocked index parameter / ii = local index parameter
			// \rho = blocked parameters / pp = local parameters / b = block size parameter
			
			int nRow_blConstr = 2+2*nInd+nParam;
			int nColumn_blConstr = 3+2*nParam+2*nInd;
			int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
			for (int l=0; l<nRow_blConstr; l++)
			 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
			for (int i=0; i<nRow_blConstr; i++)
				for (int j=0; j<nColumn_blConstr; j++)
					blockedConstr[i][j] = 0;
			
			// (First line)
			blockedConstr[0][0] = 1;
			for (int j=0; j<nInd; j++)
				blockedConstr[0][1+j]= linPart[i][j] * ((long) scale[j]);
			for (int j=0; j<nParam; j++)
				blockedConstr[0][1+2*nInd+j] = paramPart[i][j];
			blockedConstr[0][nColumn_blConstr-1] = kmin[i];
			
			// (Second lines)
			for (int j=0; j<nInd; j++) {
				blockedConstr[1+j][0] = 1;
				blockedConstr[1+j][1+nInd+j] = 1;
			}
			
			// (Third lines)
			for (int j=0; j<nInd; j++) {
				blockedConstr[1+nInd+j][0] = 1;
				blockedConstr[1+nInd+j][1+nInd+j] = -1;
				blockedConstr[1+nInd+j][1+2*nInd+2*nParam] = ((long) scale[j]);
				blockedConstr[1+nInd+j][nColumn_blConstr-1] = -1;
			}
			
			// (Fourth line)
			blockedConstr[1+2*nInd][0] = 1;
			blockedConstr[1+2*nInd][1+2*nInd+2*nParam] = 1;
			blockedConstr[1+2*nInd][nColumn_blConstr-1] = -option->minBlSizeParam;
			
			// (Fifth lines)
			for (int j=0; j<nParam; j++) {
				blockedConstr[2+2*nInd+j][1+2*nInd+nParam+j] = 1;
			}
			
			polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
			domIthList.push_back(polyRet);
		} else {
			
			// We consider all the different cases arising with the constraint (inner union)
			for (long k=kmin[i]+1; k<=kmax[i]; k++) {
				// * Case k > kmin
				
				// In Polylib format, the matrix of constraint is:
				//
				// ( eq | \alpha | ii |  \rho  |  pp  |  b   | const)    <= Row to know which column corresponds to what...
				//
				// [ 0    Q_i,*.D   0    Qp_i,*    0      0     k   ]
				// [ 1      0       Q      0     Qp_i,*  -k     q   ]
				// [ 1      0      Id      0       0      0     0   ]
				// [ 1      0     -Id      0       0     D.1   -1   ]
				// [ 1      0       0      0       0      1   -minb ]
				// First column: 0=equality / 1 = inequality
				// \alpha = blocked index parameter / ii = local index parameter
				// \rho = blocked parameters / pp = local parameters / b = block size parameter
				
				int nRow_blConstr = 3+2*nInd;
				int nColumn_blConstr = 3+2*nParam+2*nInd;
				
				int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
				for (int l=0; l<nRow_blConstr; l++)
				 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
				for (int i=0; i<nRow_blConstr; i++)
					for (int j=0; j<nColumn_blConstr; j++)
						blockedConstr[i][j] = 0;
				
				// (First line)
				for (int j=0; j<nInd; j++)
					blockedConstr[0][1+j] = linPart[i][j] * ((long) scale[j]);
				for (int j=0; j<nParam; j++) {
					blockedConstr[0][1+2*nInd+j] = paramPart[i][j];
				}
				blockedConstr[0][nColumn_blConstr-1] = k;
				
				// (Second line)
				blockedConstr[1][0] = 1;
				for (int j=0; j<nInd; j++)
					blockedConstr[1][1+nInd+j] = linPart[i][j];
				for (int j=0; j<nParam; j++)
					blockedConstr[1][1+2*nInd+nParam+j] = paramPart[i][j];
				blockedConstr[1][1+2*nInd+2*nParam] = -k;
				blockedConstr[1][nColumn_blConstr-1] = constPart[i];
				
				// (Third lines)
				for (int j=0; j<nInd; j++) {
					blockedConstr[2+j][0] = 1;
					blockedConstr[2+j][1+nInd+j] = 1;
				}
				
				// (Fourth lines)
				for (int j=0; j<nInd; j++) {
					blockedConstr[2+nInd+j][0] = 1;
					blockedConstr[2+nInd+j][1+nInd+j] = -1;
					blockedConstr[2+nInd+j][1+2*nInd+2*nParam] = ((long) scale[j]);
					blockedConstr[2+nInd+j][nColumn_blConstr-1] = -1;
				}
				
				// (Fifth line)
				blockedConstr[2+2*nInd][0] = 1;
				blockedConstr[2+2*nInd][1+2*nInd+2*nParam] = 1;
				blockedConstr[2+2*nInd][nColumn_blConstr-1] = -option->minBlSizeParam;
				
				polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
				domIthList.push_back(polyRet);
			}
			
			// * Case k = kmin
			
			// In Polylib format, the matrix of constraint is:
			//
			// ( eq | \alpha | ii | \rho   |  pp  |  b  |  const)    <= Row to know which column corresponds to what...
			//
			// [ 1   Q_i,*.D   0    Qp_i,*    0      0     kmin ]
			// [ 1      0      Id     0       0      0      0   ]
			// [ 1      0     -Id     0       0     D.1    -1   ]
			// [ 1      0      0      0       0      1    -minb ]
			// First column: 0 = equality / 1 = inequality
			// \alpha = blocked index parameter / ii = local index parameter
			// \rho = blocked parameters / pp = local parameters / b = block size parameter
			
			int nRow_blConstr = 2+2*nInd;
			int nColumn_blConstr = 3+2*nParam+2*nInd;
			int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
			for (int l=0; l<nRow_blConstr; l++)
			 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
			for (int i=0; i<nRow_blConstr; i++)
				for (int j=0; j<nColumn_blConstr; j++)
					blockedConstr[i][j] = 0;
			
			// (First line)
			blockedConstr[0][0] = 1;
			for (int j=0; j<nInd; j++)
				blockedConstr[0][1+j]= linPart[i][j] * ((long) scale[j]);
			for (int j=0; j<nParam; j++)
				blockedConstr[0][1+2*nInd+j] = paramPart[i][j];
			blockedConstr[0][nColumn_blConstr-1] = kmin[i];
			
			// (Second lines)
			for (int j=0; j<nInd; j++) {
				blockedConstr[1+j][0] = 1;
				blockedConstr[1+j][1+nInd+j] = 1;
			}
			
			// (Third lines)
			for (int j=0; j<nInd; j++) {
				blockedConstr[1+nInd+j][0] = 1;
				blockedConstr[1+nInd+j][1+nInd+j] = -1;
				blockedConstr[1+nInd+j][1+2*nInd+2*nParam] = ((long) scale[j]);
				blockedConstr[1+nInd+j][nColumn_blConstr-1] = -1;
			}
			
			// (Fourth line)
			blockedConstr[1+2*nInd][0] = 1;
			blockedConstr[1+2*nInd][1+2*nParam+2*nInd] = 1;
			blockedConstr[1+2*nInd][nColumn_blConstr-1] = -option->minBlSizeParam;
			
			polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
			domIthList.push_back(polyRet);
		}
		
		resDom.push_back(domIthList);
	}
	
	// Free temporary structures
	freePolyhedron(polyScalar_noEq);
	free(ineqEqPart);
	freeMatrix(linPart, nConstr);
	freeMatrix(paramPart, nConstr);
	free(constPart);
	
	return resDom;
}



map<polyhedronMPP*, affFuncMPP*> getRectangularTiledFunction(affFuncMPP *affScalar, int64 *scale, int64 *scaleIm, optionMPP *option) {
	assert(option->minBlSizeParam>=1);
	
	int nParam = affScalar->nParam;
	int nInd = affScalar->nInd;
	int dimOut = affScalar->dimOut;
	
	// DEBUG
	//cout << "dimOut = " << dimOut << " | nInd = " << nInd << " | nParam = " << nParam << endl;
	
	// Returned data structure
	map<polyhedronMPP*, affFuncMPP*> result;
	
	
	// Special case: image dimension is of dimension 0
	if (affScalar->dimOut==0) {
		int64** matConst = (int64**) malloc(1 * sizeof(int64*));
		matConst[0] = (int64*) malloc( (3 + 2*affScalar->nInd + 2*affScalar->nParam) * sizeof(int64));
		polyhedronMPP* nullPoly = buildPolyhedron(matConst, 1, 2*affScalar->nInd, 2*affScalar->nParam+1);
		
		int64** matAffFunc = (int64**) malloc(0 * sizeof(int64*));
		affFuncMPP* affFunc = buildAffineFunction(matAffFunc, 0, 2*affScalar->nInd, 2*affScalar->nParam+1);
		
		result.insert ( pair<polyhedronMPP*, affFuncMPP*>(nullPoly, affFunc) );
		return result;
	}
	
	// Extraction of the informations
	int64** linPart = (int64**) malloc(dimOut * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		linPart[i] = (int64*) malloc(nInd * sizeof(int64));
	int64** paramPart = (int64**) malloc(dimOut * sizeof(int64*));
	for (int i=0; i<dimOut; i++)
		paramPart[i] = (int64*) malloc(nParam * sizeof(int64));
	int64* constPart = (int64*) malloc(dimOut * sizeof(int64));
	extractAffFunc(affScalar, linPart, paramPart, constPart);
	
	
	// We compute kmin and kmax
	vector<long> kmax(affScalar->dimOut);
	vector<long> kmin(affScalar->dimOut);
	
	for (int i=0; i<dimOut; i++) {
		long tempkmin = 0;		// ( \sum_i QD-_i + \sum_i Q^p-_i )
		long tempkmax = 0;		// ( \sum_i QD+_i + \sum_i Q^p+_i )
		
		for (int j=0; j<nInd; j++) {
			if (linPart[i][j]>0)
				tempkmax += scale[j] * linPart[i][j];
			else if (linPart[i][j]<0)
				tempkmin += scale[j] * linPart[i][j];
		}
		
		for (int j=0; j<nParam; j++) {
			if (paramPart[i][j]>0)
				tempkmax += paramPart[i][j];
			else if (paramPart[i][j]<0)
				tempkmin += paramPart[i][j];
		}
		
		// kmax
		if (constPart[i]-tempkmax>=0) {
			if (option->kMinMaxOption==0)
				kmax[i] = tempkmax + (long) floor( ((double) (constPart[i] - tempkmax)) / ((double) option->minBlSizeParam) );
			else
				kmax[i] = tempkmax;
		} else
			kmax[i] = tempkmax-1;
		
		// kmin
		if (constPart[i]-tempkmin>=0)
			kmin[i] = tempkmin;
		else {
			if (option->kMinMaxOption==0)
				kmin[i] = tempkmin + (long) floor( ((double) (constPart[i] - tempkmin)) / ((double) option->minBlSizeParam) );
			else
				kmin[i] = tempkmin-1;
		}
	}
	/* DEBUG - kmin / kmax
	for (int i=0; i<dimOut; i++)
		cout << "kmin[" << i << "] = " << kmin[i] << "  |  kmax[" << i << "] = " << kmax[i] << endl;
	cout << endl;
	//*/
	
	// Now, we build our new blocked affine functions: we have a part of piece-wise affine function per \vec{k},
	//		which is: \phi(\alpha,ii) = (Q.D.\alpha + Q_p.\rho + k  ,  Q.ii + q - b.k) when  b.k <= Q.ii+q < b.(k+1)
	// Thus, we iterate on the multi-dimensional vector k
	vector<long> kCurr(dimOut);		// Multi-dimensional iterator
	for (int i=0; i<dimOut; i++)
		kCurr[i] = kmin[i];
	
	while (kCurr[dimOut-1]<=kmax[dimOut-1]) {
		// We build the piece-wise part corresponding to the vector "kCurr"
		// In Polylib format, the matrix of input constraints is:
		// 
		// ( eq | \alpha | ii | \rho  | pp |    b   |  const )    <= Row to know which column corresponds to what...
		// 
		// [ 1      0       Q     0     Qp    -D'.k       q     ]
		// [ 1      0      -Q     0    -Qp   D'.(k+1)  -D'.1-q  ]
		// [ 1      0      Id     0      0      0         0     ]
		// [ 1      0     -Id     0      0     D.1       -1     ]
		// [ 1      0       0     0      0      1       -minb   ]
		// If divisible: 
		// [ 0      0       0     0     Id      0         0     ]
		// 
		// First column: 0=equality / 1 = inequality
		// \alpha = blocked index parameter / ii = local index parameter
		// \rho = blocked parameters / pp = local parameters / b = block size parameter
		//
		int nRow_blConstr = (option->areParamDiv)? 1+2*dimOut+2*nInd+nParam : 1+2*dimOut+2*nInd;
		int nColumn_blConstr = 3+2*nParam+2*nInd;
		
		int64** inputConstrLongMat = (int64**) malloc(nRow_blConstr * sizeof(int64*));
		for (int l=0; l<nRow_blConstr; l++)
		 	inputConstrLongMat[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
		for (int i=0; i<nRow_blConstr; i++)
				for (int j=0; j<nColumn_blConstr; j++)
					inputConstrLongMat[i][j] = 0;
		
		// (First rows)
		for (int i=0; i<dimOut; i++) {
			inputConstrLongMat[i][0] = 1;
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[i][1+nInd+j] = linPart[i][j];
			for (int j=0; j<nParam;j++)
				inputConstrLongMat[i][1+2*nInd+nParam+j] = paramPart[i][j];
			inputConstrLongMat[i][1+2*nInd+2*nParam] = -kCurr[i] * scaleIm[i];
			inputConstrLongMat[i][nColumn_blConstr-1] = constPart[i];
		}
		
		// (Second rows)
		for (int i=dimOut; i<2*dimOut; i++) {
			inputConstrLongMat[i][0] = 1;
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[i][1+nInd+j] = -linPart[i-dimOut][j];
			for (int j=0; j<nParam; j++)
				inputConstrLongMat[i][1+2*nInd+nParam+j] = -paramPart[i-dimOut][j];
			inputConstrLongMat[i][1+2*nInd+2*nParam] = scaleIm[i-dimOut]* ( kCurr[i-dimOut]+1 );
			inputConstrLongMat[i][nColumn_blConstr-1] = -scaleIm[i-dimOut]-constPart[i-dimOut];
		}
		
		// (Third rows)
		for (int i=0; i<nInd; i++) {
			inputConstrLongMat[i+2*dimOut][0] = 1;
			inputConstrLongMat[i+2*dimOut][1+nInd+i] = 1;
		}
		
		// (Fourth rows)
		for (int i=0; i<nInd; i++) {
			inputConstrLongMat[i+2*dimOut+nInd][0] = 1;
			inputConstrLongMat[i+2*dimOut+nInd][1+nInd+i] = -1;
			inputConstrLongMat[i+2*dimOut+nInd][1+2*nInd+2*nParam] = (long) scale[i];
			inputConstrLongMat[i+2*dimOut+nInd][nColumn_blConstr-1] = -1;
		}
		
		// (Fifth rows)
		inputConstrLongMat[2*dimOut+2*nInd][0] = 1;
		inputConstrLongMat[2*dimOut+2*nInd][1+2*nInd+2*nParam] = 1;
		inputConstrLongMat[2*dimOut+2*nInd][nColumn_blConstr-1] = -option->minBlSizeParam;
		
		// (Sixth rows)
		if (option->areParamDiv) {
			for (int i=0; i<nParam; i++) {
				inputConstrLongMat[2*dimOut+2*nInd+i][1+2*nInd+nParam+i] = 1;
			}
		}
		
		/* DEBUG
		cout << " * inputConstrLongMat:" << endl;
		for (int i=0; i<nRow_blConstr; i++) {
			for (int j=0; j<nColumn_blConstr; j++)
				cout << inputConstrLongMat[i][j] << " ";
			cout << endl;
		}
		cout << endl;
		//*/
		polyhedronMPP* polyRet = buildPolyhedron(inputConstrLongMat, nRow_blConstr, 2*nInd, 2*nParam+1);
		
		
		// The matrix of the affine function is:
		// 
		// (    \alpha    | ii |     \rho    |   pp   |     b      | const)    <= Row to know which column corresponds to what...
		//
		// [  D'^{-1}Q.D     0   D'^{-1}.Qp      0          0         k  ]
		// [      0          Q       0           Qp    -D'^{-1}k      q  ]
		//
		int nRow_relConstr = 2*dimOut;
		int nColumn_relConstr = 2+2*nParam+2*nInd;
		
		int64** relationConstrLongMat = (int64**) malloc(nRow_relConstr * sizeof(int64*));
		for (int l=0; l<nRow_relConstr; l++)
		 	relationConstrLongMat[l] = (int64*) malloc(nColumn_relConstr * sizeof(int64));
		for (int i=0; i<nRow_relConstr; i++)
				for (int j=0; j<nColumn_relConstr; j++)
					relationConstrLongMat[i][j] = 0;
		
		int64* divConstrLongMat = (int64*) malloc(nRow_relConstr * sizeof(int64));
		for (int i=0; i<nRow_relConstr; i++)
			divConstrLongMat[i] = 1;
		
		// (First rows)
		for (int i=0; i<dimOut; i++) {
			if (option->errorIfModulo) {
				for (int j=0; j<nInd; j++) {
					long temp = (long) (linPart[i][j] * scale[j] / scaleIm[i]);
					if (temp*scaleIm[i] != linPart[i][j] * scale[j]) {
						cerr << "The resulting affine function has Z-polyhedral constraints (linear | i = " << i << " | j = " << j << " )" << endl;
						exit(-1);
					}
					relationConstrLongMat[i][j] = temp;
				}
				for (int j=0; j<nParam; j++) {
					long temp = (long) (paramPart[i][j]/scaleIm[i]);
					if (temp*scaleIm[i] != paramPart[i][j] && option-> errorIfModulo) {
						cerr << "The resulting affine function has Z-polyhedral constraints (param | i = " << i << " | j = " << j << " )" << endl;
						exit(-1);
					}
					relationConstrLongMat[i][2*nInd+j] = temp;
				}
				relationConstrLongMat[i][nColumn_relConstr-1] = kCurr[i];
			} else {
				// Value in this case:
				// [  Q.D     0   Qp      0          0         D'.k  ] / div: D'
				for (int j=0; j<nInd; j++)
					relationConstrLongMat[i][j] = linPart[i][j] * scale[j];
				for (int j=0; j<nParam; j++)
					relationConstrLongMat[i][2*nInd+j] = paramPart[i][j];
				relationConstrLongMat[i][nColumn_relConstr-1] = scaleIm[i] * kCurr[i];
				
				divConstrLongMat[i] = scaleIm[i];
			}
		}
		
		// (Second rows)
		for (int i=dimOut; i<2*dimOut; i++) {
			for (int j=0; j<nInd; j++)
				relationConstrLongMat[i][nInd+j] = linPart[i-dimOut][j];
			for (int j=0; j<nParam;j++)
				relationConstrLongMat[i][2*nInd+nParam+j] = paramPart[i-dimOut][j];
			relationConstrLongMat[i][2*nInd+2*nParam] = -scaleIm[i-dimOut]*kCurr[i-dimOut];
			relationConstrLongMat[i][nColumn_relConstr-1] = constPart[i-dimOut];
		}
		
		/* DEBUG
		cout << " * relationConstrLongMat:" << endl;
		for (int i=0; i<nRow_relConstr; i++) {
			for (int j=0; j<nColumn_relConstr; j++)
				cout << relationConstrLongMat[i][j] << " ";
			cout << endl;
		}
		cout << endl;
		cout << " * divConstrLongMat:" << endl;
		for (int i=0; i<nRow_relConstr; i++)
			cout << divConstrLongMat[i] << " ";
		cout << endl;
		//*/
		
		affFuncMPP* affFuncRet = buildRationalAffineFunction(relationConstrLongMat, divConstrLongMat, 2*dimOut, 2*nInd, 2*nParam+1);
		simplifyAffFuncMPP(affFuncRet);
		result.insert ( pair<polyhedronMPP*, affFuncMPP*>(polyRet, affFuncRet) );
		
		// We increase the multi-dimensional iterator, starting from the first dimension and propagating the overflows
		kCurr[0]++;
		for (int i=0; i<dimOut-1; i++)
			if (kCurr[i]>kmax[i]) {
				kCurr[i] = kmin[i];
				kCurr[i+1]++;
			}
	} // End of multi-dimensional loop
	
	// Free temporary structures
	freeMatrix(linPart, dimOut);
	freeMatrix(paramPart, dimOut);
	free(constPart);
	
	return result;
}


