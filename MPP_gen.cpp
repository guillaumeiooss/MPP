// Standalone C++ version of Monoparametric tiling
// Author: Guillaume Iooss

#include "MPP_gen.h"

// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------


// Compute max = max_{b} |_ num_b_max/den_b_max + num_cst_max/(b*den_cst_max) _|
double aux_max_computation(int64 num_b_max, int64 den_b_max, int64 num_cst_max, int64 den_cst_max, optionMPP* option) {
	assert(den_b_max>0);
	assert(den_cst_max>0);
	
	double max = ( ((double) num_b_max) / ((double) den_b_max) );
	if (option->kMinMaxOption==0) {
		// => Compute a kMin/kMax for every possible values of "b".
		// max = max(n_b/d_b + n_cst/(b*d_cst))
		// d_cst is non-negative, so we have to look at the sign of n_cst
		if (num_cst_max>=0) {
			// b=minb is the best we can do
			max += ( ((double) num_cst_max) / ((double) (den_cst_max * option->minBlSizeParam) ) );
		} else {
			// b as large as possible is the best we can do
			
			// => if den_b_max>1 (already rational), 
			if (den_b_max>1) {
				// put the fraction at "-\epsilon". We majorate it by "0"
				max += 0.0;
			} else {
				// Because the "-\epsilon" makes a difference when the rest is integer
				// when taking the lower-bound, we have a special case here
				max -= 1.0;
			}
		}
	} else {
		assert(option->kMinMaxOption==1);
		// => Compute a tighter kMin/kMax by assuming that "b" is always "very large".
		if (num_cst_max>0) {
			max += 0.001;			// We have an "+\epsilon" here => majorate it by "0.001"
		} else {
			if (den_b_max>1)		// (cf option->kMinMaxOption==0)
				max += 0.0;
			else
				max -= 1.0;
		}
	}
	return max;
}

// Compute min = min_{b} |_ num_b_min/den_b_min + num_cst_min/(b*den_cst_min) _|
double aux_min_computation(int64 num_b_min, int64 den_b_min, int64 num_cst_min, int64 den_cst_min, optionMPP* option) {
	assert(den_b_min>0);
	assert(den_cst_min>0);
	
	double min = ( ((double) num_b_min) / ((double) den_b_min) );
	if (option->kMinMaxOption==0) {
		// => Compute a kMin/kMax for every possible values of "b".
		// min = min(n_b/d_b + n_cst/(b*d_cst))
		// d_cst is non-negative, so we have to look at the sign of n_cst
		if (num_cst_min>=0) {
			// b as large as possible is the best we can do
			//		=> put the fraction at "+\epsilon". We minorate it by "0"
			min += 0.0;
		} else {
			// b=minb is the best we can do
			min += ( ((double) num_cst_min) / ((double) (den_cst_min * option->minBlSizeParam) ) );
		}
	} else {
		assert(option->kMinMaxOption==1);
		// => Compute a tighter kMin/kMax by assuming that "b" is "very large".
		if (num_cst_min>0) {
			min += 0.0;		// (cf option->kMinMaxOption==0)
		} else {
			// We have an "-\epsilon" here => minorate it by "0.001"
			// If the rest is integral, then the "-\epsilon" is important => special case
			if (den_b_min>1) {
				min += -0.001;	// We have an "-\epsilon" here => minorate it by "0.001"
			} else {
				min += -1.0;
			}
		}
	}
	return min;
}



// Auxilliary function for the polyhedral case
// (\forall c) k_c(alpha, i_l, p_l) = |_ Q_c.L.D^{-1}.\alpha + \frac{Q_c.i_l + Q^{(p)}_c.p_l + q_c}{b} _|
void get_kmaxkmin_poly(long* kmax, long* kmin,
			int64** linPart, int64** paramPart, int64* constPart, int nConstr, int nInd, int nParam,
			int64* ineqEqPart_shape, int64** linPart_shape, int64** paramPart_shape, int64* constPart_shape, int nConstr_shape,
			int64** int_lattice, int64* den_lattice, optionMPP* option) {
	
	// Prelim: we build the polyhedron {i_l, p_l | i_l \in shape && 0\leq p_l \leq b.\vec{1}}
	int nRow_mat_ilpl = nConstr_shape+2*nParam+1;
	int nCol_mat_ilpl = 2+nInd+nParam+1;
	int64** mat_ilpl = (int64**) malloc(nRow_mat_ilpl * sizeof(int64*));
	for (int i=0; i<nRow_mat_ilpl; i++)
		mat_ilpl[i] = (int64*) malloc(nCol_mat_ilpl * sizeof(int64));
	for (int i=0; i<nRow_mat_ilpl; i++)
		for (int j=0; j<nCol_mat_ilpl; j++)
			mat_ilpl[i][j] = 0;
	
	// (First lines): shape*b
	for (int i=0; i<nConstr_shape; i++) {
		mat_ilpl[i][0] = ineqEqPart_shape[i];						// Inequality
		for (int j=0; j<nInd; j++)
			mat_ilpl[i][1+j] = linPart_shape[i][j];					// i_l
		mat_ilpl[i][1+nInd+nParam] = paramPart_shape[i][0];			// "b"
		mat_ilpl[i][1+nInd+nParam+1] = constPart_shape[i];			// Const
	}
	
	// (Second lines): 0\leq p_l
	for (int i=0; i<nParam; i++) {
		mat_ilpl[nConstr_shape+i][0] = 1;							// Inequality
		mat_ilpl[nConstr_shape+i][1+nInd+i] = 1;					// p_l
	}
	
	// (Third lines): 0\leq b.\vec{1] - p_l 
	for (int i=0; i<nParam; i++) {
		mat_ilpl[nConstr_shape+nParam+i][0] = 1;					// Inequality
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+i] = -1;			// p_l
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+nParam] = 1;		// "b"
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+nParam+1] = -1;		// Const
	}
	
	// (Fourth line): b_min \leq b
	mat_ilpl[nConstr_shape+2*nParam][0] = 1;
	mat_ilpl[nConstr_shape+2*nParam][1+nInd+nParam] = 1;
	mat_ilpl[nConstr_shape+2*nParam][1+nInd+nParam+1] = - option->minBlSizeParam;
	
	polyhedronMPP* constr_ilpl = buildPolyhedron(mat_ilpl, nRow_mat_ilpl, nInd+nParam, 1);
	
	// For all dimensions of k...
	for (int c=0; c<nConstr; c++) {
		// 1) Getting the maximal value of k_c
		// Note that an upper-bound is actually enough (it will generate empty polyhedra in the union if not tight enough)
		
		// First part: max_{i_l,p_l} \frac{Q_c.i_l + Q^{(p)}_c.p_l + q_c}{b} ?
		//      where i_l \in shape and \vec{0}<=p_l<b.\vec{1}
		int nRow_mat_objfunc = 1;
		int nCol_mat_objfunc = nInd+nParam+1+1;
		int64** mat_objfunc = (int64**) malloc(nRow_mat_objfunc * sizeof(int64*));
		for (int i=0; i<nRow_mat_objfunc; i++)
			mat_objfunc[i] = (int64*) malloc(nCol_mat_objfunc * sizeof(int64));
		
		for (int j=0; j<nInd; j++)
			mat_objfunc[0][j] = linPart[c][j];				// i_l
		for (int j=0; j<nParam; j++)
			mat_objfunc[0][nInd+j] = paramPart[c][j];		// p_l
		mat_objfunc[0][nCol_mat_objfunc-1] = constPart[c];	// const
		
		affFuncMPP* obj_func = buildAffineFunction(mat_objfunc, 1, nInd+nParam, 1);
		
		rational64** max_1 = max(obj_func, constr_ilpl);
		// max_1 is a 1*2 matrix (first col: b / second col: const)
		
		rational64 coeff_b_max1   = max_1[0][0];
		rational64 coeff_cst_max1 = max_1[0][1];
		int64 num_b_max1 = coeff_b_max1.num;
		int64 den_b_max1 = coeff_b_max1.den;
		int64 num_cst_max1 = coeff_cst_max1.num;
		int64 den_cst_max1 = coeff_cst_max1.den;
		
		// DEBUG
		//cout << "num_b_max1 = " << num_b_max1 << endl;
		//cout << "den_b_max1 = " << den_b_max1 << endl;
		//cout << "num_cst_max1 = " << num_cst_max1 << endl;
		//cout << "den_cst_max1 = " << den_cst_max1 << endl;
		
		// max1 = max_{b} |_ num_b_max1/den_b_max1 + num_cst_max1/(b*den_cst_max1) _|
		double max1 = aux_max_computation(num_b_max1, den_b_max1, num_cst_max1, den_cst_max1, option);
		
		
		// Second part: max_{\alpha} Q_c.L.D^{-1}.\alpha ?
		//	=> Done by looking at the sign of Q_c.L (note: \vec{0} \leq \alpha < D.\vec{1})
		int64* QcL = vectMatrixMultiplication(linPart[c], nInd, int_lattice, nInd, nInd); 
		double max2 = 0.0;
		for (int i=0; i<nInd; i++) {
			if (QcL[i]>=0) {
				max2 += ((double) QcL[i] * (den_lattice[i]-1)) / ((double) den_lattice[i]);
			} else {
				max2 += 0.0;  // Negative contribution
			}
		}
		
		// Third part: |_ Part 1 + Part 2 _|, and update of kmax
		double max3 = max1 + max2;
		long max_final = floor(max3);
		kmax[c] = max_final;
		
		// DEBUG
		//cout << "max1 = " << max1 << " | max2 = " << max2 << " | max3 = " << max3 << endl << endl;
		
		
		// 2) Getting the minimal value of k_c
		
		// First part: min_{i_l,p_l} \frac{Q_c.i_l + Q^{(p)}_c.p_l + q_c}{b} ?
		//      where i_l \in shape and \vec{0}<=p_l<b.\vec{1}
		rational64** min_1 = min(obj_func, constr_ilpl);
		rational64 coeff_b_min1   = min_1[0][0];
		rational64 coeff_cst_min1 = min_1[0][1];
		int64 num_b_min1 = coeff_b_min1.num;
		int64 den_b_min1 = coeff_b_min1.den;
		int64 num_cst_min1 = coeff_cst_min1.num;
		int64 den_cst_min1 = coeff_cst_min1.den;
		
		// DEBUG
		//cout << "num_b_min1 = " << num_b_min1 << endl;
		//cout << "den_b_min1 = " << den_b_min1 << endl;
		//cout << "num_cst_min1 = " << num_cst_min1 << endl;
		//cout << "den_cst_min1 = " << den_cst_min1 << endl;
		
		assert(den_b_min1>0);
		assert(den_cst_min1>0);
		
		double min1 = aux_min_computation(num_b_min1, den_b_min1, num_cst_min1, den_cst_min1, option);
		
		// Second part: min_{\alpha} Q_c.L.D^{-1}.\alpha ?
		//	=> Done by looking at the sign of Q_c.L (note: \vec{0} \leq \alpha < D.\vec{1})
		double min2 = 0.0;
		for (int i=0; i<nInd; i++) {
			if (QcL[i]>=0) {
				min2 += 0.0;	// Positive contribution
			} else {
				min2 += ((double) QcL[i] * (den_lattice[i]-1)) / ((double) den_lattice[i]);
			}
		}
		
		// Third part: |_ Part 1 + Part 2 _|, and update of kmin
		double min3 = min1 + min2;
		long min_final = floor(min3);
		kmin[c] = min_final;
		
		// DEBUG
		//cout << "min1 = " << min1 << " | min2 = " << min2 << " | min3 = " << min3 << endl << endl;
		
		// Free temporary structures
		freeAffineFunction(obj_func);
		free(QcL);
	}
	
	// Free temporary structures
	freePolyhedron(constr_ilpl);
	
	return;
}

// Polyhedral case
list<list<polyhedronMPP*> > getTiledDomain(polyhedronMPP *polyScalar, polyhedronMPP *shape, int64** lattice, optionMPP *option) {
	// Assertions
	assert(option->minBlSizeParam>=1);
	assert(polyScalar->nInd==shape->nInd);
	assert(shape->nParam==1);			// Only the block parameter
	
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
		matConst[0] = (int64*) malloc( (3+2*nParam) * sizeof(int64));
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
	
	
	
	// Extraction of the information from the polyhedron
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
	
	// Extraction of the information from the polyhedral shape
	polyhedronMPP* shape_noEq = eliminateEqualities(shape);
	int nConstr_shape = shape_noEq->nConstr;
	
	int64* ineqEqPart_shape = (int64*) malloc(nConstr_shape * sizeof(int64));
	int64** linPart_shape = (int64**) malloc(nConstr_shape * sizeof(int64*));
	for (int i=0; i<nConstr_shape; i++)
		linPart_shape[i] = (int64*) malloc(nInd * sizeof(int64));
	int64** paramPart_shape = (int64**) malloc(nConstr_shape * sizeof(int64*));
	for (int i=0; i<nConstr_shape; i++)
		paramPart_shape[i] = (int64*) malloc(1 * sizeof(int64));
	int64* constPart_shape = (int64*) malloc(nConstr_shape * sizeof(int64));
	
	extractPoly(shape_noEq, ineqEqPart_shape, linPart_shape, paramPart_shape, constPart_shape);
	
	// Extraction of the information from the lattice of tile origin
	int64** int_lattice = (int64**) malloc(nInd * sizeof(int64*));
	for (int i=0; i<nInd; i++)
		int_lattice[i] = (int64*) malloc(nInd * sizeof(int64));
	int64* den_lattice = (int64*) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++)
			int_lattice[i][j] = lattice[i][j];
	for (int j=0; j<nInd; j++)
		den_lattice[j] = lattice[nInd][j];
	
	// \delta = ppcm(den_lattice[..])
	int64 delta = ppcm_array(den_lattice, nInd);
	int64* delta_den = (int64*) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		delta_den[i] = delta / den_lattice[i];
	
	// LDiagdelta_den = int_lattice * Diag(delta_den)
	int64** LDiagdelta_den = (int64**) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		LDiagdelta_den[i] = (int64*) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++)
			LDiagdelta_den[i][j] = int_lattice[i][j] * delta_den[j];
	
	// Computation of kmin/kmax
	long* kmax = (long*) malloc(nConstr * sizeof(long));
	long* kmin = (long*) malloc(nConstr * sizeof(long));
	for (int i=0; i<nConstr; i++) {
		kmin[i] = 0;
		kmax[i] = 0;
	}
	
	get_kmaxkmin_poly(kmax, kmin, linPart, paramPart, constPart, nConstr, nInd, nParam,
	    ineqEqPart_shape, linPart_shape, paramPart_shape, constPart_shape, nConstr_shape,
	    int_lattice, den_lattice, option);
	
	/* DEBUG - kmin / kmax
	cout << "DEBUG: kmin - kmax" << endl;
	for (int i=0; i<nConstr; i++)
		cout << "kmin[" << i << "] = " << kmin[i] << "  |  kmax[" << i << "] = " << kmax[i] << endl;
	cout << endl;
	//*/
	
	// Note: we are iterating over the multi-dimensional bounding box [|\vec{kmin}; \vec{kmax}|]
	//		=> We can improve the algorithm by only iterating over the feasible \vec{k} instead
	//			of considering each one of its component independently.
	//		Benefits: less empty polyhedron to remove from the resulting union
	//			However, it complexifies the algorithm
	int64* alpha = (int64*) malloc(nInd * sizeof(int64));
	
	// We consider the ith constraint (outer intersection)
	for (int c=0; c<nConstr; c++) {
		// Iterating over the \vec{0}\leq \alpha < \vec{den_lattice}
		for (int i=0; i<nInd; i++)
			alpha[i] = 0;
		
		while (alpha[nInd-1]<den_lattice[nInd-1]) {
			list<polyhedronMPP*> domIthList;
			// * Case k_c = k_c^{min}
			//
			// In Polylib format, the matrix of constraint is:
			//
			// Row to know which column corresponds to what...
			// ( eq |          i_b          | i_l |    p_b    | p_l |   b    |                 const                       )
			//
			// [ 1  | Q_c.L.Diag(delta_den) |  0  |  \del.R_c |  0  |   0    | \del.k_c^min - Q_c.L.Diag(delta_den).\alpha ]
			// [shEq|           0           |shlin|     0     |  0  | shpar  |                sh_const                     ]
			// [ 1  |           0           |  0  |     0     |  0  |   1    |                 -minb                       ]
			//
			// (if option->areParamDiv = true)
			// [ 0  |           0           |  0  |     0     | Id  |   0    |                   0                         ]
			//
			// First column: 0 = equality / 1 = inequality
			// i_b = blocked index parameter / i_l = local index parameter
			// p_b = blocked parameters / p_l = local parameters / b = block size parameter
			//
			// Modulo constraint:
			// ( mod | i_b | i_l | p_b | p_l | b |  const   )
			//
			// [ den | Id  |  0  |  0  |  0  | 0 | -\alpha  ]
			int nRow_blConstr = (option->areParamDiv)? 2+nConstr_shape+nParam : 2+nConstr_shape;
			int nColumn_blConstr = 3+2*nParam+2*nInd;
			int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
			for (int l=0; l<nRow_blConstr; l++)
			 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
			for (int i=0; i<nRow_blConstr; i++)
				for (int j=0; j<nColumn_blConstr; j++)
					blockedConstr[i][j] = 0;
			
			// (First line)
			int64* QcLDiagdelta_den = vectMatrixMultiplication(linPart[c], nInd, LDiagdelta_den, nInd, nInd);
			
			blockedConstr[0][0] = 1;
			for (int j=0; j<nInd; j++)
				blockedConstr[0][1+j] = QcLDiagdelta_den[j];
			for (int j=0; j<nParam; j++)
				blockedConstr[0][1+2*nInd+j] = delta * paramPart[c][j];
			int64 tempScalarProd = 0;
			for (int k=0; k<nInd; k++)
				tempScalarProd += QcLDiagdelta_den[k] * alpha[k];
			blockedConstr[0][nColumn_blConstr-1] = delta * kmin[c] - tempScalarProd;
			
			// (Second lines)
			for (int i=0; i<nConstr_shape; i++) {
				blockedConstr[1+i][0] = ineqEqPart_shape[i];
				for (int j=0; j<nInd; j++)
					blockedConstr[1+i][1+nInd+j] = linPart_shape[i][j];
				blockedConstr[1+i][1+2*nInd+2*nParam] = paramPart_shape[i][0];  // Single column
				blockedConstr[1+i][1+2*nInd+2*nParam+1] = constPart_shape[i];
			}
			
			// (Third line)
			blockedConstr[1+nConstr_shape][0] = 1;
			blockedConstr[1+nConstr_shape][1+2*nInd+2*nParam] = 1;
			blockedConstr[1+nConstr_shape][nColumn_blConstr-1] = - option->minBlSizeParam;
			
			// (Fourth line)
			if (option->areParamDiv)
				for (int i=0; i<nParam; i++)
					blockedConstr[2+nConstr_shape+i][1+2*nInd+nParam+i] = 1;
			
			
			
			// Modulo
			int nRow_modConstr = nInd;
			int nColumn_modConstr = 3+2*nParam+2*nInd;
			int64** modBlockedConstr = (int64**) malloc(nRow_modConstr * sizeof(int64*));
			for (int l=0; l<nRow_modConstr; l++)
			 	modBlockedConstr[l] = (int64*) malloc(nColumn_modConstr * sizeof(int64));
			for (int i=0; i<nRow_modConstr; i++)
				for (int j=0; j<nColumn_modConstr; j++)
					modBlockedConstr[i][j] = 0;
			
			// (First line)
			for (int i=0; i<nInd; i++) {
				modBlockedConstr[i][0] = den_lattice[i];
				modBlockedConstr[i][1+i] = 1;
				modBlockedConstr[i][nColumn_modConstr-1] = -alpha[i];
			}
			
			polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
			polyRet->nConstrMod = nInd;
			polyRet->modConstr = modBlockedConstr;
			domIthList.push_back(polyRet);
			
			
			// * Case k_c^{min} < k_c <= k_c{max}
			for (int kc = kmin[c]+1; kc<=kmax[c]; kc++) {
				// In Polylib format, the matrix of constraint is:
				//
				// Row to know which column corresponds to what...
				// ( eq |          i_b          |  i_l   |  p_b   |   p_l  |   b    |             const                       )
				//
				// [ 0  | Q_c.L.Diag(delta_den) |   0    |\del.R_c|    0   |   0    | \del.k_c - Q_c.L.Diag(delta_den).\alpha ]
				// [ 1  |           0           |\del.Q_c|   0    |\del.R_c|   A    |          \delta.q_c                     ]
				// [shEq|           0           |  shlin |   0    |    0   |  shpar |            sh_const                     ]
				// [ 1  |           0           |   0    |   0    |    0   |   1    |             -minb                       ]
				//
				// where A = Q_c.L.Diag(delta_den).\alpha - \del.k_c
				//
				// (if option->areParamDiv = true)
				// [ 0  |           0           |   0    |   0    |   Id   |   0    |               0                         ]
				//
				// First column: 0 = equality / 1 = inequality
				// i_b = blocked index parameter / i_l = local index parameter
				// p_b = blocked parameters / p_l = local parameters / b = block size parameter
				//
				// Modulo constraint:
				// ( mod | i_b | i_l | p_b | p_l | b |  const   )
				//
				// [ den | Id  |  0  |  0  |  0  | 0 | -\alpha  ]
				
				// Examine the gcd of i_b compared to the gcd of the rest of the first constraint (which is an equality)
				int64 gcd_ib = gcd_array(QcLDiagdelta_den, nInd);
				int64 gcd_pb = (nParam==0)? delta : delta * gcd_array(paramPart[c], nParam);
				int64 gcd_ibpb = gcd(gcd_ib, gcd_pb);
				
				int64 const_val_eq = delta * kc - tempScalarProd;
				int64 mod = (const_val_eq % gcd_ibpb);
				
				bool equality_unsatisfiable = !(mod==0);
				if (equality_unsatisfiable)
					continue;
				
				// Building the matrix of constraints
				int nRow_blConstr = (option->areParamDiv)? 3+nConstr_shape+nParam : 3+nConstr_shape;
				int nColumn_blConstr = 3+2*nParam+2*nInd;
				int64** blockedConstr = (int64**) malloc(nRow_blConstr * sizeof(int64*));
				for (int l=0; l<nRow_blConstr; l++)
				 	blockedConstr[l] = (int64*) malloc(nColumn_blConstr * sizeof(int64));
				for (int i=0; i<nRow_blConstr; i++)
					for (int j=0; j<nColumn_blConstr; j++)
						blockedConstr[i][j] = 0;
				
				// (First line)
				// Note: QcLDiagdelta_den already computed for the case kc=kmin
				// tempScalarProd also already computed
				for (int j=0; j<nInd; j++)
					blockedConstr[0][1+j] = QcLDiagdelta_den[j];
				for (int j=0; j<nParam; j++)
					blockedConstr[0][1+2*nInd+j] = delta * paramPart[c][j];
				blockedConstr[0][nColumn_blConstr-1] = delta * kc - tempScalarProd;
				
				// (Second line)
				blockedConstr[1][0] = 1;
				for (int j=0; j<nInd; j++)
					blockedConstr[1][1+nInd+j] = delta * linPart[c][j];
				for (int j=0; j<nParam; j++)
					blockedConstr[1][1+2*nInd+nParam+j] = delta * paramPart[c][j];
				blockedConstr[1][1+2*nInd+2*nParam] = tempScalarProd - delta * kc;
				blockedConstr[1][nColumn_blConstr-1] = delta * constPart[c];
				
				// (Third line)
				for (int i=0; i<nConstr_shape; i++) {
					blockedConstr[2+i][0] = ineqEqPart_shape[i];
					for (int j=0; j<nInd; j++)
						blockedConstr[2+i][1+nInd+j] = linPart_shape[i][j];
					blockedConstr[2+i][1+2*nInd+2*nParam] = paramPart_shape[i][0];
					blockedConstr[2+i][1+2*nInd+2*nParam+1] = constPart_shape[i];
				}
				
				// (Fourth line)
				blockedConstr[2+nConstr_shape][0] = 1;
				blockedConstr[2+nConstr_shape][1+2*nInd+2*nParam] = 1;
				blockedConstr[2+nConstr_shape][nColumn_blConstr-1] = - option->minBlSizeParam;
				
				// (Fifth line)
				if (option->areParamDiv)
					for (int i=0; i<nParam; i++)
						blockedConstr[2+nConstr_shape+i][1+2*nInd+nParam+i] = 1;
				
				// Modulo
				
				int nRow_modConstr = nInd;
				int nColumn_modConstr = 3+2*nParam+2*nInd;
				int64** modBlockedConstr = (int64**) malloc(nRow_modConstr * sizeof(int64*));
				for (int l=0; l<nRow_modConstr; l++)
				 	modBlockedConstr[l] = (int64*) malloc(nColumn_modConstr * sizeof(int64));
				for (int i=0; i<nRow_modConstr; i++)
					for (int j=0; j<nColumn_modConstr; j++)
						modBlockedConstr[i][j] = 0;
				
				// (First line)
				for (int i=0; i<nInd; i++) {
					modBlockedConstr[i][0] = den_lattice[i];
					modBlockedConstr[i][1+i] = 1;
					modBlockedConstr[i][nColumn_modConstr-1] = -alpha[i];
				}
			
				polyhedronMPP* polyRet = buildPolyhedron(blockedConstr, nRow_blConstr, 2*nInd, 2*nParam+1);
				polyRet->nConstrMod = nInd;
				polyRet->modConstr = modBlockedConstr;
				domIthList.push_back(polyRet);
			}
			free(QcLDiagdelta_den);
			
			// We increase the multi-dimensional iterator, starting from the first dimension and propagating the overflows
			alpha[0]++;
			for (int i=0; i<nInd-1; i++)
				if (alpha[i]>=den_lattice[i]) {
					alpha[i] = 0;
					alpha[i+1]++;
				}
			resDom.push_back(domIthList);
		}
	}
	
	// Free temporary structures
	freePolyhedron(polyScalar_noEq);
	free(ineqEqPart);
	freeMatrix(linPart, nConstr);
	freeMatrix(paramPart, nConstr);
	free(constPart);
	
	freePolyhedron(shape_noEq);
	free(ineqEqPart_shape);
	freeMatrix(linPart_shape, nConstr_shape);
	freeMatrix(paramPart_shape, nConstr_shape);
	free(constPart_shape);
	
	freeMatrix(int_lattice, nInd);
	free(den_lattice);
	free(delta_den);
	freeMatrix(LDiagdelta_den, nInd);
	
	free(kmax);
	free(kmin);
	free(alpha);
	
	return resDom;
}


// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------

// Auxilliary function for the affine function case
void get_kmaxkmin_func(long* kmax, long* kmin, long* kImmax, long* kImmin,
		int64** linPart, int64** paramPart, int64* constPart, int nDimOut, int nInd, int nParam,
		int64* ineqEqPart_shape, int64** linPart_shape, int64** paramPart_shape, int64* constPart_shape, int nConstr_shape,
		int64* ineqEqPart_shapeIm, int64** linPart_shapeIm, int64** paramPart_shapeIm, int64* constPart_shapeIm, int nConstr_shapeIm,
		int64** int_lattice, int64* den_lattice, int64 delta,
		int64** int_latticeIm, int64* den_latticeIm, int64 deltaIm,
		optionMPP* option) {
	
	// PART 1: kIm(\alpha', i'_l) = |_ L'.D'^{-1} . \alpha' + i'_l/b _|
	//				where 0\leq L'.D^{-1}.\alpha'< \delta'.\vec{1} && i'_l \in shapeIm
	
	// We pre-build the polyhedron { i'_l | i'_l \in shapeIm }, b being the only parameter
	int nCol_constr1Im = 3 + nDimOut;
	int64** matConstr1Im = (int64**) malloc(nConstr_shapeIm * sizeof(int64*));
	for (int i=0; i<nConstr_shapeIm ; i++)
		matConstr1Im[i] = (int64*) malloc(nCol_constr1Im * sizeof(int64));
	
	for (int i=0; i<nConstr_shapeIm; i++) {
		matConstr1Im[i][0] = ineqEqPart_shapeIm[i];				// Ineq/eq
		for (int j=0; j<nDimOut; j++)
			matConstr1Im[i][1+j] = linPart_shapeIm[i][j];		// Lin
		matConstr1Im[i][1+nDimOut] = paramPart_shapeIm[i][0];	// b
		matConstr1Im[i][2+nDimOut] = constPart_shapeIm[i];		// Const
	}
	polyhedronMPP* constr1Im = buildPolyhedron(matConstr1Im, nConstr_shapeIm, nDimOut, 1);
	
	
	// First part of kIm: i'_l/b where i'_l \in constr1Im
	// For the objective function, we get the min/max of each components of i'_l separately
	//		by examining each of its dimension separately
	double* kIm_min_1 = (double*) malloc(nDimOut * sizeof(double));
	double* kIm_max_1 = (double*) malloc(nDimOut * sizeof(double));
	
	for (int c=0; c<nDimOut; c++) {
		// Objective function for the c-th dimension: (\vec{i} \mapsto \vec{e_c}.\vec{i})
		int nRow_mat_objfunc = 1;
		int nCol_mat_objfunc = nDimOut+1+1;
		int64** mat_objfunc = (int64**) malloc(nRow_mat_objfunc * sizeof(int64*));
		for (int i=0; i<nRow_mat_objfunc; i++)
			mat_objfunc[i] = (int64*) malloc(nCol_mat_objfunc * sizeof(int64));
		for (int j=0; j<nCol_mat_objfunc; j++)
			mat_objfunc[0][j] = 0;
		
		mat_objfunc[0][c] = 1;							// Lin
		affFuncMPP* obj_func = buildAffineFunction(mat_objfunc, 1, nInd, 1);
		
		rational64** maxIm_1 = max(obj_func, constr1Im);
		
		rational64 coeff_b_maxIm1   = maxIm_1[0][0];
		rational64 coeff_cst_maxIm1 = maxIm_1[0][1];
		int64 num_b_maxIm1 = coeff_b_maxIm1.num;
		int64 den_b_maxIm1 = coeff_b_maxIm1.den;
		int64 num_cst_maxIm1 = coeff_cst_maxIm1.num;
		int64 den_cst_maxIm1 = coeff_cst_maxIm1.den;
		
		kIm_max_1[c] = aux_max_computation(num_b_maxIm1, den_b_maxIm1, num_cst_maxIm1, den_cst_maxIm1, option);
		
		
		rational64** minIm_1 = min(obj_func, constr1Im);
		
		rational64 coeff_b_minIm1   = minIm_1[0][0];
		rational64 coeff_cst_minIm1 = minIm_1[0][1];
		int64 num_b_minIm1 = coeff_b_minIm1.num;
		int64 den_b_minIm1 = coeff_b_minIm1.den;
		int64 num_cst_minIm1 = coeff_cst_minIm1.num;
		int64 den_cst_minIm1 = coeff_cst_minIm1.den;
		
		kIm_min_1[c] = aux_min_computation(num_b_minIm1, den_b_minIm1, num_cst_minIm1, den_cst_minIm1, option);
		
		// Free temporary structures
		freeAffineFunction(obj_func);
	}
	
	// Second part: L'.D'^{-1}.\alpha' where 0\leq L'.D^{-1}.\alpha'< \delta'.\vec{1} 
	// Because of the definition of \alpha', kIm_min_2 = 0 and kIm_max_2 = delta'
	double* kIm_min_2 = (double*) malloc(nDimOut * sizeof(double));
	double* kIm_max_2 = (double*) malloc(nDimOut * sizeof(double));
	
	for (int c=0; c<nDimOut; c++) {
		kIm_min_2[c] = 0.0;
		kIm_max_2[c] = (double) deltaIm;
	}
	
	// Final part: fill the values of kImmax and kImmin
	for (int c=0; c<nDimOut; c++) {
		long max_final = floor(kIm_max_1[c] + kIm_max_2[c]);
		kImmax[c] = max_final;
		
		long min_final = floor(kIm_min_1[c] + kIm_min_2[c]);
		kImmin[c] = min_final;
	}
	
	// Free temporary structures
	freePolyhedron(constr1Im);
	free(kIm_min_1);
	free(kIm_max_1);
	free(kIm_min_2);
	free(kIm_max_2);
	
	
	// ------------------------------------------------------------------------
	// PART 2: k = |_ Q.L.D^{-1} . \alpha + (Q.i_l+R.p_l.q)/b _|
	
	// Prelim: we build the polyhedron {i_l, p_l | i_l \in shape && 0\leq p_l \leq b.\vec{1}}
	int nRow_mat_ilpl = nConstr_shape+2*nParam+1;
	int nCol_mat_ilpl = 2+nInd+nParam+1;
	int64** mat_ilpl = (int64**) malloc(nRow_mat_ilpl * sizeof(int64*));
	for (int i=0; i<nRow_mat_ilpl; i++)
		mat_ilpl[i] = (int64*) malloc(nCol_mat_ilpl * sizeof(int64));
	for (int i=0; i<nRow_mat_ilpl; i++)
		for (int j=0; j<nCol_mat_ilpl; j++)
			mat_ilpl[i][j] = 0;
	
	// (First lines): shape*b
	for (int i=0; i<nConstr_shape; i++) {
		mat_ilpl[i][0] = ineqEqPart_shape[i];						// Inequality
		for (int j=0; j<nInd; j++)
			mat_ilpl[i][1+j] = linPart_shape[i][j];					// i_l
		mat_ilpl[i][1+nInd+nParam] = paramPart_shape[i][0];			// "b"
		mat_ilpl[i][1+nInd+nParam+1] = constPart_shape[i];			// Const
	}
	
	// (Second lines): 0\leq p_l
	for (int i=0; i<nParam; i++) {
		mat_ilpl[nConstr_shape+i][0] = 1;							// Inequality
		mat_ilpl[nConstr_shape+i][1+nInd+i] = 1;					// p_l
	}
	
	// (Third lines): 0\leq b.\vec{1] - p_l 
	for (int i=0; i<nParam; i++) {
		mat_ilpl[nConstr_shape+nParam+i][0] = 1;					// Inequality
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+i] = -1;			// p_l
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+nParam] = 1;		// "b"
		mat_ilpl[nConstr_shape+nParam+i][1+nInd+nParam+1] = -1;		// Const
	}
	
	// (Fourth line): b_min \leq b
	mat_ilpl[nConstr_shape+2*nParam][0] = 1;
	mat_ilpl[nConstr_shape+2*nParam][1+nInd+nParam] = 1;
	mat_ilpl[nConstr_shape+2*nParam][1+nInd+nParam+1] = - option->minBlSizeParam;
	
	polyhedronMPP* constr_ilpl = buildPolyhedron(mat_ilpl, nRow_mat_ilpl, nInd+nParam, 1);
	
	
	// First part of k: (Q.i_l+R.p_l.q)/b
	double* k_min_1 = (double*) malloc(nDimOut * sizeof(double));
	double* k_max_1 = (double*) malloc(nDimOut * sizeof(double));
	
	// For all dimensions of k...
	for (int c=0; c<nDimOut; c++) {
		// Building the corresponding objective function
		int nRow_mat_objfunc = 1;
		int nCol_mat_objfunc = nInd+nParam+1+1;
		int64** mat_objfunc = (int64**) malloc(nRow_mat_objfunc * sizeof(int64*));
		for (int i=0; i<nRow_mat_objfunc; i++)
			mat_objfunc[i] = (int64*) malloc(nCol_mat_objfunc * sizeof(int64));
		
		for (int j=0; j<nInd; j++)
			mat_objfunc[0][j] = linPart[c][j];				// i_l
		for (int j=0; j<nParam; j++)
			mat_objfunc[0][nInd+j] = paramPart[c][j];		// p_l
		mat_objfunc[0][nCol_mat_objfunc-1] = constPart[c];	// const
		
		affFuncMPP* obj_func = buildAffineFunction(mat_objfunc, 1, nInd+nParam, 1);
		
		rational64** max_1 = max(obj_func, constr_ilpl);
		// max_1 is a 1*2 matrix (first col: b / second col: const)
		
		rational64 coeff_b_max1   = max_1[0][0];
		rational64 coeff_cst_max1 = max_1[0][1];
		int64 num_b_max1 = coeff_b_max1.num;
		int64 den_b_max1 = coeff_b_max1.den;
		int64 num_cst_max1 = coeff_cst_max1.num;
		int64 den_cst_max1 = coeff_cst_max1.den;
		
		double max1 = aux_max_computation(num_b_max1, den_b_max1, num_cst_max1, den_cst_max1, option);
		k_max_1[c] = max1;
		
		rational64** min_1 = min(obj_func, constr_ilpl);
		
		rational64 coeff_b_min1   = min_1[0][0];
		rational64 coeff_cst_min1 = min_1[0][1];
		int64 num_b_min1 = coeff_b_min1.num;
		int64 den_b_min1 = coeff_b_min1.den;
		int64 num_cst_min1 = coeff_cst_min1.num;
		int64 den_cst_min1 = coeff_cst_min1.den;
		
		double min1 = aux_min_computation(num_b_min1, den_b_min1, num_cst_min1, den_cst_min1, option);
		k_min_1[c] = min1;
		
		// Free temporary structures
		freeAffineFunction(obj_func);
	}
	
	
	// Second part: Q. (L.D^{-1}.\alpha)  where 0\leq L.D^{-1}.\alpha < \delta.vec{1}
	//		=> we just have to find the max and min of Q.\lambda, where 0\leq \lambda < \delta.vec{1}
	double* k_min_2 = (double*) malloc(nDimOut * sizeof(double));
	double* k_max_2 = (double*) malloc(nDimOut * sizeof(double));
	
	for (int c=0; c<nDimOut; c++) {
		double max2 = 0.0;
		for (int i=0; i<nInd; i++) {
			if (linPart[i]>=0) {
				max2 += (double) linPart[c][i] * (delta-1);
			} else {
				max2 += 0.0;  // Negative contribution
			}
		}
		k_max_2[c] = max2;
	}
	
	for (int c=0; c<nDimOut; c++) {
		double min2 = 0.0;
		for (int i=0; i<nInd; i++) {
			if (linPart[i]>=0) {
				min2 += 0.0;	// Positive contribution
			} else {
				min2 += (double) linPart[c][i] * (delta-1);
			}
		}
		k_min_2[c] = min2;
	}
	
	// Final part: fill the values of kImmax and kImmin
	for (int c=0; c<nDimOut; c++) {
		long max_final = floor(k_max_1[c] + k_max_2[c]);
		kmax[c] = max_final;
		
		long min_final = floor(k_min_1[c] + k_min_2[c]);
		kmin[c] = min_final;
	}
	
	// Free temporary structures
	free(k_min_1);
	free(k_max_1);
	free(k_min_2);
	free(k_max_2);
	
	return;
}


// Affine function case
map<polyhedronMPP*, affFuncMPP*> getTiledFunction(affFuncMPP *affScalar,
			polyhedronMPP *shape, int64** lattice,
			polyhedronMPP *shapeIm, int64** latticeIm,
			optionMPP *option) {
	
	// Assertions
	assert(option->minBlSizeParam>=1);
	assert(affScalar->nInd==shape->nInd);
	assert(affScalar->dimOut==shapeIm->nInd);
	assert(shape->nParam==1);			// Only the block parameter
	assert(shapeIm->nParam==1);			// Only the block parameter
	
	int nParam = affScalar->nParam;
	int nInd = affScalar->nInd;
	int nDimOut = affScalar->dimOut;
	
	// DEBUG
	//cout << "nInd = " << nInd << " | nParam = " << nParam << " | nDimOut = " << nDimOut << endl;
	
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
	int64** linPart = (int64**) malloc(nDimOut * sizeof(int64*));
	for (int i=0; i<nDimOut; i++)
		linPart[i] = (int64*) malloc(nInd * sizeof(int64));
	int64** paramPart = (int64**) malloc(nDimOut * sizeof(int64*));
	for (int i=0; i<nDimOut; i++)
		paramPart[i] = (int64*) malloc(nParam * sizeof(int64));
	int64* constPart = (int64*) malloc(nDimOut * sizeof(int64));
	extractAffFunc(affScalar, linPart, paramPart, constPart);
	
	// Extraction of the information from the polyhedral shape for the input space
	polyhedronMPP* shape_noEq = eliminateEqualities(shape);
	int nConstr_shape = shape_noEq->nConstr;
	
	int64* ineqEqPart_shape = (int64*) malloc(nConstr_shape * sizeof(int64));
	int64** linPart_shape = (int64**) malloc(nConstr_shape * sizeof(int64*));
	for (int i=0; i<nConstr_shape; i++)
		linPart_shape[i] = (int64*) malloc(nInd * sizeof(int64));
	int64** paramPart_shape = (int64**) malloc(nConstr_shape * sizeof(int64*));
	for (int i=0; i<nConstr_shape; i++)
		paramPart_shape[i] = (int64*) malloc(1 * sizeof(int64));
	int64* constPart_shape = (int64*) malloc(nConstr_shape * sizeof(int64));
	extractPoly(shape_noEq, ineqEqPart_shape, linPart_shape, paramPart_shape, constPart_shape);
	
	// Extraction of the information from the lattice of tile origin for the input space
	int64** int_lattice = (int64**) malloc(nInd * sizeof(int64*));
	for (int i=0; i<nInd; i++)
		int_lattice[i] = (int64*) malloc(nInd * sizeof(int64));
	int64* den_lattice = (int64*) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++)
			int_lattice[i][j] = lattice[i][j];
	for (int j=0; j<nInd; j++)
		den_lattice[j] = lattice[nInd][j];
	
	
	// Extraction of the information from the polyhedral shape for the output space
	polyhedronMPP* shapeIm_noEq = eliminateEqualities(shapeIm);
	int nConstr_shapeIm = shapeIm_noEq->nConstr;
	
	int64* ineqEqPart_shapeIm = (int64*) malloc(nConstr_shapeIm * sizeof(int64));
	int64** linPart_shapeIm = (int64**) malloc(nConstr_shapeIm * sizeof(int64*));
	for (int i=0; i<nConstr_shapeIm; i++)
		linPart_shapeIm[i] = (int64*) malloc(nDimOut * sizeof(int64));
	int64** paramPart_shapeIm = (int64**) malloc(nConstr_shapeIm * sizeof(int64*));
	for (int i=0; i<nConstr_shapeIm; i++)
		paramPart_shapeIm[i] = (int64*) malloc(1 * sizeof(int64));
	int64* constPart_shapeIm = (int64*) malloc(nConstr_shapeIm * sizeof(int64));
	extractPoly(shapeIm_noEq, ineqEqPart_shapeIm, linPart_shapeIm, paramPart_shapeIm, constPart_shapeIm);
	
	// Extraction of the information from the lattice of tile origin for the output space
	int64** int_latticeIm = (int64**) malloc(nDimOut * sizeof(int64*));
	for (int i=0; i<nDimOut; i++)
		int_latticeIm[i] = (int64*) malloc(nDimOut * sizeof(int64));
	int64* den_latticeIm = (int64*) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++)
			int_latticeIm[i][j] = latticeIm[i][j];
	for (int j=0; j<nDimOut; j++)
		den_latticeIm[j] = latticeIm[nInd][j];
	
	
	// Pre-computation of D.L^{-1} and D'.L'^{-1}
	rational64** invL = (rational64**) malloc(nInd * sizeof(rational64*));
	for (int i=0; i<nInd; i++)
		invL[i] = (rational64*) malloc(nInd * sizeof(rational64));
	inverseDet(int_lattice, invL, nInd); // Does not use the determinant of L here
	
	rational64** DLinv = (rational64**) malloc(nInd * sizeof(rational64*));
	for (int i=0; i<nInd; i++)
		DLinv[i] = (rational64*) malloc(nInd * sizeof(rational64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++) {
			rational64 temp; temp.num = den_lattice[i]; temp.den = 1;
			DLinv[i][j] = multiplyRational(temp, invL[i][j]);
		}
	freeMatrix(invL, nInd);
	
	rational64** invLIm = (rational64**) malloc(nDimOut * sizeof(rational64*));
	for (int i=0; i<nDimOut; i++)
		invLIm[i] = (rational64*) malloc(nDimOut * sizeof(rational64));
	inverseDet(int_latticeIm, invLIm, nDimOut); // Does not use the determinant of LIm here
	
	rational64** DLinvIm = (rational64**) malloc(nDimOut * sizeof(rational64*));
	for (int i=0; i<nDimOut; i++)
		DLinvIm[i] = (rational64*) malloc(nDimOut * sizeof(rational64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++) {
			rational64 temp; temp.num = den_latticeIm[i]; temp.den = 1;
			DLinvIm[i][j] = multiplyRational(temp, invLIm[i][j]);
		}
	
	rational64** LinvDIm = (rational64**) malloc(nDimOut * sizeof(rational64*));
	for (int i=0; i<nDimOut; i++)
		LinvDIm[i] = (rational64*) malloc(nDimOut * sizeof(rational64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++) {
			rational64 temp; temp.num = den_latticeIm[j]; temp.den=1;
			LinvDIm[i][j] = multiplyRational(invLIm[i][j], temp);
		}
	
	freeMatrix(invLIm, nDimOut);
	
	
	// Computation of \delta and \delta', such that $\delta.DLinv$ and $\delta'.DLinvIm$ are integral
	int64* DLinvDen = (int64*) malloc(nInd * nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++)
			DLinvDen[i*nInd+j] = DLinv[i][j].den;
	int64 delta = ppcm_array(DLinvDen, nInd * nInd);
	free(DLinvDen);
	
	int64* DLinvDenIm = (int64*) malloc(nDimOut * nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++)
			DLinvDenIm[i*nDimOut+j] = DLinvIm[i][j].den;
	int64 deltaIm = ppcm_array(DLinvDenIm, nDimOut * nDimOut);
	free(DLinvDenIm);
	
	// Computation of \delta.D.L^{-1}
	int64** deltaDLinv = (int64**) malloc(nInd * sizeof(int64*));
	for (int i=0; i<nInd; i++)
		deltaDLinv[i] = (int64*) malloc(nInd * sizeof(int64));
	for (int i=0; i<nInd; i++)
		for (int j=0; j<nInd; j++)
			deltaDLinv[i][j] = (delta / DLinv[i][j].den) * DLinv[i][j].num;
	
	// Computation of \delta'.D'.L'^{-1}
	int64** deltaDLinvIm = (int64**) malloc(nDimOut * sizeof(int64*));
	for (int i=0; i<nDimOut; i++)
		deltaDLinvIm[i] = (int64*) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++)
			deltaDLinvIm[i][j] = (deltaIm / DLinvIm[i][j].den) * DLinvIm[i][j].num;
	
	
	// Computation of kmin/kmax and kImmin/kImmax
	long* kmax = (long*) malloc(nDimOut * sizeof(long));
	long* kmin = (long*) malloc(nDimOut * sizeof(long));
	for (int i=0; i<nDimOut; i++) {
		kmin[i] = 0;
		kmax[i] = 0;
	}
	
	long* kImmax = (long*) malloc(nDimOut * sizeof(long));
	long* kImmin = (long*) malloc(nDimOut * sizeof(long));
	for (int i=0; i<nDimOut; i++) {
		kImmax[i] = 0;
		kImmin[i] = 0;
	}
	
	get_kmaxkmin_func(kmax, kmin, kImmax, kImmin,
	    linPart, paramPart, constPart, nDimOut, nInd, nParam,
	    ineqEqPart_shape, linPart_shape, paramPart_shape, constPart_shape, nConstr_shape,
	    ineqEqPart_shapeIm, linPart_shapeIm, paramPart_shapeIm, constPart_shapeIm, nConstr_shapeIm,
	    int_lattice, den_lattice, delta,
	    int_latticeIm, den_latticeIm, deltaIm,
	    option);
	
	/* DEBUG - kmin / kmax
	cout << "DEBUG: kmin - kmax" << endl;
	for (int i=0; i<nDimOut; i++)
		cout << "kmin[" << i << "] = " << kmin[i] << "  |  kmax[" << i << "] = " << kmax[i] << endl;
	cout << endl;
	//*/
	
	// Note: we are iterating over the multi-dimensional bounding box [|\vec{kmin}; \vec{kmax}|]
	//		=> We can improve the algorithm by only iterating over the feasible \vec{k} instead
	//			of considering each one of its component independently.
	//		Benefits: less empty polyhedron to remove from the resulting union
	//			However, it complexifies the algorithm
	
	int64* kCurr = (int64*) malloc(nDimOut * sizeof(int64));
	int64* alpha = (int64*) malloc(nInd * sizeof(int64));
	int64* kImCurr = (int64*) malloc(nDimOut * sizeof(int64));
	int64* alphaIm = (int64*) malloc(nDimOut * sizeof(int64));
	
	// Level 1 of iteration: kCurr
	for (int i=0; i<nDimOut; i++)
		kCurr[i] = kmin[i];
	while (kCurr[nDimOut-1]<=kmax[nDimOut-1]) {
	
	// Level 2 of iteration: kImCurr
	for (int i=0; i<nDimOut; i++)
		kImCurr[i] = kImmin[i];
	while (kImCurr[nDimOut-1]<=kImmax[nDimOut-1]) {
	
	// Level 3 of iterator: alpha, where \vec{0} \leq L.D^{-1}.\vec{\alpha} < \delta.\vec{1}
	for (int i=0; i<nInd; i++)
		alpha[i] = ;		// TODO
	while () {				// TODO
	
	for (int i=0; i<nDimOut; i++)
		alphaIm[i] = ;		// TODO
	while() {				// TODO
		
		// We build the branch of the piecewise affine function corresponding to (kCurr, kImCurr, alpha, alphaIm)
		
		// Precomputation of the following products:
		//		- LDinv = L.D^{-1}
		//		- QLDinv = Q.L.D^{-1}
		//		- DLinvImQLDinv = D'.L'^{-1}.Q.L.D^{-1}
		rational64** LDinv = (rational64**) malloc(nInd * sizeof(rational64*));
		for (int i=0; i<nInd; i++)
			LDinv[i] = (rational64*) malloc(nInd * sizeof(rational64));
		for (int i=0; i<nInd; i++)
			for (int j=0; j<nInd; j++) {
				rational64 temp; temp.num = int_lattice[i][j]; temp.den = den_lattice[i];
				LDinv[i][j] = simplify(temp);
			}
		rational64** Qrat = toRationalMatrix(linPart, nDimOut, nInd);
		rational64** QLDinv = matrixMultiplication(Qrat, nDimOut, nInd, LDinv, nInd, nInd);
		rational64** DLinvImQLDinv = matrixMultiplication(DLinvIm, nDimOut, nDimOut, QLDinv, nDimOut, nInd);
		
		freeMatrix(Qrat, nDimOut);
		
		// Precomputation of LDinvIm = L'.D'^{-1}
		rational64** LDinvIm = (rational64**) malloc(nDimOut * sizeof(rational64*));
		for (int i=0; i<nDimOut; i++)
			LDinvIm[i] = (rational64*) malloc(nDimOut * sizeof(rational64));
		for (int i=0; i<nDimOut; i++)
			for (int j=0; j<nDimOut; j++) {
				rational64 temp; temp.num = int_latticeIm[i][j]; temp.den = den_latticeIm[i];
				LDinvIm[i][j] = simplify(temp);
			}
		
		rational64* ratAlpha = toRationalVector(alpha, nInd);
		
		// Precomputation of LDinv.\alpha for the second row
		rational64* LDinvalpha = matVectProduct(LDinv, nInd, nInd, ratAlpha, nInd);
		
		
		
		// Note: \vec{i_b} = \vec{\alpha} mod A <=> (\exist \lambda) \vec{i_b} = \vec{\alpha} + A.\lambda
		//		If A^{-1} = 1/det . Ainv where "Ainv" is integral and "det" an integer, then
		//			<=> (\exist \lambda) Ainv.\vec{i_b} = Ainv.\vec{\alpha} + det.\lambda
		//			<=> Ainv.\vec{i_b} = Ainv.\vec{\alpha} mod (det.\vec{1})
		
		
		// In Polylib format, the matrix of input constraints is:
		// 
		// ( eq/ineq| i_b |  i_l  |  p_b  |    p_l    |     b     |  const )    <= Row to know which column corresponds to what...
		// 
		// [ shEq   |  0  |   shlin   | 0 |     0     |   shpar   |       sh_const     ]		// i_l \in shape
		// [ shEqIm |  0  | shlinIm.Q | 0 | shlinIm.R |     X     | shconstIm+shlinIm.q]		// i'_l \in shapeIm
		// [   1    |  0  |     Q     | 0 |     R     | Q.L.D^{-1}.\alpha-k   |   q    ]
		// [   1    |  0  |    -Q     | 0 |    -R     | k+1-Q.L.D^{-1}.\alpha|  -q     ]
		// [  detM  | Minv |    0     | 0 |     0     |     0     |   -Minv.\alpha     ]
		// [ detM'  | Minv'|    0     | 0 |     0     |     0     |  -Minv'.\alpha'    ]
		// where:
		//		* X = shparIm + shlinIm(QLDinv.\alpha - LinvDIm.\alpha'+k'-k)
		//		* A = \delta.D.L^{-1}
		//		* det = det(A)  and Auninv = det * A^{-1}
		//      * M = deltaDLinv and M^{-1} = 1/detM * Minv
		//      * M' = \delta'.D'.L'^{-1} and M'^{-1} = 1/detM' * Minv'
		//
		// Note:
		// \alpha = i_b mod (\delta . D. L^{-1})      <=> Minv.\alpha - Minv.i_b = 0 mod detM
		// \alpha' = i'_b mod (\delta' . D'. L'^{-1}) <=> Minv'.\alpha' - Minv'.i_b' = 0 mod detM'
		//		
		// First column: 0=equality / 1 = inequality
		// \alpha = blocked index parameter / ii = local index parameter
		// \rho = blocked parameters / pp = local parameters / b = block size parameter
		
		int nRow_blConstr = nConstr_shape + nConstr_shapeIm + 2*nDimOut + nInd + nDimOut;
		int nCol_blConstr = 2*nInd+2*nParam+3;
		
		int64** inputConstrLongMat = (int64**) malloc(nRow_blConstr * sizeof(int64*));
		for (int l=0; l<nRow_blConstr; l++)
		 	inputConstrLongMat[l] = (int64*) malloc(nCol_blConstr * sizeof(int64));
		for (int i=0; i<nRow_blConstr; i++)
				for (int j=0; j<nCol_blConstr; j++)
					inputConstrLongMat[i][j] = 0;
		
		// (First rows)
		for (int i=0; i<nConstr_shape; i++) {
			inputConstrLongMat[i][0] = ineqEqPart_shape[i];							// Eq/Ineq
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[i][1+nInd+j] = linPart_shape[i][j];				// i_l
			inputConstrLongMat[i][1+2*nInd+2*nParam] = paramPart_shape[i][0];		// b
			inputConstrLongMat[i][nCol_blConstr-1] = constPart_shape[i];			// Const
		}
		
		// Precomputation for the second row
		int64** shlinImQ = matrixMultiplication(linPart_shapeIm, nConstr_shapeIm, nDimOut, linPart, nDimOut, nInd);
		int64** shlinImR = matrixMultiplication(linPart_shapeIm, nConstr_shapeIm, nDimOut, paramPart, nDimOut, nParam);
		
		// tempVectb = QLDinv.\alpha - LinvDIm.\alpha'+k'-k
		rational64* tempVectb = (rational64*) malloc(nDimOut * sizeof(rational64));
		for (int i=0; i<nDimOut; i++) {
			rational64 temp; temp.num = kImCurr[i] - kCurr[i]; temp.den = 1;
			tempVectb[i] = temp;
		}
		rational64* QLDinvalpha = matVectProduct(QLDinv, nDimOut, nInd, ratAlpha, nInd);
		rational64* LinvDImAlphaIm = matVectProduct(LinvDIm, nDimOut, nDimOut, alphaIm, nInd);
		for (int i=0; i<nDimOut; i++)
			tempVectb[i] = addRational(tempVectb[i], subRational(QLDinvalpha[i], LinvDImAlphaIm[i]));
		free(LinvDImAlphaIm);
		
		// (Second rows)
		for (int i=0; i<nConstr_shapeIm; i++) {
			inputConstrLongMat[nConstr_shape+i][0] = ineqEqPart_shapeIm[i];						// Eq/Ineq
			
			// Because the contribution for "b" might be a rational, we start computing it first, then
			//		multiply all the coeffs by its denominator
			// b : shparIm + shlinIm.tempVectb
			rational64* ratLinShapeLine = toRationalVector(linPart_shapeIm[i], nDimOut);
			rational64 tempbRat = dotProduct(ratLinShapeLine, tempVectb, nDimOut);
			free(ratLinShapeLine);
			
			rational64 shparImRat; shparImRat.num = paramPart_shapeIm[i]; shparImRat.den = 1;
			rational64 tempRatb = addRational(tempbRat, shparImRat);
			
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[nConstr_shape+i][1+nInd+j] = tempRatb.den * shlinImQ[i][j];			// i_l
			for (int j=0; j<nParam; j++)
				inputConstrLongMat[nConstr_shape+i][1+2*nInd+nParam+j] = tempRatb.den *shlinImR[i][j];	// p_l
			
			inputConstrLongMat[nConstr_shape+i][1+2*nInd+2*nParam] = tempRatb.num;						// b
			
			int64 shlinImq = dotProduct(constPart_shape, nDimOut, constPart, nDimOut);					// Const
			inputConstrLongMat[nConstr_shape+i][nCol_blConstr-1] = (constPart_shapeIm + shlinImq) * tempRatb.den;
		}
		
		freeMatrix(shlinImQ, nConstr_shapeIm);
		freeMatrix(shlinImR, nConstr_shapeIm);
		free(tempVectb);
		
		
		
		// (Third rows)
		for (int i=0; i<nDimOut; i++) {
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+i][0] = 1;														// Eq/Ineq
			
			// Because the contribution for "b" might be a rational, we start computing it first, then
			//		multiply all the coeffs by its denominator
			// b: Q.L.D^{-1}.\alpha-k
			rational64 tempRatb = subRational(QLDinvalpha[i], toRational(kCurr[i]));
			
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[nConstr_shape+nConstr_shapeIm+i][1+nInd+j] = linPart[i][j] * tempRatb.den;				// i_l
			for (int j=0; j<nParam; j++)
				inputConstrLongMat[nConstr_shape+nConstr_shapeIm+i][1+2*nInd+nParam+j] = paramPart[i][j] * tempRatb.den;	// p_l
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+i][1+2*nInd+2*nParam] = tempRatb.num;							// b
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+i][nCol_blConstr-1] = constPart[i] * tempRatb.den;				// Const
		}
		
		// (Fourth rows)
		for (int i=0; i<nDimOut; i++) {
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+nDimOut+i][0] = 1;													// Eq/Ineq
			
			// Because the contribution for "b" might be a rational, we start computing it first, then
			//		multiply all the coeffs by its denominator
			// b: k+1-Q.L.D^{-1}.\alpha
			rational64 tempRatb = subRational(toRational(kCurr[i]+1), QLDinvalpha[i]);
			
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[nConstr_shape+nConstr_shapeIm+nDimOut+i][1+nInd+i] = -linPart[i][j] * tempRatb.den;			// i_l
			for (int j=0; j<nParam; j++)
				inputConstrLongMat[nConstr_shape+nConstr_shapeIm+nDimOut+i][1+2*nInd+nParam+i] = -paramPart[i][j]*tempRatb.den;	// p_l
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+nDimOut+i][1+2*nInd+2*nParam] = tempRatb.num;						// b
			inputConstrLongMat[nConstr_shape+nConstr_shapeIm+nDimOut+i][nCol_blConstr-1] = - constPart[i] * tempRatb.den;		// Const
		}
		
		free(QLDinvalpha);
		
		
		// Precomputation for the fifth row
		rational64** MinvDivdetM = (rational64**) malloc(nInd * sizeof(rational64*));
		for (int i=0; i<nInd; i++)
			MinvDivdetM[i] = (rational64*) malloc(nInd * sizeof(rational64));
		int64 detM = inverseDet(deltaDLinv, MinvDivdetM, nInd);
		
		int64** Minv = (int64**) malloc(nInd * sizeof(int64*));
		for (int i=0; i<nInd; i++)
			Minv[i] = (int64*) malloc(nInd * sizeof(int64));
		for (int i=0; i<nInd; i++)
			for (int j=0; j<nInd; j++) {
				rational64 temp = MinvDivdetM[i][j];
				temp = multiplyRational(temp, toRational(detM));
				assert(temp.den==1);				// By construction of detM and MinvDivdetM
				Minv[i][j] = temp.num;
			}
		freeMatrix(MinvDivdetM, nInd);
		
		
		// (Fifth rows)
		int offset_5r = nConstr_shape + nConstr_shapeIm + 2*nDimOut;
		for (int i=0; i<nInd; i++) {
			inputConstrLongMat[offset_5r+i][0] = detM;								// Eq/Ineq
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[offset_5r+i][1+j] = Minv[i][j];					// i_b
			
			int64 temp = dotProduct(Minv[i], nInd, alpha, nInd);
			inputConstrLongMat[offset_5r+i][nCol_blConstr-1] = -temp;				// Const
		}
		
		freeMatrix(Minv, nInd);
		
		
		// Precomputation for the sixth row
		rational64** MinvDivdetMIm = (rational64**) malloc(nDimOut * sizeof(rational64*));
		for (int i=0; i<nDimOut; i++)
			MinvDivdetMIm[i] = (rational64*) malloc(nDimOut * sizeof(rational64));
		int64 detMIm = inverseDet(deltaDLinvIm, MinvDivdetMIm, nDimOut);
		
		int64** MinvIm = (int64**) malloc(nDimOut * sizeof(int64*));
		for (int i=0; i<nDimOut; i++)
			MinvIm[i] = (int64*) malloc(nDimOut * sizeof(int64));
		for (int i=0; i<nDimOut; i++)
			for (int j=0; j<nDimOut; j++) {
				rational64 temp = MinvDivdetMIm[i][j];
				temp = multiplyRational(temp, toRational(detMIm));
				assert(temp.den==1);				// By construction of detMIm and MinvDivdetMIm
				MinvIm[i][j] = temp.num;
			}
		freeMatrix(MinvDivdetMIm, nDimOut);
		
		// (Sixth rows)
		for (int i=0; i<nDimOut; i++) {
			inputConstrLongMat[offset_5r+nInd+i][0] = detMIm;						// Eq/Ineq
			for (int j=0; j<nInd; j++)
				inputConstrLongMat[offset_5r+nInd+i][1+j] = MinvIm[i][j];			// i_b
			int64 temp = dotProduct(MinvIm[i], nInd, alphaIm, nInd);
			inputConstrLongMat[offset_5r+nInd+i][nCol_blConstr-1] = -temp;			// Const
		}
		
		freeMatrix(MinvIm, nDimOut);
		
		
		polyhedronMPP* polyRet = buildPolyhedron(inputConstrLongMat, nRow_blConstr, 2*nInd, 2*nParam+1);
		
		
		
		// The matrix of the affine function is:
		// 
		// (    \alpha     |ii |   \rho    |pp |     b         | const)    <= Row to know which column corresponds to what...
		// [ DLinvImQLDinv | 0 | DLinvIm.R | 0 |  0  | \alpha'-DLinvImQLDinv.\alpha + DLinvIm.(k - k') ]
		// [       0       | Q |     0     | R | Q.LDinv.\alpha - LDinvIm.\alpha'+k'-k    |     q      ]
		//	where:
		//		- DLinvIm = D'.L'^{-1}
		//		- LDinv = L.D^{-1}
		//		- DLinvImQLDinv = DLinvIm.Q.LDinv
		//		- LDinvIm = L'.D'^{-1}
		
		int nRow_relConstr = 2*nDimOut;
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
		// [ DLinvImQLDinv | 0 | DLinvIm.R | 0 |  0  | \alpha'-DLinvImQLDinv.\alpha + DLinvIm.(k - k') ]
		for (int i=0; i<nDimOut; i++) {
			// Doing the computations on the temporary corresponding rational row
			rational64* tempRatRow = (rational64*) malloc(nColumn_relConstr * sizeof(rational64));
			for (int j=0; j<nInd; j++)						// \alpha
				tempRatRow[j] = DLinvImQLDinv[i][j];
			
			for (int j=0; j<nParam; j++) {					// \rho
				// Computing in place DLinvIm*R
				rational64 temp = toRational(0);
				for (int k=0; k<nDimOut; k++)
					temp = addRational(temp, multiplyRational(DLinvIm[i][k], toRational(paramPart[k][j])));
				tempRatRow[2*nInd+j] = temp;
			}
			
			// constElem = \alpha' - DLinvImQLDinv.\alpha + DLinvIm.(k - k') 
			rational64 constElem = toRational(alphaIm[i]);
			
			rational64 elemConstAlpha = dotProduct(DLinvImQLDinv[i], ratAlpha, nInd);
			constElem = subRational(constElem, elemConstAlpha);
			
			int64* kminuskIm = (int64*) malloc(nDimOut * sizeof(int64*));
			for (int i=0; i<nDimOut; i++)
				kminuskIm[i] = kCurr[i] - kImCurr[i];
			rational64* ratKminuskIm = toRationalVector(kminuskIm, nDimOut);
			rational64 elemConstkCurr = dotProduct(DLinvIm[i], ratKminuskIm, nDimOut);
			
			constElem = addRational(constElem, elemConstkCurr);
			tempRatRow[nColumn_relConstr-1] = constElem;
			
			
			// Divisor management (aka, going back to int64)
			for (int j=0; j<nColumn_relConstr; j++)
				tempRatRow[j] = simplify(tempRatRow[j]);
			
			int64 commonDiv = 1;
			for (int j=0; j<nColumn_relConstr; j++)
				commonDiv = ppcm(commonDiv, tempRatRow[j].den);
			
			divConstrLongMat[i] = commonDiv;
			for (int j=0; j<nColumn_relConstr; j++)
				relationConstrLongMat[i][j] = (commonDiv/tempRatRow[j].den) * tempRatRow[j].num;
			// Guarranty that the division on the last line is integral
			
			free(tempRatRow);
			free(kminuskIm);
			free(ratKminuskIm);
		}
		
		
		
		// (Second rows)
		// [       0       | Q |     0     | R | Q.LDinv.\alpha - LDinvIm.\alpha'+k'-k    |     q      ]
		for (int i=0; i<nDimOut; i++) {
			// Doing the computations on the temporary corresponding rational row
			rational64* tempRatRow = (rational64*) malloc(nColumn_relConstr * sizeof(rational64));
		
			for (int j=0; j<nInd; j++)
				tempRatRow[nInd+j] = toRational(linPart[i][j]);
			for (int j=0; j<nParam; j++)
				tempRatRow[2*nInd+nParam+j] = toRational(paramPart[i][j]);
			
			// b : Q.LDinv.\alpha - LDinvIm.\alpha' + k'-k
			rational64 ratbElem = toRational(kImCurr[i] - kCurr[i]);
			
			rational64* ratLinPartRow = toRationalVector(linPart[i], nInd);
			rational64 ratbElemAlpha = dotProduct(ratLinPartRow, LDinvalpha, nInd);
			ratbElem = addRational(ratbElem, ratbElemAlpha);
			
			rational64* ratAlphaIm = toRationalVector(alphaIm, nDimOut);
			rational64 ratbElemAlphaIm = dotProduct(LDinvIm[i], ratAlphaIm, nDimOut);
			ratbElem = subRational(ratbElem, ratbElemAlphaIm);
			
			tempRatRow[2*nInd+2*nParam] = ratbElem;
			
			tempRatRow[nColumn_relConstr-1] = toRational(constPart[i]);
			
			
			// Divisor management (aka, going back to int64)
			for (int j=0; j<nColumn_relConstr; j++)
				tempRatRow[j] = simplify(tempRatRow[j]);
			
			int64 commonDiv = 1;
			for (int j=0; j<nColumn_relConstr; j++)
				commonDiv = ppcm(commonDiv, tempRatRow[j].den);
			
			divConstrLongMat[i] = commonDiv;
			for (int j=0; j<nColumn_relConstr; j++)
				relationConstrLongMat[nDimOut+i][j] = (commonDiv/tempRatRow[j].den) * tempRatRow[j].num;
			// Guarranty that the division on the last line is integral
			
			free(ratLinPartRow);
			free(ratAlphaIm);
			free(tempRatRow);
		}
		
		
		affFuncMPP* affFuncRet = buildRationalAffineFunction(relationConstrLongMat, divConstrLongMat, 2*nDimOut, 2*nInd, 2*nParam+1);
		simplifyAffFuncMPP(affFuncRet);
		
		result.insert ( pair<polyhedronMPP*, affFuncMPP*>(polyRet, affFuncRet) );
		
		// Free temporary data structures
		freeMatrix(LDinv, nInd);
		freeMatrix(QLDinv, nDimOut);
		freeMatrix(DLinvImQLDinv, nDimOut);
		freeMatrix(LDinvIm, nDimOut);
		free(LDinvalpha);
		free(ratAlpha);
		// End of construction of the branch
		
		
		// TODO: alphaIm++
		
		
	} // End of while loop on alphaIm
	
	// TODO alpha++
	
	
	} // End of while loop on alpha
	
	// "kImCurr++", with overflow propagation
	kImCurr[0]++;
	for (int i=0; i<nDimOut-1; i++)
		if (kImCurr[i]>kImmax[i]) {
			kImCurr[i] = kImmin[i];
			kImCurr[i+1]++;
		}
	} // End of while loop on kImCurr
	
	// "kCurr++", with overflow propagation
	kCurr[0]++;
	for (int i=0; i<nDimOut-1; i++)
		if (kCurr[i]>kmax[i]) {
			kCurr[i] = kmin[i];
			kCurr[i+1]++;
		}
	} //  // End of while loop on kCurr
	
	
	// Free temporary data structures
	freeMatrix(linPart, nDimOut);
	freeMatrix(paramPart, nDimOut);
	free(constPart);
	
	freePolyhedron(shape_noEq);
	free(ineqEqPart_shape);
	freeMatrix(linPart_shape, nConstr_shape);
	freeMatrix(paramPart_shape, nConstr_shape);
	free(constPart_shape);
	
	freeMatrix(int_lattice, nInd);
	free(den_lattice);
	freeMatrix(DLinv, nInd);
	freeMatrix(deltaDLinv, nInd);
	
	freePolyhedron(shapeIm_noEq);
	free(ineqEqPart_shapeIm);
	freeMatrix(linPart_shapeIm, nConstr_shapeIm);
	freeMatrix(paramPart_shapeIm, nConstr_shapeIm);
	free(constPart_shapeIm);
	
	freeMatrix(int_latticeIm, nDimOut);
	free(den_latticeIm);
	freeMatrix(DLinvIm, nDimOut);
	freeMatrix(deltaDLinvIm, nDimOut);
	freeMatrix(LinvDIm, nDimOut);
	
	free(kmax);
	free(kmin);
	free(kImmax);
	free(kImmin);
	
	free(kCurr);
	free(kImCurr);
	free(alpha);
	free(alphaIm);
	
	return result;
}


