// Standalone C++ version of Monoparametric tiling
// Author: Guillaume Iooss

#include "MPP_gen.h"

// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------

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
		
		assert(den_b_max1>0);
		assert(den_cst_max1>0);
		
		// max1 = max_{b} |_ num_b_max1/den_b_max1 + num_cst_max1/(b*den_cst_max1) _|
		double max1 = ( ((double) num_b_max1) / ((double) den_b_max1) );
		if (option->kMinMaxOption==0) {
			// => Compute a kMin/kMax for every possible values of "b".
			// max1 = max(n_b/d_b + n_cst/(b*d_cst))
			// d_cst is non-negative, so we have to look at the sign of n_cst
			if (num_cst_max1>=0) {
				// b=minb is the best we can do
				max1 += ( ((double) num_cst_max1) / ((double) (den_cst_max1 * option->minBlSizeParam) ) );
			} else {
				// b as large as possible is the best we can do
				
				// => if den_b_max1>1 (already rational), 
				if (den_b_max1>1) {
					// put the fraction at "-\epsilon". We majorate it by "0"
					max1 += 0.0;
				} else {
					// Because the "-\epsilon" makes a difference when the rest is integer
					// when taking the lower-bound, we have a special case here
					max1 -= 1.0;
				}
			}
		} else {
			assert(option->kMinMaxOption==1);
			// => Compute a tighter kMin/kMax by assuming that "b" is always "very large".
			if (num_cst_max1>0) {
				max1 += 0.001;			// We have an "+\epsilon" here => majorate it by "0.001"
			} else {
				if (den_b_max1>1)		// (cf option->kMinMaxOption==0)
					max1 += 0.0;
				else
					max1 -= 1.0;
			}
		}
		
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
		
		double min1 = ( ((double) num_b_min1) / ((double) den_b_min1) );
		if (option->kMinMaxOption==0) {
			// => Compute a kMin/kMax for every possible values of "b".
			// min1 = min(n_b/d_b + n_cst/(b*d_cst))
			// d_cst is non-negative, so we have to look at the sign of n_cst
			if (num_cst_min1>=0) {
				// b as large as possible is the best we can do
				//		=> put the fraction at "+\epsilon". We minorate it by "0"
				min1 += 0.0;
			} else {
				// b=minb is the best we can do
				min1 += ( ((double) num_cst_min1) / ((double) (den_cst_min1 * option->minBlSizeParam) ) );
			}
		} else {
			assert(option->kMinMaxOption==1);
			// => Compute a tighter kMin/kMax by assuming that "b" is "very large".
			if (num_cst_min1>0) {
				min1 += 0.0;		// (cf option->kMinMaxOption==0)
			} else {
				// We have an "-\epsilon" here => minorate it by "0.001"
				// If the rest is integral, then the "-\epsilon" is important => special case
				if (den_b_min1>1) {
					min1 += -0.001;	// We have an "-\epsilon" here => minorate it by "0.001"
				} else {
					min1 += -1.0;
				}
			}
		}
		
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
	}
	
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
	if (polyScalar->nInd==0) {
		int64** matConst = (int64**) malloc(1 * sizeof(int64*));
		matConst[0] = (int64*) malloc( (3+2*polyScalar->nParam) * sizeof(int64));
		polyhedronMPP* nullPoly = buildPolyhedron(matConst, 1, 0, 1+2*nParam);
		
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
void get_kmaxkmin_func(long* kmax, long* kmin,
		int64** linPart, int64** paramPart, int64* constPart, int nDimOut, int nInd, int nParam,
		int64* ineqEqPart_shape, int64** linPart_shape, int64** paramPart_shape, int64* constPart_shape, int nConstr_shape,
		int64* ineqEqPart_shapeIm, int64** linPart_shapeIm, int64** paramPart_shapeIm, int64* constPart_shapeIm, int nConstr_shapeIm,
		int64** int_lattice, int64* den_lattice, int64** int_latticeIm, int64* den_latticeIm, optionMPP* option) {
	
	// TODO
	
	return;
}



// Affine function case
map<polyhedronMPP*, affFuncMPP*> getRectangularTiledFunction(affFuncMPP *affScalar,
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
		int_lattice[i] = (int64*) malloc(nDimOut * sizeof(int64));
	int64* den_latticeIm = (int64*) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++)
			int_latticeIm[i][j] = latticeIm[i][j];
	for (int j=0; j<nDimOut; j++)
		den_latticeIm[j] = latticeIm[nInd][j];
	
	
	// \deltaIm = ppcm(den_latticeIm[..])
	int64 deltaIm = ppcm_array(den_latticeIm, nDimOut);
	int64* delta_denIm = (int64*) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		delta_denIm[i] = deltaIm / den_latticeIm[i];
	
	// LDiagdelta_denIm = int_latticeIm * Diag(delta_denIm)
	int64** LDiagdelta_denIm = (int64**) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		LDiagdelta_denIm[i] = (int64*) malloc(nDimOut * sizeof(int64));
	for (int i=0; i<nDimOut; i++)
		for (int j=0; j<nDimOut; j++)
			LDiagdelta_denIm[i][j] = int_latticeIm[i][j] * delta_denIm[j];
	
	
	// Computation of kmin/kmax
	long* kmax = (long*) malloc(nDimOut * sizeof(long));
	long* kmin = (long*) malloc(nDimOut * sizeof(long));
	for (int i=0; i<nDimOut; i++) {
		kmin[i] = 0;
		kmax[i] = 0;
	}
	
	get_kmaxkmin_func(kmax, kmin, linPart, paramPart, constPart, nDimOut, nInd, nParam,
	    ineqEqPart_shape, linPart_shape, paramPart_shape, constPart_shape, nConstr_shape,
	    ineqEqPart_shapeIm, linPart_shapeIm, paramPart_shapeIm, constPart_shapeIm, nConstr_shapeIm,
	    int_lattice, den_lattice, int_latticeIm, den_latticeIm, option);
	
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
	
	
	
	
	// TODO: iteration over the k + construction of the matrices
	
	
	
	
	
	
	
	
	
	
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
	free(delta_den);
	freeMatrix(LDiagdelta_den, nInd);
	
	freePolyhedron(shapeIm_noEq);
	free(ineqEqPart_shapeIm);
	freeMatrix(linPart_shapeIm, nConstr_shapeIm);
	freeMatrix(paramPart_shapeIm, nConstr_shapeIm);
	free(constPart_shapeIm);
	
	freeMatrix(int_latticeIm, nDimOut);
	free(den_latticeIm);
	free(delta_denIm);
	freeMatrix(LDiagdelta_denIm, nDimOut);
	
	free(kmax);
	free(kmin);
	
	// TODO: add iterator from the main loop
	
	return result;
}


