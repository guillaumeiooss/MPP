
#include "lexmin.h"

// Pretty-printer
void printMatrix(rational64** mat, int nRow, int nCol) {
	for (int i=0; i<nRow; i++) {
		cout << "[ ";
		for (int j=0; j<nCol; j++) {
			cout << mat[i][j].num << "/" << mat[i][j].den << " ";
		}
		cout << "]" << endl;
	}
	cout << endl;
	
	return;
}

/* ------------------------------------------------------------------- */

PipVector* copy_pip_vector(PipVector* v) {
	Entier* nnum = (Entier*) malloc(sizeof(Entier) * v->nb_elements);
	Entier* nden = (Entier*) malloc(sizeof(Entier) * v->nb_elements);
	
	for (int i=0; i<v->nb_elements; i++) {
		nnum[i] = v->the_vector[i];
		nden[i] = v->the_deno[i];
	}
	
	PipVector* nv = new PipVector();
	nv->nb_elements = v->nb_elements;
	nv->the_vector = nnum;
	nv->the_deno = nden;
	
	return nv;
}


PipList* copy_pip_list(PipList* l) {
	if (l==NULL)
		return NULL;
	
	PipList* nl = new PipList();
	nl->vector = copy_pip_vector(l->vector);
	nl->next = copy_pip_list(l->next);
	
	return nl;
}


// Sub-function of "merge_sols"
PipVector* merge_row_sols(PipVector* row_then, PipVector* row_else, bool bmax) {
	assert(row_then->nb_elements = row_else->nb_elements);
	
	// In our context, we are only looking at *integral* solutions
	//	=> We can assume than the denominator is always a vector of 0s (apparently encoded like that)
	Entier* num_then = row_then->the_vector;
	Entier* den_then = row_then->the_deno;
	Entier* num_else = row_else->the_vector;
	Entier* den_else = row_else->the_deno;
	
	int nb_elem = row_then->nb_elements;
	
	Entier* num_merge = (Entier*) malloc(sizeof(Entier) * nb_elem);
	Entier* den_merge = (Entier*) malloc(sizeof(Entier) * nb_elem);
	
	for (int i=0; i<nb_elem; i++) {
		float val_then = ((float) num_then[i]) / ((float) den_then[i]);
		float val_else = ((float) num_else[i]) / ((float) den_else[i]);
		
		// This works because the parameters of the polyhedron
		//      are inside the positive quadrant
		if (bmax) {
			if (val_then >= val_else) {
				num_merge[i] = num_then[i];
				den_merge[i] = den_then[i];
			} else {
				num_merge[i] = num_else[i];
				den_merge[i] = den_else[i];
			}
		} else {
			if (val_then >= val_else) {
				num_merge[i] = num_else[i];
				den_merge[i] = den_else[i];
			} else {
				num_merge[i] = num_then[i];
				den_merge[i] = den_then[i];
			}
		}
	}
	
	// Building the solution
	PipVector* row_merge = new PipVector();
	row_merge->nb_elements = nb_elem;
	row_merge->the_vector = num_merge;
	row_merge->the_deno = den_merge;
	
	return row_merge;
}


// Merge the solution of two branches of a quast
PipList* merge_sols(PipList* sol_then, PipList* sol_else, bool bmax) {
	assert((sol_then==NULL && sol_else==NULL) || (sol_then!=NULL && sol_else!=NULL));
	
	// Termination condition
	if (sol_then==NULL && sol_else==NULL) {
		return (PipList*) NULL;
	}
	
	// Recursion
	PipList* next_then = sol_then->next;
	PipList* next_else = sol_else->next;
	PipList* merge_next = merge_sols(next_then, next_else, bmax);
	
	// Merging current row
	PipVector* row_then = sol_then->vector;
	PipVector* row_else = sol_else->vector;
	PipVector* row_merge = merge_row_sols(row_then, row_else, bmax);
	
	// Building final result
	PipList* ret_list = new PipList();
	ret_list->vector = row_merge;
	ret_list->next = merge_next;
	
	return ret_list;
}



// Reduce all branches into a single one which is an upper-bound of all of them
PipList* flatten_quast(PipQuast* quast, bool bmax) {
	PipList* sol;
	if (quast->condition==NULL) {
		// We are on a leaf of the quast
		sol = copy_pip_list(quast->list);
	} else {
		// Recursions on next_then and next_else
		PipList* sol_then = flatten_quast(quast->next_then, bmax);
		PipList* sol_else = flatten_quast(quast->next_else, bmax);
		
		// We merge the two branches
		sol = merge_sols(sol_then, sol_else, bmax);
		
		// Free the temp structures
		pip_list_free(sol_then);
		pip_list_free(sol_else);
	}
	return sol;
}

/* ------------------------------------------------------- */

// Using piplib to find the lexmin of a polyhedron, while flattening 
// We assume that the parameters of the polyhedron are living on the positive quadrant
//	(ie, all parameter are positive)
// Initial inspiration: Pluto code which uses pip (l1077 - "constraints.c")
rational64** lexminmax(polyhedronMPP *poly, PipMatrix* context, bool bmax) {
	assert(poly->nParam+2 == (int) context->NbColumns);
	
	// 1) Building pipMat
	// Inspiration: "pluto_constraints_to_pip_matrix" / "pip_matrix_populate"
	// Note: the matrix is exactly the same than the one inside poly,
	//		except that there is no notion of program parameter
	int nRows = poly->nConstr;
	int nCols = poly->nInd + poly->nParam + 2;
	
	// Transform this matrix into a PipMatrix
	PipMatrix* pipMat = pip_matrix_alloc(nRows, nCols);
	Entier* p = pipMat->p_Init;
    for (unsigned int i=0; i<pipMat->NbRows; i++)  {
        for (unsigned int j=0; j<pipMat->NbColumns;j++)   {
            *(p++) = poly->polyScalar[i][j];
        }
    }
	
	// 2) Call to the core function
	int bignum = -1;					// No big parameter
	
	PipOptions* pipOptions = pip_options_init();
	pipOptions->Nq          = 0;		// Integral solutions
	pipOptions->Simplify    = 1;		// Simplify solution (needed?)
	pipOptions->Maximize    = bmax?1:0;	// Looking for maximum
	pipOptions->Urs_parms   =-1;		// No constr on coeff of params
	pipOptions->Urs_unknowns=-1;		// No constr on coeff of indices
	
	PipQuast* quastSol = pip_solve(pipMat, context, bignum, pipOptions);
	
	// DEBUG
	//cout << "QUAST = " << endl; pip_quast_print(stdout, quastSol, 0);
	
	// 3) Manage the different branches by merging them into a single one
	PipList* listSol = flatten_quast(quastSol, bmax);
	
	// DEBUG
	//cout << "LIST = " << endl; pip_list_print(stdout, listSol, 0);
	
	
	// 4) Convert back into the matrix of coefficients on the parameters
	int nrow_sol = poly->nInd;
	int ncol_sol = poly->nParam+1;
	
	// DEBUG
	//cout << "ncol_sol = " << ncol_sol << endl;
	//cout << "listSol->vector->nb_elements = " << listSol->vector->nb_elements << endl;
	
	// Check that there is no existential parameters
	if (nrow_sol>0)
		assert(listSol->vector->nb_elements==ncol_sol);
	
	rational64** sol = (rational64**) malloc(nrow_sol * sizeof(rational64*));
	for (int i=0; i<nrow_sol; i++)
		sol[i] = (rational64*) malloc(ncol_sol * sizeof(rational64));
	
	PipList* itSol = listSol;
	for (int i=0; i<nrow_sol; i++) {
		if (i>0)
			itSol = itSol->next;
		
		for (int j=0; j<ncol_sol; j++) {
			sol[i][j].num = itSol->vector->the_vector[j];
			if (itSol->vector->the_deno[j]==0)
				sol[i][j].den = 1;
			else
				sol[i][j].den = itSol->vector->the_deno[j];
		}
	}
	assert(itSol->next==NULL);
	
    // Free temporary structures
    pip_matrix_free(pipMat);
    pip_matrix_free(context);
    pip_options_free(pipOptions);
    pip_quast_free(quastSol);
    pip_list_free(listSol);
	
    return sol;
}

/* ------------------------------------------------------- */
// Option to put the context as a parameter

rational64** lexminmax(polyhedronMPP *poly, bool bmax) {
	int nParam = poly->nParam;
	
	// Context: all parameters must be positive
	// (else, contraction of quast will not work)
	PipMatrix* context = pip_matrix_alloc(nParam, nParam+2);
	for (int i=0; i<nParam; i++) {
		context->p[i][0] = 1;
		context->p[i][1+i] = 1;
	}
	
	return lexminmax(poly, context, bmax);
}


rational64** lexminmax(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context, bool bmax) {
	PipMatrix* pip_context = pip_matrix_alloc(nrow_context, ncol_context);
	for (int i=0; i<nrow_context; i++) {
		for (int j=0; j<ncol_context; j++)
			pip_context->p[i][j] = context[i][j];
	}
	
	return lexminmax(poly, pip_context, bmax);
}

/* ------------------------------------------------------- */
// Interface functions

rational64** lexmax(polyhedronMPP *poly) {
	return lexminmax(poly, true);
}

rational64** lexmin(polyhedronMPP *poly) {
	return lexminmax(poly, false);
}


rational64** lexmax(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context) {
	return lexminmax(poly, context, nrow_context, ncol_context, true);
}

rational64** lexmin(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context) {
	return lexminmax(poly, context, nrow_context, ncol_context, false);
}


