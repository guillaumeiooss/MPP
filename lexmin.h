/* Interface to a lexmin function (PIP, ISL).
   We took the same signature than Pluto (cf "constraints.h"
      => "pluto_constraints_lexmin(const PlutoConstraints *cst, int negvar) ).
   
   The implementation was strongly inspired by some methods from Pluto's "contraints.c"
   source file */
#include "linAlg.h"
#include "include/piplib/piplib64.h"

#ifndef __LEXMIN_MPP_H
#define __LEXMIN_MPP_H

struct rational64 {
	int64 num;
	int64 den;
};

// Pretty-printer
void printMatrix(rational64** mat, int nRow, int nCol);


// Lexmax functions (using Piplib behind)
rational64** lexmax(polyhedronMPP *poly);
rational64** lexmax(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context);

rational64** lexmin(polyhedronMPP *poly);
rational64** lexmin(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context);


#endif

