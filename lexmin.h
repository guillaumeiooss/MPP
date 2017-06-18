/* Interface to a lexmin function (PIP, ISL).
   We took the same signature than Pluto (cf "constraints.h"
      => "pluto_constraints_lexmin(const PlutoConstraints *cst, int negvar) ).
   
   The implementation was strongly inspired by some methods from Pluto's "contraints.c"
   source file */
#include "linAlg.h"
#include "include/piplib/piplib64.h"

#ifndef __LEXMIN_MPP_H
#define __LEXMIN_MPP_H


// Lexmax functions (using Piplib behind)
// Note that in our context, we just need an upper-bound/lower-bound of the resulting QUAST.
//		Thus, we simplify the QUAST given by Piplib into a single (rational) affine function
rational64** lexmax(polyhedronMPP *poly);
rational64** lexmax(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context);

rational64** lexmin(polyhedronMPP *poly);
rational64** lexmin(polyhedronMPP *poly, int64** context, int nrow_context, int ncol_context);


// Adaptation of the lexmax/min functions to max/min functions
rational64** argmax(affFuncMPP *obj, polyhedronMPP *dom);
rational64** max(affFuncMPP *obj, polyhedronMPP *dom);

rational64** argmin(affFuncMPP *obj, polyhedronMPP *dom);
rational64** min(affFuncMPP *obj, polyhedronMPP *dom);

#endif

