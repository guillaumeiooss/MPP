// Standalone C++ version of Monoparametric tiling
// Author: Guillaume Iooss

#ifndef __MPP_RECT_H
#define __MPP_RECT_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "linAlg.h"

using namespace std;


// Structure for the options of the monoparametric tiling transformation, plus default values
//		kMinMaxOption = 0 or 1 ==> Select how conservative we are in our hypothesis when computing kmin and kmax
//		areParamDiv            ==> Are the program parameter multiple of the block size parameter b
//		minBlSizeParam         ==> Minimal value of the block size parameter b
//
// 	kMinMaxOption=0 => Compute a kMin/kMax for every possible values of "b".
// 	kMinMaxOption=1 => Compute a tighter kMin/kMax by assuming that "b" is "big enough".
struct optionMPP {
  int kMinMaxOption = 0;
  bool areParamDiv = false;
  int minBlSizeParam = 4;
  bool errorIfModulo = false;
};


// Monoparametric tiling applied to a polyhedron
//	- The input domain is "polyScalar"
//	- "scale" defines the (rectangular) shape of the tile
//	- "option" are the options of the monoparametric tiling transformation
// The output of the transformation is an intersection (first "list") of union (second "list") of polyhedra. The simplification into a single union of polyhedra is doable, but would introduce many empty polyhedron and would require more polyhedral machinery.
list<list<polyhedronMPP*> > getRectangularTiledDomain(polyhedronMPP *polyScalar, int64 *scale, optionMPP *option);



// Monoparametric tiling applied to an affine function
//	- The input affine function is "affScalar"
//	- "scale" defines the (rectangular) shape of the tile for the input space of the function
//	- "scaleIm" defines the (rectangular) shape of the tile for the output space of the function
//	- "option" are the options of the monoparametric tiling transformation
// The output of the transformation is a piecewise affine function: the branches are stored in a "map",
//		the first element of the pair is the condition of the branch, and
//		the second element of the pair is the value of the function on the branch.
map<polyhedronMPP*, affFuncMPP*> getRectangularTiledFunction(affFuncMPP *affScalar, int64 *scale, int64 *scaleIm, optionMPP *option);



// Monoparametric tiling applied to a polyhedron
//	- "hyperplanes" (the additional argument compared to getRectangularTiledDomain) is defining the tile shape
//		using a set of unimodular hyperplanes vectors.
list<list<polyhedronMPP*> > getRectangularCoBTiledFunction(polyhedronMPP* polyScalar, int64 **hyperplanes,
		int64 *scale, optionMPP* option);

#endif

