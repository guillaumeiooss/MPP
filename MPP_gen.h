// Standalone C++ version of Monoparametric tiling
// Author: Guillaume Iooss

#ifndef __MPP_GEN_H
#define __MPP_GEN_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "linAlg.h"
#include "lexmin.h"
#include "MPP_rect.h"

using namespace std;


// Monoparametric tiling applied to a polyhedron
//	- The input domain is "polyScalar"
//	- "shape" defines the (rectangular) shape of the tile
//	- "lattice" defines the lattice of tiles origins:
//		Each column of this lattice corresponds to one of the vector of the basis
//		Because some vectors might not be integral (ex: diamond tiling), an additional row is added at the bottom
//			of the matrix, representing the denominator
// The output of the transformation is an intersection (first "list") of union (second "list") of polyhedra. The simplification into a single union of polyhedra is doable, but would introduce many empty polyhedron and would require more polyhedral machinery.
list<list<polyhedronMPP*> > getTiledDomain(polyhedronMPP *polyScalar, polyhedronMPP *shape, int64** lattice, optionMPP *option);


#endif

