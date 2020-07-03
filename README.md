# MonoParametric Partitioning (MPP) library

# Overall description

This library implements the basic mathematical operations behind the MonoParametric Partitioning (MPP) transformation.
1. **MPP_rect.h/.cpp** corresponds to the (simpler) case where the tile shapes are rectangle.
2. **MPP_gen.h/.cpp** corresponds to the general case where the tile shapes can be anything polyhedral which tesselates the space.

In both case, we provides an operation for polyhedron, and for affine functions.

### Requirement

Piplib (http://www.piplib.org/) needs to be installed (version 1.4.0 works).



# Encoding of the partitioned polyhedron and affine function

### Polylib convention

We are using the Polylib convention to encode the coefficients of our polyhedron and affine functions.

For a polyhedron:
 1. The first column encodes if the constraints is an equality (`0`) of an inequality (`1`).
 2. The next columns correspond to the indices of the polyhedron (partitioned: blocked indices, then local indices)
 3. The next columns correspond to the parameters of the polyhedron (partitioned: blocked parameters, then local parameters, then block size parameter `b`)
 4. The last column corresponds to the constant coefficient.

For an affine function, the encoding is identical, but without the first column.


### (Custom) Extension to manage modulo constraints and rational expressions

Because the partitioned polyhedron we obtain might have additional modulo constraints, we have adapted the Polylib convention to encode them. We store the modulo coefficient in the first column (replacing the usual `0`). Because this coefficient must be strictly greater than 1, we can easily differentiate it from a common equality or inequality.

For example, the constraint `(i - N - 3) % 2 = 0` is encoded as the row `(2  1  -1  -3)`.

Also, in the case of general tiling for affine functions, the value of the partitioned affine function might contains rational expressions. In that case, we compute the smallest common denominator of the whole row and add it on its right.

For example, the value `( 2.i +N ) /2` inside an affine function is encoded as the row `(2  1  0) / 2`.



# Rectangular monoparametric tiling

### Polyhedron

The function ``getRectangularTiledDomain`` takes as an input:
 1. A polyhedron to be partitioned (typically, an iteration domain)
 2. The shape of the rectangular tiles, defined through their sizes. The size of this array must correspond to the number of dimensions of the polyhedron. For example, `[1, 2]` corresponds to a `b*2b` 2D rectangular tile.
 3. Some options to simplify the result.

It returns a list of list of polyhedron: the outer list corresponds to an intersection, the inner one corresponds to an union.

### Affine function

The function ``getRectangularTiledFunction`` takes as an input:
 1. An affine function to be partitioned (typically, a dependence function)
 2. The shape of the rectangular tiling of the input space of the affine function. The size of this array must correspond to the number of input dimension of the affine function.
 3. The shape of the rectangular tiling of the output space of the affine function. The size of this array must correspond to the number of output dimension of the affine function.
 4. Some options to simplify the result.

It returns a piecewise affine function, defined by using a map: the condition (polyhedron) is associated with the value of the function (affine function).


### Extension to unimodular parallelogram tiling

We can use the rectangular tiling for parallelogram tiling, if the hyperplanes defining the boundaries of this parallelogram form a unimodular matrix.

The function ``getRectangularCoBTiledDomain`` takes an additional argument (compared to ``getRectangularTiledDomain``) which are the matrix whose columns are the normal vectors of the hyperplanes. This matrix must be unimodular.

### Options
 1. `kMinMaxOption` : selects how conservative we are in the computation of kmin and kmax (the iterator over all the different tile shapes). If it is 0, we assume that the block size parameter `b` can be any positive value. If it is 1, we assume that `b` is big enough and we can skip the corner cases.
 2. `areParamDiv` : if it is set to true, we assume that the block size parameter `b` divides exactly the parameters of the polyhedron/affine function. This option is useful to avoid generating boundary tiles.
 3. `minBlSizeParam` : set a minimal value for the block size parameter `b`.
 4. `errorIfModulo` : triggers an exception if the resulting piecewise affine function contains modulo conditions. This is only used by the rectangular case and for affine functions.




# General monoparametric tiling

### Polyhedron

The function ``getTiledDomain`` takes as an input:
 1. A polyhedron to be partitioned
 2. The tile shape used for tiling. This shape must be monoparametric. This means that it must correspond to a non-parametric shape, scaled by a factor `b` (the tile size parameter). To obtain this shape, you can build it from its constraints (remember: it must only have one parameter). Alternatively, the rectangular and parallelogram shape have dedicated functions inside **linAlg.h/.cpp** to build them.
 3. The lattice of tile origins. This lattice defines how the tiles are positioned between each other. Each column of this matrix corresponds to one of the vector of the basis of the lattice. Because this vector might not be integral (ex: non-unimodular parallelogram tiling, cf diamond tiling), an additional row is added at the bottom of the matrix, representing the common denominator of the whole column. The condition is that the lattice combined with the tile shape must tessellate the space (i.e., cover the whole space with no overlapping [1]).
 4. Some options to simplify the result.

It (also) returns a list of list of polyhedron: the outer list corresponds to an intersection, the inner one corresponds to an union.

[1] Technically, our proof and implementation should support overlapped tiling. However, this is not fully tested yet.


### Affine functions

The function ``getTiledFunction`` takes as an input:
 1. An affine function to be partitioned
 2. The tile shape used for the tiling of the input space of the function
 3. The lattice of tile origins of the input space tiling
 4. The tile shape used for the tiling of the output space of the function
 5. The lattice of tile origins of the output space tiling
 6. Some options to simplify the result (note: `errorIfModulo` does not have any effect here).

It returns a piecewise affine function, defined by using a map: the condition (polyhedron) is associated with the value of the function (affine function).



# Other files

The file **linAlg.h/.cpp** contains the implementation of various basic linear algebraic operations (up to matrix inversion) and the implementation of the polyhedral and affine function data structures used internally (mostly a matrix of coefficients, similar to the Polylib format).

The file **lexmin.h/.cpp** manages the interface with the Piplib library, and exposes the lexmin/max and min/max functions using the structures defined in **linAlg.h**.

The file **test_MPP.cpp** contains some examples of use of the library (and their corresponding result).



# References

### Rectangular case:
@inproceedings{Iooss2014impact,
	author = {Iooss, Guillaume and Rajopadhye, Sanjay and Alias, Christophe and Zou, Yun},
	title = {Constant Aspect Ratio Tiling},
	booktitle = {Proceedings of the 
	    4th International Workshop on Polyhedral Compilation Techniques},
	editor = {Rajopadhye, Sanjay and Verdoolaege, Sven},
	year   = 2014,
	month  = Jan,
	address = {Vienna, Austria}
}

### Extension to any tile shape, but assuming integral tile origins:
@phdthesis{Iooss2016PhD,
  TITLE = {Detection of Linear Algebra Operations in Polyhedral Programs},
  AUTHOR = {Iooss, Guillaume},
  URL = {https://tel.archives-ouvertes.fr/tel-01370553},
  NUMBER = {2016LYSEN019},
  SCHOOL = {Université de Lyon},
  YEAR = {2016},
  MONTH = Jul,
  KEYWORDS = {Polyhedral model ; Tiling ; Program equivalence ; Template recognition ; BLAS ; Modèle polyédrique ; Tuilage ; équivalence de programme ; Reconnaissance de template},
  TYPE = {Theses},
  PDF = {https://tel.archives-ouvertes.fr/tel-01370553/file/IOOSS_Guillaume_2016LYSEN019_These.pdf},
  HAL_ID = {tel-01370553},
  HAL_VERSION = {v1},
}

### General proof:
@unpublished{iooss:hal-02493164,
  TITLE = {Monoparametric Tiling of Polyhedral Programs},
  AUTHOR = {Iooss, Guillaume and Alias, Christophe and Rajopadhye, Sanjay},
  URL = {https://hal.inria.fr/hal-02493164},
  NOTE = {working paper or preprint},
  YEAR = {2020},
  MONTH = Feb,
  KEYWORDS = {Tiling ; Program Transformation ; Polyhedral Model ; Compilation},
  PDF = {https://hal.inria.fr/hal-02493164/file/MPP_Hal_27_02_20.pdf},
  HAL_ID = {hal-02493164},
  HAL_VERSION = {v1},
}


# Contact

In case of bug/question, feel free to contact me at: guillaume [dot] iooss [at] gmail.com

Last update of the documentation: 3 July 2020.
