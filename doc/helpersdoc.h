/*
    Copyright (c) 2005 Paul Smith

	Permission is granted to copy, distribute and/or modify this document under
	the terms of the GNU Free Documentation License, Version 1.2 or any later
	version published by the Free Software Foundation; with no Invariant
	Sections, no Front-Cover Texts, and no Back-Cover Texts.

    You should have received a copy of the GNU Free Documentation License
    License along with this library; if not, write to the Free Software
    Foundation, Inc.
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

*/
// A proxy version of the numhelpers class,
// cleaned up to present a comprehensible
// version of the interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>

/// All classes and functions are within this namespace
namespace TooN
{

/// @name Helper functions.
/// Defined in <code>TooN/helpers.h</code>
//@{

/**
Normalise a vector (make the sum of the squares of the elements one). The
result is a vector with the same direction, but whose length (2-norm) is one.
The normalised vector is
\f$ \hat{\underline{v}} = \dfrac{\underline{v}}{\|\underline{v}\|} =
  \dfrac{\underline{v}}{sqrt{\sum_i v_i}}\f$.
This function operates on the input vector.
@code
Vector<3> a = 0,3,4;
normalize(a);  // now a = [0.0 0.6 0.8]
@endcode
@relates Vector
*/
template<int Size>
void normalize(Vector<Size>& v);

/**
Divide all the elements in a vector by the last element. This function
operates on the input vector.
@code
Vector<3> a = 2,3,4;
normalize_last(a);  // now a = [0.5 0.75 1.0]
@endcode
@relates Vector
*/
template<int Size>
void normalize_last(Vector<Size>& v);

/**
Normalise the vector ignoring the last element. The modified vector will be
in the same direction, and the length of the sub-vector made up of all the
vectors but the last one will be one. This function operates on the input
vector.
@code
Vector<3> a = 3,4,2;
normalize_last(a);  // now a = [0.6 0.8 0.4]
@endcode
@relates Vector
This is useful if the vector in question represent the equation of a plane
in homogeneous co-ordinates. The equation of a plane \f$\underline{r} \cdot
\underline{n} = d\f$ can also be written in homogeneous co-ordinates as
\f[\begin{bmatrix}\underline{r} \\
1\end{bmatrix}\cdot\begin{bmatrix}\underline{n} \\ -d\end{bmatrix}  = 0\f]
The dot product of the plane vector
\f$\begin{bmatrix}\underline{n} & -d\end{bmatrix}^T\f$  with a homogeneous point
then gives the distance from the point to the plane, if the normal vector
\f$\underline{n}\f$ is of unit length. This function allows this vector to
be normalised.
*/
template<int Size>
void normalize_but_last(Vector<Size>& v);

/**
Project a homogeneous vector down to a non-homogeneous one. This divide all
the elements in a vector by the last element and returns all the elements apart
from the last one (in a homogeneous vector, the last element represents the
scale factor; a non-homogeneous vector is assumed to have a scale factor of
one).
@code
Vector<3> a = 2,3,4;
Vector<2> b = project(a);  // now b = [0.5 0.75]
@endcode
@relates Vector
*/
template<int Size>
Vector<Size-1> project(Vector<Size>& v);

/**
Convert a non-homogeneous vector into a homogeneous vector. This returns the
same vector, augmented by a one more element with value 1.0.
@code
Vector<3> a = 2,3,4;
Vector<4> b = unproject(a);  // now b = [2, 3, 4, 1]
@endcode
@relates Vector
*/
template<int Size>
Vector<Size+1> unproject(Vector<Size>& v);

/**
Treat an array of doubles as a vector. This avoids having to construct a
Vector if it is not needed. The array should be the same length as Size.
@relates Vector
*/
template<int Size> 
Vector<Size>&  as_vector(double* data);

/**
Treat an array of doubles as a vector. This avoids having to construct a
Vector if it is not needed. The array should be the same length as Size.
@relates Vector
*/
template<int Size> 
Vector<Size>&  as_vector(const double* data);

/**
Set a matrix to the Identity (or some multiple of). This replaces
<code>m</code> with a matrix with <code>factor</code> (= 1 by default) down the
diagonal and zeros elsewhere. This function is only defined for square,
fixed-size matrices.
*/
@relates Matrix
template <int Size> 
void Identity(Matrix<Size,Size>&m, const double factor=1);

/**
Make a matrix symmetrical. This leaves the diagonal elements unchanged and
averages the other elements across the diagonal
i.e. \f$m_{ij}^\text{new} = m_{ji}^\text{new} = (m_{ij}^\text{old} +
m_{ji}^\text{old}) / 2 \f$.
This function is only defined for square, fixed-size matrices.
@relates Matrix
*/
template <int Size> 
void Symmetrize(Matrix<Size,Size>& m);

/**
Transpose a matrix. The preferred means of transposing a matrix is to use
Matrix::T() which is very efficient (it just defines a different memory layout).
This means that it is not possible to say <code>M = M.T();</code> (and it is not
usually necessary to do so). If this operation is required, this function can be
used. This function is only defined for square, fixed-size matrices.
@relates Matrix
*/
template <int Size> 
void Transpose(Matrix<Size,Size>& m);

/**
Set all of the elements of a vector to zero
@code
Vector<3> a = 1,2,3;
Zero(a);   // now a = [0 0 0];
@endcode
@relates Vector
*/
template <int Size> 
void Zero(Vector<Size>& v);


//@}


}

#endif
