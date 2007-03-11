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
// A proxy version of the LU class,
// cleaned up to present a comprehensible
// version of the LU interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>
#include <lapack.h>
#include <TooN/TooN.h>

/// All classes and functions are within this namespace
namespace TooN
{
/**
@class LU LUdoc.h TooN/LU.h
Performs %LU decomposition and back substitutes to solve equations.
The %LU decomposition is the fastest way of solving the equation 
\f$M\underline{x} = \underline{c}\f$m, but it becomes unstable when
\f$M\f$ is (nearly) singular (in which cases the SymEigen or SVD decompositions
are better). It decomposes a matrix \f$M\f$ into
\f[M = L \times U\f]
where \f$L\f$ is a lower-diagonal matrix with unit diagonal and \f$U\f$ is an 
upper-diagonal matrix. The library only supports the decomposition of square matrices.
It can be used as follows to solve the \f$M\underline{x} = \underline{c}\f$ problem as follows:
@code
// construct M
double d1[][] = {{1,2,3},{4,5,6},{7,8,10}};
Matrix<3> M(d1);
// construct c
Vector<3> c = 2,3,4;
// create the LU decomposition of M
LU<3> luM(M);
// compute x = M^-1 * c
Vector<3> x = luM.backsub(c);
@endcode
The convention LU<> (=LU<-1>) is used to create an LU decomposition whose size is 
determined at runtime.
@ingroup gDecomps
**/

template <int Size>
class LU 
{
public:
	/// Construct the %LU decomposition of a matrix. This initialises the class, and
	/// performs the decomposition immediately.
	LU(const Matrix<Size,Size>& M);
  
	/// Perform the %LU decompsition of another matrix.
  	void compute(const Matrix<Size,Size>& M);

	/// Calculate result of multiplying the inverse of M by another matrix. For a matrix \f$A\f$, this
	/// calculates \f$M^{-1}A\f$ by back substitution (i.e. without explictly calculating the inverse).
  	template<int Cols>
  	Matrix<Size, Cols> backsub(const Matrix<Size, Cols>& rhs);
  
	/// Calculate result of multiplying the inverse of M by a vector. For a vector \f$b\f$, this
	/// calculates \f$M^{-1}b\f$ by back substitution (i.e. without explictly calculating the inverse).
	Vector<Size> backsub(const Vector<Size>& rhs);

	/// Calculate inverse of the matrix. This is not usually needed: if you need the inverse just to 
	/// multiply it by a matrix or a vector, use one of the backsub() functions, which will be faster.
	Matrix<Size,Size> get_inverse();
  
	/// Calculate the determinant of the matrix
	double determinant() const;


	/// Returns the L and U matrices. The permutation matrix is not returned.
	/// Since L is lower-triangular (with unit diagonal)
	/// and U is upper-triangular, these are returned conflated into one matrix, where the 
	/// diagonal and above parts of the matrix are U and the below-diagonal part, plus a unit diagonal, 
	/// are L.
  	Matrix<Size,Size>& get_lu();

};

}

#endif
