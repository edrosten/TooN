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
// A proxy version of the SVD class,
// cleaned up to present a comprehensible
// version of the SVD interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>
#include <lapack.h>
#include <TooN/toon.h>

/// All classes and functions are within this namespace
namespace TooN 
{
/**
@class SVD SVDdoc.h TooN/SVD.h
Performs %SVD and back substitute to solve equations.
Singular value decompositions are more robust than LU decompositions in the face of 
singular or nearly singular matrices. They decompose a matrix (of any shape) \f$M\f$ into:
\f[M = U \times D \times V^T\f]
where \f$D\f$ is a diagonal matrix of positive numbers whose dimension is the minimum 
of the dimensions of \f$M\f$. If \f$M\f$ is tall and thin (more rows than columns) 
then \f$U\f$ has the same shape as \f$M\f$ and \f$V\f$ is square (vice-versa if \f$M\f$ 
is short and fat). The columns of \f$U\f$ and the rows of \f$V\f$ are orthogonal 
and of unit norm (so one of them lies in SO(N)). The inverse of \f$M\f$ (or pseudo-inverse 
if \f$M\f$ is not square) is then given by
\f[M^{\dagger} = V \times D^{-1} \times U^T\f]
 
If \f$M\f$ is nearly singular then the diagonal matrix \f$D\f$ has some small values 
(relative to its largest value) and these terms dominate \f$D^{-1}\f$. To deal with 
this problem, the inverse is conditioned by setting a maximum ratio 
between the largest and smallest values in \f$D\f$ (passed as the <code>condition</code>
parameter to the various functions). Any values which are too small 
are set to zero in the inverse (rather than a large number)
 
It can be used as follows to solve the \f$M\underline{x} = \underline{c}\f$ problem as follows:
@code
// construct M
double d1[][] = {{1,2,3},{4,5,6},{7,8,10}};
Matrix<3> M(d1);
// construct c
 Vector<3> c;
c = 2,3,4;
// create the SVD decomposition of M
SVD<3> svdM(M);
// compute x = M^-1 * c
Vector<3> x = svdM.backsub(c);
 @endcode

SVD<> (= SVD<-1>) can be used to create an SVD whose size is determined at run-time.
@ingroup gDecomps
**/
template <int Rows, int Cols>
class SVD
{
public:
	/// Default constructor. Does nothing.
	SVD(){}
  
	/// Construct the %SVD decomposition of a matrix. This initialises the class, and
	/// performs the decomposition immediately.
  	SVD(const Matrix<Rows,Cols>& M);

	/// Compute the %SVD decomposition of M, typically used after the default constructor
	void compute(const Matrix<Rows,Cols>& M);

	/// Calculate result of multiplying the (pseudo-)inverse of M by another matrix. 
	/// For a matrix \f$A\f$, this calculates \f$M^{\dagger}A\f$ by back substitution 
	/// (i.e. without explictly calculating the (pseudo-)inverse). 
	/// See the detailed description for a description of condition variables.
    Matrix<Cols,RHS> backsub(const Matrix<Rows,RHS>& rhs, const double condition = 1e9);
    
	/// Calculate result of multiplying the (pseudo-)inverse of M by a vector. 
	/// For a vector \f$b\f$, this calculates \f$M^{\dagger}b\f$ by back substitution 
	/// (i.e. without explictly calculating the (pseudo-)inverse). 
	/// See the detailed description for a description of condition variables.
    Vector<Cols> backsub(const Vector<Rows>& v, const double condition = 1e9);

	/// Calculate (pseudo-)inverse of the matrix. This is not usually needed: 
	/// if you need the inverse just to multiply it by a matrix or a vector, use 
	/// one of the backsub() functions, which will be faster.
	/// See the detailed description of the pseudo-inverse and condition variables.
	Matrix<Cols,Rows> get_pinv(const double condition = 1e9);
};

}

#endif
