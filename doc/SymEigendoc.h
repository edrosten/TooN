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
// A proxy version of the SymEigen class,
// cleaned up to present a comprehensible
// version of the SymEigen interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>
#include <lapack.h>
#include <TooN/toon.h>

/// All classes and functions are within this namespace
namespace TooN 
{

/**
@class SymEigen SymEigendoc.h TooN/SymEigen.h
Performs eigen decomposition of a matrix.
Real symmetric (and hence square matrices) can be decomposed into
\f[M = U \times \Lambda \times U^T\f]
where \f$U\f$ is an orthogonal matrix (and hence \f$U^T = U^{-1}\f$) whose columns
are the eigenvectors of \f$M\f$ and \f$\Lambda\f$ is a diagonal matrix whose entries
are the eigenvalues of \f$M\f$. These quantities are often of use directly, and can
be obtained as follows:
@code
// construct M
double d1[][] = {{1,2,3},{2,5,6},{3,6,7}};
Matrix<3> M(d1);
// create the eigen decomposition of M
SymEigen<3> eigM(M);
// print the smallest eigenvalue
cout << eigM.get_evalues()[0] << endl;
// print the associated eigenvector
cout << eigM.get_evectors()[0] << endl;
@endcode

This decomposition is very similar to the SVD (q.v.), and can be used to solve
equations using backsub() or get_pinv(), with the same treatment of condition numbers.

SymEigen<> (= SymEigen<-1>) can be used to create an eigen decomposition whose size is determined at run-time.
@ingroup gDecomps
**/
template <int Size>
class SymEigen
{
public:
	/// Default constructor. Does nothing.
	SymEigen(){}
  
	/// Construct the eigen decomposition of a matrix. This initialises the class, and
	/// performs the decomposition immediately.
  	SymEigen(const Matrix<Size,Size>& M);

	/// Perform the eigen decomposition of a matrix.
  	void compute(const Matrix<Size,Size>& M);

	/// Calculate result of multiplying the (pseudo-)inverse of M by another matrix. 
	/// For a matrix \f$A\f$, this calculates \f$M^{\dagger}A\f$ by back substitution 
	/// (i.e. without explictly calculating the (pseudo-)inverse). 
	/// See the SVD detailed description for a description of condition variables.
    Matrix<Size,RHS> backsub(const Matrix<Size,RHS>& rhs, const double condition = 1e9);
    
	/// Calculate result of multiplying the (pseudo-)inverse of M by a vector. 
	/// For a vector \f$b\f$, this calculates \f$M^{\dagger}b\f$ by back substitution 
	/// (i.e. without explictly calculating the (pseudo-)inverse). 
	/// See the SVD detailed description for a description of condition variables.
    Vector<Size> backsub(const Vector<Size>& v, const double condition = 1e9);

	/// Calculate (pseudo-)inverse of the matrix. This is not usually needed: 
	/// if you need the inverse just to multiply it by a matrix or a vector, use 
	/// one of the backsub() functions, which will be faster.
	/// See the SVD detailed description for a description of the pseudo-inverse 
	/// and condition variables.
	Matrix<Size,Size> get_pinv(const double condition = 1e9);

	/// Calculates the reciprocals of the eigenvalues of the matrix.
	/// The vector <code>invdiag</code> lists the eigenvalues in order, from
	/// the largest (i.e. smallest reciprocal) to the smallest.
	/// These are also the diagonal values of the matrix \f$Lambda^{-1}\f$. 
	/// Any eigenvalues which are too small are set to zero (see the SVD 
	/// detailed description for a description of the and condition variables).
	void get_inv_diag(Vector<Size>& invdiag, double condition=1e9);
  
	/// Returns the eigenvectors of the matrix.
	/// This returns \f$U^T\f$, so that the rows of the matrix are the eigenvectors,
	/// which can be extracted using usual Matrix::operator[]() subscript operator.
	/// They are returned in order of the size of the corresponding eigenvalue, i.e.
	/// the vector with the largest eigenvalue is first.
	Matrix<Size,Size,RowMajor>& get_evectors() {return my_evectors;}

	/// Returns the eigenvectors of the matrix.
	/// This returns \f$U^T\f$, so that the rows of the matrix are the eigenvectors,
	/// which can be extracted using usual Matrix::operator[]() subscript operator.
	/// They are returned in order of the size of the corresponding eigenvalue, i.e.
	/// the vector with the largest eigenvalue is first.
	const Matrix<Size,Size,RowMajor>& get_evectors();

	/// Returns the eigenvalues of the matrix.
	/// The eigenvalues are listed in order, from the largest to the smallest.
	/// These are also the diagonal values of the matrix \f$\Lambda\f$. 
	Vector<Size>& get_evalues();

	/// Returns the eigenvalues of the matrix.
	/// The eigenvalues are listed in order, from the largest to the smallest.
	/// These are also the diagonal values of the matrix \f$\Lambda\f$. 
	const Vector<Size>& get_evalues() const;

};

}

#endif
