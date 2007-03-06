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
// A proxy version of the Cholesky class,
// cleaned up to present a comprehensible
// version of the Cholesky interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS


#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <limits>

/// All classes and functions are within this namespace
namespace TooN
{
/**
@class Cholesky Choleskydoc.h TooN/Cholesky.h
Decomposes a positive-semidefinite symmetric matrix A (such as a covariance) into L*D*L^T, where L is lower-triangular and D is diagonal.
Also can compute A = S*S^T, with S lower triangular.  The LDL^T form is faster to compute than the class Cholesky decomposition.
The decomposition can be used to compute A^-1*x, A^-1*M, M*A^-1*M^T, and A^-1 itself, though the latter rarely needs to be explicitly represented.
Also efficiently computes det(A) and rank(A).
It can be used as follows:
@code
// Declare some matrices.
Matrix<3> A = ...; // we'll pretend it is pos-def
Matrix<2,3> M;
Matrix<2> B;
Vector<3> y = (make_Vector, 2,3,4);
// create the Cholesky decomposition of A
Cholesky<3> chol(A);
// compute x = A^-1 * y
Vector<3> x = cholA.inverse_times(y);
// Identical to above
x = cholA.backsub(y);
// compute B = M*A^-1*M^T
B = cholA.transform_inverse(M);
//compute A^-1
Matrix<3> Ainv = cholA.get_inverse();
Matrix<3> C = ... // again, C is pos-def
//compute the 'square-root' of C
Matrix<3> L = Cholesky<3>::sqrt(C);
@endcode
@ingroup gDecomps
**/

template <int N>
class Cholesky {
public:
    /// Construct the Cholesky-ish decomposition of a matrix. This initialises the class, and
    /// performs the decomposition immediately.
    /// Run time is O(N^3)
    template<class Accessor> Cholesky(const FixedMatrix<N,N,Accessor>& m);
    
    /// Compute the LDL^T decomposition of another matrix.
    /// Run time is O(N^3)
    template<class Accessor> void compute(const FixedMatrix<N,N,Accessor>& m);

    /// Get the rank of the matrix
    int get_rank() const;

    /// Get the determinant
    double get_determinant() const;

    /// Get the internal representation of L and D.
    /// L has 1's on the diagonal, so here D is stored in its diagonal.
    /// L is lower triangular; the upper-triangular contents are not initialized.
    /// Not yet present in Cholesky<-1>
    const Matrix<N>& get_L_D() const { return L; }

    /// Compute the lower-triangular 'square-root' L of A, s.t. A = L*L^T
    /// Run time is O(N^3)
    /// Not yet present in Cholesky<-1>
    template <class A1, class A2> static void sqrt(const FixedMatrix<N,N,A1>& A, FixedMatrix<N,N,A2>& L);

    /// Compute the lower-triangular 'square-root' L, of A, s.t. A = L*L^T
    /// Run time is O(N^3)
    /// Not yet present in Cholesky<-1>
    template <class A1> static Matrix<N> sqrt(const FixedMatrix<N,N,A1>& A);

    /// Compute the lower-triangular 'square-root' L of the decomposed matrix, s.t. A = L*L^T
    /// L is stored in M
    /// Not yet present in Cholesky<-1>
    template <class A> void get_sqrt(FixedMatrix<N,N,A>& M) const;

    /// Compute the lower-triangular 'square-root' L of the decomposed matrix, s.t. A = L*L^T
    /// Not yet present in Cholesky<-1>
    Matrix<N> get_sqrt() const;
	
    /// Compute the lower-triangular 'square-root' L of the decomposed matrix, s.t. A = L*L^T
    /// @deprecated Use get_sqrt, or the static sqrt function
    Matrix<N> get_L() const __attribute__ ((deprecated));
	
    /// Compute the lower-triangular 'square-root' L of the inverse of the decomposed matrix, s.t. A^-1 = L*L^T
    /// Not yet present in Cholesky<-1>
    template <class A> void get_inv_sqrt(FixedMatrix<N,N,A>& M) const;
	
    /// Compute the Mahalanobis squared magnitude of x, d = x*A^-1*x
    /// Run time is O(N^2)
    /// Not yet present in Cholesky<-1>
    template <class A> double mahalanobis(const FixedVector<N,A>& v) const;

    /// Compute J*A^-1*J^T, and store in T using F::eval (for instance, F might be util::PlusAssign to add the result to T)
    /// Run time is O(MN^2 + NM^2)
    /// Not yet present in Cholesky<-1>
    template <class F, int M, class A1, class A2> void transform_inverse(const FixedMatrix<M,N,A1>& J, FixedMatrix<M,M,A2>& T) const;

    /// Compute J*A^-1*J^T, and store in T
    /// Run time is O(MN^2 + NM^2)
    /// Not yet present in Cholesky<-1>
    template <int M, class A1, class A2> void transform_inverse(const FixedMatrix<M,N,A1>& J, FixedMatrix<M,M,A2>& T) const;

    /// Compute J*A^-1*J^T
    /// Run time is O(MN^2 + NM^2)
    /// Not yet present in Cholesky<-1>
    template <int M, class A> Matrix<M> transform_inverse(const FixedMatrix<M,N,A>& J) const;

    /// Compute x = A^-1*v
    /// Run time is O(N^2)
    /// Not yet present in Cholesky<-1>
    template <class A1, class A2> inline void inverse_times(const FixedVector<N,A1>& v, FixedVector<N,A2>& x) const;

    /// Compute A^-1*v
    /// Run time is O(N^2)
    /// Not yet present in Cholesky<-1>
    template <class Accessor> inline Vector<N> inverse_times(const FixedVector<N,Accessor>& v) const;

    /// Compute x = A^-1*v
    /// Run time is O(N^2)
    template <class Accessor> inline Vector<N> backsub(const FixedVector<N,Accessor>& v) const;
      
    /// Compute A^-1*B
    /// Run time is O(MN^2)
    /// Not yet present in Cholesky<-1>
    template <class A, int M> inline Matrix<N,M> inverse_times(const FixedMatrix<N,M,A>& B);

    /// Compute A^-1 and store in M
    /// Run time is O(N^3)
    template <class A> void get_inverse(FixedMatrix<N,N,A>& M) const;

    /// Compute A^-1
    /// Run time is O(N^3)
    Matrix<N> get_inverse() const;
	
    /// Update the decomposition for A' = A + V*V^T
    /// Run time is O(MN^2), but for small N, it's faster to compute A' and recompute the decomposition
    /// Not yet present in Cholesky<-1>
    template <int M, class A>  void update(const FixedMatrix<N,M,A>& V);

    /// Update the decomposition for A' = A + v*v^T
    /// Run time is O(N^2), but for small N, it's faster to compute A' and recompute the decomposition
    /// Not yet present in Cholesky<-1>
    template <class A>  void update(const FixedVector<N,A>& v);
};

};

#endif
