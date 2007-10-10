
/*                       
	 Copyright (C) 2005 Tom Drummond

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __LAPACK_H
#define __LAPACK_H

// LAPACK and BLAS routines

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  // inverse of a triangular matrix * a vector (BLAS level 2)
  void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* lda, double* B, int* ldb);

  // SVD of a general matrix of doubles
  void dgesvd_(const char* JOBU, const char* JOBVT, int* M, int *N, double* A, int* lda,
	       double* S, double *U, int* ldu, double* VT, int* ldvt,
	       double* WORK, int* lwork, int* INFO);

  // Eigen decomposition of a symmetric matrix of doubles
  void dsyev_(const char* JOBZ, const char* UPLO, int* N, double* A, int* lda, double* W,
	      double* WORK, int* LWORK, int* INFO);

    // Cholesky decomposition
    void dpotrf_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO);

    // Cholesky solve AX=B given decomposition
    void dpotrs_(const char* UPLO, const int* N, const int* NRHS, const double* A, const int* LDA, double* B, const int* LDB, int* INFO);

    // Cholesky inverse given decomposition
    void dpotri_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO);
}


#ifndef TOON_NO_NAMESPACE
}
#endif 



#endif
