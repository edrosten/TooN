//-*- c++ -*-

//     Copyright (C) 2002 Tom Drummond (twd20@eng.cam.ac.uk)

//     This library is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation; either
//     version 2.1 of the License, or (at your option) any later version.

//     This library is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.

//     You should have received a copy of the GNU Lesser General Public
//     License along with this library; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef __LAPACK_H
#define __LAPACK_H

// LAPACK and BLAS routines

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

  // inverse of a triangular matrix * a vector (BLAS level 2)
  void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* lda, double* B, int* ldb);

  // SVD of a general matrix of doubles
  void dgesvd_(char* JOBU, char* JOBVT, int* M, int *N, double* A, int* lda,
	       double* S, double *U, int* ldu, double* VT, int* ldvt,
	       double* WORK, int* lwork, int* INFO);

  // Eigen decomposition of a symmetric matrix of doubles
  void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* lda, double* W,
	      double* WORK, int* LWORK, int* INFO);

}






#endif
