// -*- c++ -*-

// Copyright (C) 2005,2009,2010 Tom Drummond (twd20@cam.ac.uk), E. Rosten

//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND OTHER CONTRIBUTORS ``AS IS''
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR OTHER CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.

#ifndef TOON_INCLUDE_LAPCK_H
#define TOON_INCLUDE_LAPCK_H

#include <TooN/TooN.h>

// LAPACK and BLAS routines
namespace TooN {

	extern "C" {
		// LU decomoposition of a general matrix
		void dgetrf_(FortranInteger* M, FortranInteger *N, double* A, FortranInteger* lda, FortranInteger* IPIV, FortranInteger* INFO);
		void sgetrf_(FortranInteger* M, FortranInteger *N, float* A, FortranInteger* lda, FortranInteger* IPIV, FortranInteger* INFO);

		// generate inverse of a matrix given its LU decomposition
		void dgetri_(FortranInteger* N, double* A, FortranInteger* lda, FortranInteger* IPIV, double* WORK, FortranInteger* lwork, FortranInteger* INFO);
		void sgetri_(FortranInteger* N, float* A, FortranInteger* lda, FortranInteger* IPIV, float* WORK, FortranInteger* lwork, FortranInteger* INFO);

		// inverse of a triangular matrix * a vector (BLAS level 2)
		void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, FortranInteger* M, FortranInteger* N, double* alpha, double* A, FortranInteger* lda, double* B, FortranInteger* ldb);
		void strsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, FortranInteger* M, FortranInteger* N, float* alpha, float* A, FortranInteger* lda, float* B, FortranInteger* ldb);
  

		// SVD of a general matrix
		void dgesvd_(const char* JOBU, const char* JOBVT, FortranInteger* M, FortranInteger *N, double* A, FortranInteger* lda,
					 double* S, double *U, FortranInteger* ldu, double* VT, FortranInteger* ldvt,
					 double* WORK, FortranInteger* lwork, FortranInteger* INFO);

		void sgesvd_(const char* JOBU, const char* JOBVT, FortranInteger* M, FortranInteger *N, float* A, FortranInteger* lda,
					 float* S, float *U, FortranInteger* ldu, float* VT, FortranInteger* ldvt,
					 float* WORK, FortranInteger* lwork, FortranInteger* INFO);

		// Eigen decomposition of a symmetric matrix
		void dsyev_(const char* JOBZ, const char* UPLO, FortranInteger* N, double* A, FortranInteger* lda, double* W, double* WORK, FortranInteger* LWORK, FortranInteger* INFO);
		void ssyev_(const char* JOBZ, const char* UPLO, FortranInteger* N, float* A, FortranInteger* lda, float* W, float* WORK, FortranInteger* LWORK, FortranInteger* INFO);

		// Eigen decomposition of a non-symmetric matrix
		void dgeev_(const char* JOBVL, const char* JOBVR, FortranInteger* N, double* A, FortranInteger* lda, double* WR, double* WI, double* VL, FortranInteger* LDVL, double* VR, FortranInteger* LDVR , double* WORK, FortranInteger* LWORK, FortranInteger* INFO);
		void sgeev_(const char* JOBVL, const char* JOBVR, FortranInteger* N, float* A, FortranInteger* lda, float* WR, float* WI, float* VL, FortranInteger* LDVL, float* VR, FortranInteger* LDVR , float* WORK, FortranInteger* LWORK, FortranInteger* INFO);

		// Cholesky decomposition
		void dpotrf_(const char* UPLO, const FortranInteger* N, double* A, const FortranInteger* LDA, FortranInteger* INFO);
		void spotrf_(const char* UPLO, const FortranInteger* N, float* A, const FortranInteger* LDA, FortranInteger* INFO);

		// Cholesky solve AX=B given decomposition
		void dpotrs_(const char* UPLO, const FortranInteger* N, const FortranInteger* NRHS, const double* A, const FortranInteger* LDA, double* B, const FortranInteger* LDB, FortranInteger* INFO);
		void spotrs_(const char* UPLO, const FortranInteger* N, const FortranInteger* NRHS, const float* A, const FortranInteger* LDA, float* B, const FortranInteger* LDB, FortranInteger* INFO);

		// Cholesky inverse given decomposition
		void dpotri_(const char* UPLO, const FortranInteger* N, double* A, const FortranInteger* LDA, FortranInteger* INFO);
		void spotri_(const char* UPLO, const FortranInteger* N, float* A, const FortranInteger* LDA, FortranInteger* INFO);
		
		// Computes a QR decomposition of a general rectangular matrix with column pivoting
		void sgeqp3_(FortranInteger* M, FortranInteger* N, float* A, FortranInteger* LDA, FortranInteger* JPVT, float* TAU, float* WORK, FortranInteger* LWORK, FortranInteger* INFO );
		void dgeqp3_(FortranInteger* M, FortranInteger* N, double* A, FortranInteger* LDA, FortranInteger* JPVT, double* TAU, double* WORK, FortranInteger* LWORK, FortranInteger* INFO );
		
		//Reconstruct Q from a QR decomposition
		void sorgqr_(FortranInteger* M,FortranInteger* N,FortranInteger* K, float* A, FortranInteger* LDA, float* TAU, float* WORK, FortranInteger* LWORK, FortranInteger* INFO );
		void dorgqr_(FortranInteger* M,FortranInteger* N,FortranInteger* K, double* A, FortranInteger* LDA, double* TAU, double* WORK, FortranInteger* LWORK, FortranInteger* INFO );
	}


	//////////////////////////////////////////////////////////////////////////////////
	// C++ overloaded functions to access single and double precision automatically //
	//////////////////////////////////////////////////////////////////////////////////

	inline void getrf_(FortranInteger* M, FortranInteger *N, float* A, FortranInteger* lda, FortranInteger* IPIV, FortranInteger* INFO){
		sgetrf_(M, N, A, lda, IPIV, INFO);
	}

	inline void getrf_(FortranInteger* M, FortranInteger *N, double* A, FortranInteger* lda, FortranInteger* IPIV, FortranInteger* INFO){
		dgetrf_(M, N, A, lda, IPIV, INFO);
	}

	inline void trsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG, FortranInteger* M, FortranInteger* N, float* alpha, float* A, FortranInteger* lda, float* B, FortranInteger* ldb) { 
		strsm_(const_cast<char*>(SIDE), const_cast<char*>(UPLO), const_cast<char*>(TRANSA), const_cast<char*>(DIAG), M, N, alpha, A, lda, B, ldb);
	}

	inline void trsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG, FortranInteger* M, FortranInteger* N, double* alpha, double* A, FortranInteger* lda, double* B, FortranInteger* ldb) {
		dtrsm_(const_cast<char*>(SIDE), const_cast<char*>(UPLO), const_cast<char*>(TRANSA), const_cast<char*>(DIAG), M, N, alpha, A, lda, B, ldb);
	}

	inline void getri_(FortranInteger* N, double* A, FortranInteger* lda, FortranInteger* IPIV, double* WORK, FortranInteger* lwork, FortranInteger* INFO){
		dgetri_(N, A, lda, IPIV, WORK, lwork, INFO);
	}

	inline void getri_(FortranInteger* N, float* A, FortranInteger* lda, FortranInteger* IPIV, float* WORK, FortranInteger* lwork, FortranInteger* INFO){
		sgetri_(N, A, lda, IPIV, WORK, lwork, INFO);
	}

	inline void potrf_(const char * UPLO, const FortranInteger* N, double* A, const FortranInteger* LDA, FortranInteger* INFO){
		dpotrf_(UPLO, N, A, LDA, INFO);
	}

	inline void potrf_(const char * UPLO, const FortranInteger* N, float* A, const FortranInteger* LDA, FortranInteger* INFO){
		spotrf_(UPLO, N, A, LDA, INFO);
	}

	// SVD
	inline void gesvd_(const char* JOBU, const char* JOBVT, FortranInteger* M, FortranInteger *N, double* A, FortranInteger* lda,
				double* S, double *U, FortranInteger* ldu, double* VT, FortranInteger* ldvt,
				double* WORK, FortranInteger* lwork, FortranInteger* INFO){
		dgesvd_(JOBU, JOBVT, M, N, A, lda, S, U, ldu, VT, ldvt, WORK, lwork, INFO);
	}

	inline void gesvd_(const char* JOBU, const char* JOBVT, FortranInteger* M, FortranInteger *N, float* A, FortranInteger* lda,
					 float* S, float *U, FortranInteger* ldu, float* VT, FortranInteger* ldvt,
					 float* WORK, FortranInteger* lwork, FortranInteger* INFO){
		sgesvd_(JOBU, JOBVT, M, N, A, lda, S, U, ldu, VT, ldvt, WORK, lwork, INFO);
	}

	// Cholesky solve AX=B given decomposition
	inline void potrs_(const char* UPLO, const FortranInteger* N, const FortranInteger* NRHS, const double* A, const FortranInteger* LDA, double* B, const FortranInteger* LDB, FortranInteger* INFO){
		dpotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
	}

	inline void potrs_(const char* UPLO, const FortranInteger* N, const FortranInteger* NRHS, const float* A, const FortranInteger* LDA, float* B, const FortranInteger* LDB, FortranInteger* INFO){
		spotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
	}

	// Cholesky inverse given decomposition
	inline void potri_(const char* UPLO, const FortranInteger* N, double* A, const FortranInteger* LDA, FortranInteger* INFO){
		dpotri_(UPLO, N, A, LDA, INFO);
	}

	inline void potri_(const char* UPLO, const FortranInteger* N, float* A, const FortranInteger* LDA, FortranInteger* INFO){
		spotri_(UPLO, N, A, LDA, INFO);
	}

	inline void syev_(const char* JOBZ, const char* UPLO, FortranInteger* N, double* A, FortranInteger* lda, double* W, double* WORK, FortranInteger* LWORK, FortranInteger* INFO){
		dsyev_(JOBZ, UPLO, N, A, lda, W, WORK, LWORK, INFO);
	}
	inline void syev_(const char* JOBZ, const char* UPLO, FortranInteger* N, float* A, FortranInteger* lda, float* W, float* WORK, FortranInteger* LWORK, FortranInteger* INFO){
		ssyev_(JOBZ, UPLO, N, A, lda, W, WORK, LWORK, INFO);
	}

	//QR decomposition
	inline void geqp3_(FortranInteger* M, FortranInteger* N, float* A, FortranInteger* LDA, FortranInteger* JPVT, float* TAU, float* WORK, FortranInteger* LWORK, FortranInteger* INFO )
	{
		sgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO);
	}

	inline void geqp3_(FortranInteger* M, FortranInteger* N, double* A, FortranInteger* LDA, FortranInteger* JPVT, double* TAU, double* WORK, FortranInteger* LWORK, FortranInteger* INFO )
	{
		dgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO);
	}
	
	inline void orgqr_(FortranInteger* M,FortranInteger* N,FortranInteger* K, float* A, FortranInteger* LDA, float* TAU, float* WORK, FortranInteger* LWORK, FortranInteger* INFO )
	{
		sorgqr_(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
	}

	inline void orgqr_(FortranInteger* M,FortranInteger* N,FortranInteger* K, double* A, FortranInteger* LDA, double* TAU, double* WORK, FortranInteger* LWORK, FortranInteger* INFO )
	{
		dorgqr_(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
	}

	//Non symmetric (general) eigen decomposition
	inline void geev_(const char* JOBVL, const char* JOBVR, FortranInteger* N, double* A, FortranInteger* lda, double* WR, double* WI, double* VL, FortranInteger* LDVL, double* VR, FortranInteger* LDVR , double* WORK, FortranInteger* LWORK, FortranInteger* INFO){
		dgeev_(JOBVL, JOBVR, N,  A,  lda,  WR,  WI,  VL,  LDVL,  VR,  LDVR ,  WORK,  LWORK,  INFO);
	}

	inline void geev_(const char* JOBVL, const char* JOBVR, FortranInteger* N, float* A,  FortranInteger* lda, float* WR,  float* WI,  float* VL,  FortranInteger* LDVL, float* VR,  FortranInteger* LDVR , float* WORK,  FortranInteger* LWORK, FortranInteger* INFO){
		sgeev_(JOBVL, JOBVR, N,  A,  lda,  WR,  WI,  VL,  LDVL,  VR,  LDVR ,  WORK,  LWORK,  INFO);
	}
}
#endif
