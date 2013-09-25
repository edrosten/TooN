// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

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

#ifndef TOON_INCLUDE_LU_H
#define TOON_INCLUDE_LU_H

#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>

namespace TooN {
/**
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
  Matrix<3> M;
  M[0] = makeVector(1,2,3);
  M[1] = makeVector(3,2,1);
  M[2] = makeVector(1,0,1);
  // construct c
  Vector<3> c = makeVector(2,3,4);
  // create the LU decomposition of M
  LU<3> luM(M);
  // compute x = M^-1 * c
  Vector<3> x = luM.backsub(c);
@endcode
The convention LU<> (=LU<-1>) is used to create an LU decomposition whose size is 
determined at runtime.
@ingroup gDecomps
**/
template <int Size=-1, class Precision=double>
class LU {
	public:

	/// Construct the %LU decomposition of a matrix. This initialises the class, and
	/// performs the decomposition immediately.
	template<int S1, int S2, class Base>
	LU(const Matrix<S1,S2,Precision, Base>& m)
	:my_lu(m.num_rows(),m.num_cols()),my_IPIV(m.num_rows()){
		compute(m);
	}
	
	/// Perform the %LU decompsition of another matrix.
	template<int S1, int S2, class Base>
	void compute(const Matrix<S1,S2,Precision,Base>& m){
		//check for consistency with Size
		SizeMismatch<Size, S1>::test(my_lu.num_rows(),m.num_rows());
		SizeMismatch<Size, S2>::test(my_lu.num_rows(),m.num_cols());
	
		//Make a local copy. This is guaranteed contiguous
		my_lu=m;
		FortranInteger lda = m.num_rows();
		FortranInteger M = m.num_rows();
		FortranInteger N = m.num_rows();

		getrf_(&M,&N,&my_lu[0][0],&lda,&my_IPIV[0],&my_info);

		if(my_info < 0){
			std::cerr << "error in LU, INFO was " << my_info << std::endl;
		}
	}

	/// Calculate result of multiplying the inverse of M by another matrix. For a matrix \f$A\f$, this
	/// calculates \f$M^{-1}A\f$ by back substitution (i.e. without explictly calculating the inverse).
	template <int Rows, int NRHS, class Base>
	Matrix<Size,NRHS,Precision> backsub(const Matrix<Rows,NRHS,Precision,Base>& rhs){
		//Check the number of rows is OK.
		SizeMismatch<Size, Rows>::test(my_lu.num_rows(), rhs.num_rows());
	
		Matrix<Size, NRHS, Precision> result(rhs);

		FortranInteger M=rhs.num_cols();
		FortranInteger N=my_lu.num_rows();
		double alpha=1;
		FortranInteger lda=my_lu.num_rows();
		FortranInteger ldb=rhs.num_cols();
		trsm_("R","U","N","N",&M,&N,&alpha,&my_lu[0][0],&lda,&result[0][0],&ldb);
		trsm_("R","L","N","U",&M,&N,&alpha,&my_lu[0][0],&lda,&result[0][0],&ldb);

		// now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
		for(int i=N-1; i>=0; i--){
			const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
			for(int j=0; j<NRHS; j++){
				Precision temp = result[i][j];
				result[i][j] = result[swaprow][j];
				result[swaprow][j] = temp;
			}
		}
		return result;
	}

	/// Calculate result of multiplying the inverse of M by a vector. For a vector \f$b\f$, this
	/// calculates \f$M^{-1}b\f$ by back substitution (i.e. without explictly calculating the inverse).
	template <int Rows, class Base>
	Vector<Size,Precision> backsub(const Vector<Rows,Precision,Base>& rhs){
		//Check the number of rows is OK.
		SizeMismatch<Size, Rows>::test(my_lu.num_rows(), rhs.size());
	
		Vector<Size, Precision> result(rhs);

		FortranInteger M=1;
		FortranInteger N=my_lu.num_rows();
		double alpha=1;
		FortranInteger lda=my_lu.num_rows();
		FortranInteger ldb=1;
		trsm_("R","U","N","N",&M,&N,&alpha,&my_lu[0][0],&lda,&result[0],&ldb);
		trsm_("R","L","N","U",&M,&N,&alpha,&my_lu[0][0],&lda,&result[0],&ldb);

		// now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
		for(int i=N-1; i>=0; i--){
			const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
			Precision temp = result[i];
			result[i] = result[swaprow];
			result[swaprow] = temp;
		}
		return result;
	}

	/// Calculate inverse of the matrix. This is not usually needed: if you need the inverse just to 
	/// multiply it by a matrix or a vector, use one of the backsub() functions, which will be faster.
	Matrix<Size,Size,Precision> get_inverse(){
		Matrix<Size,Size,Precision> Inverse(my_lu);
		FortranInteger N = my_lu.num_rows();
		FortranInteger lda=my_lu.num_rows();
		FortranInteger lwork=-1;
		Precision size;
		getri_(&N, &Inverse[0][0], &lda, &my_IPIV[0], &size, &lwork, &my_info);
		lwork=FortranInteger(size);
		Precision* WORK = new Precision[lwork];
		getri_(&N, &Inverse[0][0], &lda, &my_IPIV[0], WORK, &lwork, &my_info);
		delete [] WORK;
		return Inverse;
	}

	/// Returns the L and U matrices. The permutation matrix is not returned.
	/// Since L is lower-triangular (with unit diagonal)
	/// and U is upper-triangular, these are returned conflated into one matrix, where the 
	/// diagonal and above parts of the matrix are U and the below-diagonal part, plus a unit diagonal, 
	/// are L.
	const Matrix<Size,Size,Precision>& get_lu()const {return my_lu;}
	
	private:
	inline int get_sign() const {
		int result=1;
		for(int i=0; i<my_lu.num_rows()-1; i++){
			if(my_IPIV[i] > i+1){
				result=-result;
			}
		}
		return result;
	}
	public:

	/// Calculate the determinant of the matrix
	inline Precision determinant() const {
		Precision result = get_sign();
		for (int i=0; i<my_lu.num_rows(); i++){
			result*=my_lu(i,i);
		}
		return result;
	}
	
	/// Get the LAPACK info
	int get_info() const { return my_info; }

 private:

	Matrix<Size,Size,Precision> my_lu;
	FortranInteger my_info;
	Vector<Size, FortranInteger> my_IPIV;	//Convenient static-or-dynamic array of ints :-)

};
}
	

#endif
