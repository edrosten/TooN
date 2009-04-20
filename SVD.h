// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk)
//
// This file is part of the TooN Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

#ifndef __SVD_H
#define __SVD_H

#include <TooN/TooN.h>
#include <TooN/lapack.h>

namespace TooN {

	// TODO - should this depend on precision?
static const double condition_no=1e9; // GK HACK TO GLOBAL

// Use <-1> template for TooN2 SVD

template<int Rows=Dynamic, int Cols=Rows, typename Precision=DefaultPrecision>
class SVD {
public:
	// this is the size of the diagonal
	// NB works for semi-dynamic sizes because -1 < +ve ints
	static const int Min_Dim = Rows<Cols?Rows:Cols;
	
	/// default constructor for Rows>0 and Cols>0
	SVD() {}

	/// constructor for Rows=-1 or Cols=-1 (or both)
	SVD(int rows, int cols)
		: my_copy(rows,cols),
		  my_diagonal(std::min(rows,cols)),
		  my_square(std::min(rows,cols), std::min(rows,cols))
	{}

	
	template <int R2, int C2, typename P2, typename B2>
	SVD(const Matrix<R2,C2,P2,B2>& m)
		: my_copy(m),
		  my_diagonal(std::min(m.num_rows(),m.num_cols())),
		  my_square(std::min(m.num_rows(),m.num_cols()),std::min(m.num_rows(),m.num_cols()))
	{
		do_compute();
	}

	template <int R2, int C2, typename P2, typename B2>
	void compute(const Matrix<R2,C2,P2,B2>& m){
		my_copy=m;
		do_compute();
	}
		
	void do_compute(){
		Precision* const a = my_copy.my_data;
		int lda = my_copy.num_cols();
		int m = my_copy.num_cols();
		int n = my_copy.num_rows();
		Precision* const uorvt = my_square.my_data;
		Precision* const s = my_diagonal.my_data;
		int ldu;
		int ldvt = lda;
		int LWORK;
		int INFO;
		char JOBU;
		char JOBVT;

		if(is_vertical()){ // u is a
			JOBU='O';
			JOBVT='S';
			ldu = lda;
		} else { // vt is a
			JOBU='S';
			JOBVT='O';
			ldu = my_square.num_cols();
		}

		Precision* wk;

		Precision size;
		LWORK = -1;

		// arguments are scrambled because we use rowmajor and lapack uses colmajor
		// thus u and vt play each other's roles.
		dgesvd_( &JOBVT, &JOBU, &m, &n, a, &lda, s, uorvt,
				 &ldvt, uorvt, &ldu, &size, &LWORK, &INFO);
	
		LWORK = (long int)(size);
		wk = new Precision[LWORK];

		dgesvd_( &JOBVT, &JOBU, &m, &n, a, &lda, s, uorvt,
				 &ldvt, uorvt, &ldu, wk, &LWORK, &INFO);
	
		delete[] wk;
	}

	bool is_vertical(){ return (my_copy.num_rows() >= my_copy.num_cols()); }

	int min_dim(){ return std::min(my_copy.num_rows(), my_copy.num_cols()); }

	Matrix<Rows,Min_Dim,Precision,RowMajor>& get_U(){
		if(is_vertical()){
			return my_copy;
		} else {
			return my_square;
		}
	}

	Vector<Min_Dim,Precision>& get_diagonal(){ return my_diagonal; }

	Matrix<Min_Dim,Cols,Precision,RowMajor>& get_VT(){
		if(is_vertical()){
			return my_square;
		} else {
			return my_copy;
		}
	}


	template <int Rows2, int Cols2, typename P2, typename B2>
	Matrix<Cols,Cols2, typename Internal::MultiplyType<Precision,P2>::type >
	backsub(const Matrix<Rows2,Cols2,P2,B2>& rhs, const Precision condition=condition_no)
	{
		Vector<Min_Dim> inv_diag(min_dim());
		get_inv_diag(inv_diag,condition);
		return (get_VT().T() * diagmult(inv_diag, (get_U().T() * rhs)));
	}

	template <int Size, typename P2, typename B2>
	Vector<Cols, typename Internal::MultiplyType<Precision,P2>::type >
	backsub(const Vector<Size,P2,B2>& rhs, const Precision condition=condition_no)
	{
		Vector<Min_Dim> inv_diag(min_dim());
		get_inv_diag(inv_diag,condition);
		return (get_VT().T() * diagmult(inv_diag, (get_U().T() * rhs)));
	}

	Matrix<Cols,Rows> get_pinv(const Precision condition = condition_no){
		Vector<Min_Dim> inv_diag(min_dim());
		get_inv_diag(inv_diag,condition);
		return diagmult(get_VT().T(),inv_diag) * get_U().T();
	}

	/// calculates the product of the singular values
	/// for square matrices this is the determinant
	Precision determinant() {
		Precision result = my_diagonal[0];
		for(int i=1; i<my_diagonal.size(); i++){
			result *= my_diagonal[i];
		}
		return result;
	}
	
	int rank(const Precision condition = condition_no) {
		if (my_diagonal[0] == 0) return 0;
		int result=1;
		for(int i=0; i<min_dim(); i++){
			if(my_diagonal[i] * condition <= my_diagonal[0]){
				result++;
			}
		}
		return result;
	}

	void get_inv_diag(Vector<Min_Dim>& inv_diag, const Precision condition){
		for(int i=0; i<min_dim(); i++){
			if(my_diagonal[i] * condition <= my_diagonal[0]){
				inv_diag[i]=0;
			} else {
				inv_diag[i]=static_cast<Precision>(1)/my_diagonal[i];
			}
		}
	}

private:
	Matrix<Rows,Cols,Precision,RowMajor> my_copy;
	Vector<Min_Dim,Precision> my_diagonal;
	Matrix<Min_Dim,Min_Dim,Precision,RowMajor> my_square; // square matrix (U or V' depending on the shape of my_copy)
};






/// version of SVD forced to be square
/// princiapally here to allow use in WLS
template<int Size, typename Precision>
struct SQSVD : public SVD<Size, Size, Precision> {
	// forward all constructors to SVD
	SQSVD() {}
	SQSVD(int size) : SVD<Size,Size,Precision>(size, size) {}
	
	template <int R2, int C2, typename P2, typename B2>
	SQSVD(const Matrix<R2,C2,P2,B2>& m) : SVD<Size,Size,Precision>(m) {}
};


}


#endif
