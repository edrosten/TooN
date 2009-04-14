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

#ifndef __SYMEIGEN_H
#define __SYMEIGEN_H

#include <iostream>
#include <cassert>
#include <TooN/lapack.h>

#include <TooN/TooN.h>

namespace TooN {

static const double symeigen_condition_no=1e9;

template <int Size> struct ComputeSymEigen {

	template<typename P, typename B>
	static inline void compute(const Matrix<Size,Size,P, B>& m, Matrix<Size,Size,P> & evectors, Vector<Size, P>& evalues) {
		evectors = m;
		int N = evalues.size();
		int lda = evalues.size();
		int info;
		int lwork=-1;
		P size;

		// find out how much space fortran needs
		syev_((char*)"V",(char*)"U",&N,&evectors[0][0],&lda,&evalues[0], &size,&lwork,&info);
		lwork = int(size);
		Vector<Dynamic, P> WORK(lwork);

		// now compute the decomposition
		syev_((char*)"V",(char*)"U",&N,&evectors[0][0],&lda,&evalues[0], &WORK[0],&lwork,&info);

		if(info!=0){
			std::cerr << "In SymEigen<"<<Size<<">: " << info 
					<< " off-diagonal elements of an intermediate tridiagonal form did not converge to zero." << std::endl
					<< "M = " << m << std::endl;
		}
	}
};

template <> struct ComputeSymEigen<2> {

	template<typename P, typename B>
	static inline void compute(const Matrix<2,2,P,B>& m, Matrix<2,2,P>& eig, Vector<2, P>& ev) {
		double trace = m[0][0] + m[1][1];
		double det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
		double disc = trace*trace - 4 * det;
		assert(disc>=0);
		double root_disc = sqrt(disc);
		ev[0] = 0.5 * (trace - root_disc);
		ev[1] = 0.5 * (trace + root_disc);
		double a = m[0][0] - ev[0];
		double b = m[0][1];
		double magsq = a*a + b*b;
		if (magsq == 0) {
			eig[0][0] = 1.0;
			eig[0][1] = 0;
		} else {
			eig[0][0] = -b;
			eig[0][1] = a;
			eig[0] *= 1.0/sqrt(magsq);
		}
		eig[1][0] = -eig[0][1];
		eig[1][1] = eig[0][0];
	}
};
    
template <int Size=Dynamic, typename Precision = double>
class SymEigen {
public:
  inline SymEigen(){}

  template<int R, int C, typename B>
  inline SymEigen(const Matrix<R, C, Precision, B>& m) : my_evectors(m.num_rows(), m.num_cols()), my_evalues(m.num_rows()) {
    compute(m);
  }

  template<int R, int C, typename B>
  inline void compute(const Matrix<R,C,Precision,B>& m){
		SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());
		SizeMismatch<R, Size>::test(m.num_rows(), my_evectors.num_rows());
		ComputeSymEigen<Size>::compute(m, my_evectors, my_evalues);
  }
  
	template <int S, typename P, typename B>
	Vector<Size, Precision> backsub(const Vector<S,P,B>& rhs) const {
		return (my_evectors.T() * diagmult(get_inv_diag(symeigen_condition_no),(my_evectors * rhs)));
	}

	template <int R, int C, typename P, typename B>
	Matrix<Size,C, Precision> backsub(const Matrix<R,C,P,B>& rhs) const {
		return (my_evectors.T() * diagmult(get_inv_diag(symeigen_condition_no),(my_evectors * rhs)));
	}

	Matrix<Size, Size, Precision> get_pinv(const double condition=symeigen_condition_no) const {
		return my_evectors.T() * diagmult(get_inv_diag(condition),my_evectors);
	}

	Vector<Size, Precision> get_inv_diag(const double condition) const {
		Precision max_diag = -my_evalues[0] > my_evalues[my_evalues.size()-1] ? -my_evalues[0]:my_evalues[my_evalues.size()-1];
		Vector<Size, Precision> invdiag(my_evalues.size());
		for(int i=0; i<my_evalues.size(); i++){
			if(fabs(my_evalues[i]) * condition > max_diag) {
				invdiag[i] = 1/my_evalues[i];
			} else {
				invdiag[i]=0;
			}
		}
		return invdiag;
	}
  
	Matrix<Size,Size,Precision>& get_evectors() {return my_evectors;}
	const Matrix<Size,Size,Precision>& get_evectors() const {return my_evectors;}
	Vector<Size, Precision>& get_evalues() {return my_evalues;}
	const Vector<Size, Precision>& get_evalues() const {return my_evalues;}
  
	bool is_posdef() const {
	  for (int i = 0; i < my_evalues.size(); ++i) {
		  if (my_evalues[i] <= 0.0)
			  return false;
	  }
	  return true;
	}
  
	bool is_negdef() const {
	  for (int i = 0; i < my_evalues.size(); ++i) {
		  if (my_evalues[i] >= 0.0)
			  return false;
	  }
	  return true;
	}

	Precision get_determinant () const {
	  Precision det = 1.0;
	  for (int i = 0; i < my_evalues.size(); ++i) {
		  det *= my_evalues[i];
	  }
	  return det;
	}

private:
	// eigen vectors laid out row-wise so evectors[i] is the ith evector
	Matrix<Size,Size,Precision> my_evectors;

	Vector<Size, Precision> my_evalues;
};

}

#endif
