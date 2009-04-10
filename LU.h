// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)
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

#ifndef TOON_INCLUDE_LU_H
#define TOON_INCLUDE_LU_H

#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>

namespace TooN {

template <int Size=-1, class Precision=double>
class LU {
  public:

  template<int S1, int S2, class Base>
  LU(const Matrix<S1,S2,Precision, Base>& m)
  :my_lu(m.num_rows(),m.num_cols()),my_IPIV(m.num_rows()){
    compute(m);
  }
  
  template<int S1, int S2, class Base>
  void compute(const Matrix<S1,S2,Precision,Base>& m){
    //check for consistency with Size
    SizeMismatch<Size, S1>::test(my_lu.num_rows(),m.num_rows());
    SizeMismatch<Size, S2>::test(my_lu.num_rows(),m.num_cols());
	
    //Make a local copy. This is guaranteed contiguous
    my_lu=m;
    int lda = m.num_rows();
    int M = m.num_rows();
    int N = m.num_rows();

    getrf_(&M,&N,&my_lu[0][0],&lda,&my_IPIV[0],&my_info);

    if(my_info < 0){
      std::cerr << "error in LU, INFO was " << my_info << std::endl;
    }
  }

  template <int Rows, int NRHS, class Base>
  Matrix<Size,NRHS,Precision> backsub(const Matrix<Rows,NRHS,Precision,Base>& rhs){
    //Check the number of rows is OK.
    SizeMismatch<Size, Rows>::test(my_lu.num_rows(), rhs.num_rows());
	
    Matrix<Size, NRHS, Precision> result(rhs);

    int M=rhs.num_cols();
    int N=my_lu.num_rows();
    double alpha=1;
    int lda=my_lu.num_rows();
    int ldb=rhs.num_cols();
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
  
  template <int Rows, class Base>
  Vector<Size,Precision> backsub(const Vector<Rows,Precision,Base>& rhs){
    //Check the number of rows is OK.
    SizeMismatch<Size, Rows>::test(my_lu.num_rows(), rhs.size());
	
    Vector<Size, Precision> result(rhs);

    int M=1;
    int N=my_lu.num_rows();
    double alpha=1;
    int lda=my_lu.num_rows();
    int ldb=1;
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

  Matrix<Size,Size,Precision> get_inverse(){
    Matrix<Size,Size,Precision> Inverse(my_lu);
    int N = my_lu.num_rows();
    int lda=my_lu.num_rows();
    int lwork=-1;
    Precision size;
    getri_(&N, &Inverse[0][0], &lda, &my_IPIV[0], &size, &lwork, &my_info);
    lwork=int(size);
    Precision* WORK = new Precision[lwork];
    getri_(&N, &Inverse[0][0], &lda, &my_IPIV[0], WORK, &lwork, &my_info);
    delete [] WORK;
    return Inverse;
  }

  Matrix<Size,Size,Precision>& get_lu(){return my_lu;}
  const Matrix<Size,Size,Precision>& get_lu()const {return my_lu;}

  inline int get_sign() const {
    int result=1;
    for(int i=0; i<my_lu.num_rows()-1; i++){
      if(my_IPIV[i] > i+1){
	result=-result;
      }
    }
    return result;
  }

  inline Precision determinant() const {
    Precision result = get_sign();
    for (int i=0; i<my_lu.num_rows(); i++){
      result*=my_lu(i,i);
    }
    return result;
  }

  int get_info() const { return my_info; }

 private:

  Matrix<Size,Size,Precision> my_lu;
  int my_info;
  Vector<Size, int> my_IPIV;  //Convenient static-or-dynamic array of ints :-)

};
}
  

#endif
