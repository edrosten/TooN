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


#ifndef __SYMEIGEN_H
#define __SYMEIGEN_H

#include <iostream>
#include <lapack.h>

#include <TooN/toon.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

static const double symeigen_condition_no=1e9;

template <int Size=-1>
class SymEigen {
public:
  inline SymEigen(){}


  template<class Accessor>
  inline SymEigen(const FixedMatrix<Size,Size,Accessor>& m){
    compute(m);
  }

  template<class Accessor>
  inline void compute(const FixedMatrix<Size,Size,Accessor>& m){
    my_evectors = m;
    int N = Size;
    int lda = Size;
    int info;
    int lwork=-1;
    double size;

    // find out how much space fortran needs
    dsyev_("V","U",&N,my_evectors.get_data_ptr(),&lda,my_evalues.get_data_ptr(),
	   &size,&lwork,&info);
    lwork = int(size);
    double* WORK = new double[lwork];

    // now compute the decomposition
    dsyev_("V","U",&N,my_evectors.get_data_ptr(),&lda,my_evalues.get_data_ptr(),
	   WORK,&lwork,&info);
    delete [] WORK;
    if(info!=0){
      std::cerr << "info was not zero in SymEigen - it was " << info << std::endl;
    }
  }
  
  template <class Accessor>
  Vector<Size> backsub(const FixedVector<Size,Accessor>& rhs){
    Vector<Size> invdiag;
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }

  template <class Accessor>
  Vector<> backsub(const DynamicVector<Accessor>& rhs){
    Vector<Size> invdiag;
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }


  template <int NRHS, class Accessor>
  Matrix<Size,NRHS> backsub(const FixedMatrix<Size,NRHS,Accessor>& rhs){
    Vector<Size> invdiag;
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }

  template <class Accessor>
  Matrix<> backsub(const DynamicMatrix<Accessor>& rhs){
    Vector<Size> invdiag;
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }


  Matrix<Size> get_pinv(const double condition=symeigen_condition_no){
    Vector<Size> invdiag;
    get_inv_diag(invdiag,condition);
    return my_evectors.T() * diagmult(invdiag,my_evectors);
  }


  void get_inv_diag(Vector<Size>& invdiag, double condition){
    double max_diag = -my_evalues[0] > my_evalues[Size-1] ? -my_evalues[0]:my_evalues[Size-1];
    for(int i=0; i<Size; i++){
      if(fabs(my_evalues[i]) * condition > max_diag) {
	invdiag[i] = 1/my_evalues[i];
      } else {
	invdiag[i]=0;
      }
    }
  }
  
  inline Matrix<Size,Size,RowMajor>& get_evectors() {return my_evectors;}
  inline const Matrix<Size,Size,RowMajor>& get_evectors() const {return my_evectors;}
  inline Vector<Size>& get_evalues() {return my_evalues;}
  inline const Vector<Size>& get_evalues() const {return my_evalues;}

private:
  // eigen vectors laid out row-wise so evectors[i] is the ith evector
  Matrix<Size,Size,RowMajor> my_evectors;

  Vector<Size> my_evalues;
};


template <>
class SymEigen<> {
public:

  inline SymEigen(int size) : 
    my_evectors(size,size),
    my_evalues(size) {}

  template<class Accessor>
  inline SymEigen(const DynamicMatrix<Accessor>& m) : 
    my_evectors(m.num_rows(), m.num_cols()),
    my_evalues(m.num_rows())
  {
    assert(m.num_rows() == m.num_cols());
    compute(m);
  }

  template<class Accessor>
  inline void compute(const DynamicMatrix<Accessor>& m){
    my_evectors = m;
    int N = m.num_cols();
    int lda = m.num_cols();
    int info;
    int lwork=-1;
    double size;

    // find out how much space fortran needs
    dsyev_("V","U",&N,my_evectors.get_data_ptr(),&lda,my_evalues.get_data_ptr(),
	   &size,&lwork,&info);
    lwork = int(size);
    double* WORK = new double[lwork];

    // now compute the decomposition
    dsyev_("V","U",&N,my_evectors.get_data_ptr(),&lda,my_evalues.get_data_ptr(),
	   WORK,&lwork,&info);
    delete [] WORK;
    if(info!=0){
      std::cerr << "info was not zero in SymEigen - it was " << info << std::endl;
    }
  }
  
  template <int Size, class Accessor>
  Vector<Size> backsub(const FixedVector<Size,Accessor>& rhs){
    Vector<> invdiag(my_evalues.size());
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }

  template <class Accessor>
  Vector<> backsub(const DynamicVector<Accessor>& rhs){
    Vector<> invdiag(my_evalues.size());
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }


  template <int NRHS, class Accessor, int Size>
  Matrix<Size,NRHS> backsub(const FixedMatrix<Size,NRHS,Accessor>& rhs){
    Vector<> invdiag(my_evalues.size());
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }

  template <class Accessor>
  Matrix<> backsub(const DynamicMatrix<Accessor>& rhs){
    Vector<> invdiag(my_evalues.size());
    get_inv_diag(invdiag,symeigen_condition_no);
    return (my_evectors.T() * diagmult(invdiag,(my_evectors * rhs)));
  }




  Matrix<> get_pinv(const double condition=symeigen_condition_no){
    Vector<> invdiag(my_evalues.size());
    get_inv_diag(invdiag,condition);
    return my_evectors.T() * diagmult(invdiag,my_evectors);
  }


  void get_inv_diag(Vector<>& invdiag, double condition){
    double max_diag = -my_evalues[0] > my_evalues[invdiag.size()-1] ? -my_evalues[0]:my_evalues[invdiag.size()-1];
    for(int i=0; i<invdiag.size(); i++){
      if(fabs(my_evalues[i]) * condition > max_diag) {
	invdiag[i] = 1/my_evalues[i];
      } else {
	invdiag[i]=0;
      }
    }
  }
  

  
  inline Matrix<-1,-1,RowMajor>& get_evectors() {return my_evectors;}
  inline const Matrix<-1,-1,RowMajor>& get_evectors() const {return my_evectors;}
  inline Vector<>& get_evalues() {return my_evalues;}
  inline const Vector<>& get_evalues() const {return my_evalues;}

private:
  // eigen vectors laid out row-wise so evectors[i] is the ith evector
  Matrix<-1,-1,RowMajor> my_evectors;
  Vector<> my_evalues;
};


#ifndef TOON_NO_NAMESPACE
}
#endif 



#endif
