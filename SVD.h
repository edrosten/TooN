
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
#ifndef __SVD_H
#define __SVD_H

#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

enum{Horizontal,Vertical};

template <int Rows, int Cols, int HV>
class HV_SVD;

// (for Rows > Cols)
template <int Rows, int Cols>
class HV_SVD <Rows,Cols,Vertical> {
 public:
  static const int Inter = (Rows>Cols?Cols:Rows);

  HV_SVD(){}

  template<class Accessor>
  HV_SVD(const FixedMatrix<Rows,Cols,Accessor>& m){
    compute(m);
  }

  template<class Accessor>
  void compute(const FixedMatrix<Rows,Cols,Accessor>& m){
    my_U=m;
    int width=Cols;
    int height=Rows;
    int lwork=-1;
    int info;
    double size;
    // find out FORTRAN space requirements
    dgesvd_(const_cast<char*>("S"),const_cast<char*>("O"),&width,&height,my_U.get_data_ptr(),&width,my_diagonal.get_data_ptr(),my_VT.get_data_ptr(),&width,my_U.get_data_ptr(),&width,&size,&lwork,&info);
    lwork = int(size);
    double * WORK = new double[lwork];
    dgesvd_(const_cast<char*>("S"),const_cast<char*>("O"),&width,&height,my_U.get_data_ptr(),&width,my_diagonal.get_data_ptr(),my_VT.get_data_ptr(),&width,my_U.get_data_ptr(),&width,WORK,&lwork,&info);
    delete [] WORK;
    if(info!=0){
      std::cerr << "error - info was " << info << std::endl;
    }
  }

  Matrix<Rows,Inter,RowMajor>& get_U(){return my_U;}
  Vector<Inter>& get_diagonal(){return my_diagonal;}
  Matrix<Inter,Cols,RowMajor>& get_VT(){return my_VT;}

 protected:
  Matrix<Rows,Inter,RowMajor> my_U;
  Vector<Inter> my_diagonal;
  Matrix<Inter,Cols,RowMajor> my_VT;
};


// (for Rows < Cols)
template <int Rows, int Cols>
class HV_SVD <Rows,Cols,Horizontal> {
public:
  static const int Inter = (Rows>Cols?Cols:Rows);

  HV_SVD(){}

  template<class Accessor>
  HV_SVD(const FixedMatrix<Rows,Cols,Accessor>& m){
    compute(m);
  }

  template<class Accessor>
  void compute(const FixedMatrix<Rows,Cols,Accessor>& m){
    my_VT=m;
    int width=Cols;
    int height=Rows;
    int lwork=-1;
    int info;
    double size;
    // find out FORTRAN space requirements
    dgesvd_(const_cast<char*>("O"),const_cast<char*>("S"),&width,&height,my_VT.get_data_ptr(),&width,my_diagonal.get_data_ptr(),my_VT.get_data_ptr(),&width,my_U.get_data_ptr(),&height,&size,&lwork,&info);
    lwork = int(size);
    double * WORK = new double[lwork];
    dgesvd_(const_cast<char*>("O"),const_cast<char*>("S"),&width,&height,my_VT.get_data_ptr(),&width,my_diagonal.get_data_ptr(),my_VT.get_data_ptr(),&width,my_U.get_data_ptr(),&height,WORK,&lwork,&info);
    delete [] WORK;
    if(info!=0){
      std::cerr << "error - info was " << info << std::endl;
    }
  }

  Matrix<Rows,Inter,RowMajor>& get_U(){return my_U;}
  Vector<Inter>& get_diagonal(){return my_diagonal;}
  Matrix<Inter,Cols,RowMajor>& get_VT(){return my_VT;}

 protected:
  Matrix<Rows,Inter> my_U;
  Vector<Inter> my_diagonal;
  Matrix<Inter,Cols> my_VT;
};


static const double condition_no=1e9; // GK HACK TO GLOBAL

template <int Rows=-1, int Cols=Rows>
class SVD : public HV_SVD< Rows, Cols, (Rows<Cols?Horizontal:Vertical) > {
  public:
  SVD(){}



  template<class Accessor>
  SVD(const FixedMatrix<Rows,Cols,Accessor>& m): HV_SVD< Rows, Cols, (Rows<Cols?Horizontal:Vertical) > (m) {}

  template <int RHS, class Accessor>
    Matrix<Cols,RHS> backsub(const FixedMatrix<Rows,RHS,Accessor>& rhs, const double condition=condition_no){
    get_inv_diag(condition);
    return (this->my_VT.T() * diagmult(this->inv_diag, (this->my_U.T() * rhs)));
  }

  template<class Accessor>
    Matrix<> backsub(const DynamicMatrix<Accessor>& rhs, const double condition=condition_no){
    get_inv_diag(condition);
    return (this->my_VT.T() * diagmult(this->inv_diag, (this->my_U.T() * rhs)));
  }

  template <class Accessor>
    Vector<Cols> backsub(const FixedVector<Rows,Accessor>& v, const double condition=condition_no){
    get_inv_diag(condition);
    return (this->my_VT.T() * diagmult(this->inv_diag, (this->my_U.T() * v)));
  }

  template <class Accessor>
    Vector<> backsub(const DynamicVector<Accessor>& v, const double condition=condition_no){
    get_inv_diag(condition);
    return (this->my_VT.T() * diagmult(this->inv_diag, (this->my_U.T() * v)));
  }
 

  Matrix<Cols,Rows> get_pinv(const double condition = condition_no){
    get_inv_diag(condition);
    return diagmult(this->my_VT.T(),this->inv_diag) * this->my_U.T();
  }


  void get_inv_diag(const double condition){
    for(int i=0; i<HV_SVD< Rows, Cols, (Rows<Cols?Horizontal:Vertical) >::Inter; i++){
      if(this->my_diagonal[i] * condition <= this->my_diagonal[0]){
	inv_diag[i]=0;
      } else {
	inv_diag[i]=1.0/this->my_diagonal[i];
      }
    }
  }

 private:

  Vector<HV_SVD< Rows, Cols, (Rows<Cols?Horizontal:Vertical) >::Inter> inv_diag;

};


template<>
class SVD<-1> {
public:
  template <class Accessor>
  SVD(const MatrixBase<Accessor>& m) : my_orig(m),
				       my_height(m.num_rows()),
				       my_width(m.num_cols()),
				       my_min_dim(my_height<my_width?my_height:my_width),
				       my_diagonal(my_min_dim),
				       my_square(my_min_dim,my_min_dim) {
    int lwork=-1;
    int info;
    double size;

    // compute the storage requirements
    if(is_vertical()){
      dgesvd_(const_cast<char*>("S"),const_cast<char*>("O"),&my_width,&my_height,my_orig.get_data_ptr(),&my_width,my_diagonal.get_data_ptr(),my_square.get_data_ptr(),&my_width,my_orig.get_data_ptr(),&my_width,&size,&lwork,&info);
    } else {
      dgesvd_(const_cast<char*>("O"),const_cast<char*>("S"),&my_width,&my_height,my_orig.get_data_ptr(),&my_width,my_diagonal.get_data_ptr(),my_orig.get_data_ptr(),&my_width,my_square.get_data_ptr(),&my_height,&size,&lwork,&info);
    }
    lwork = int(size);
    double * WORK = new double[lwork];
    if(is_vertical()){
      dgesvd_(const_cast<char*>("S"),const_cast<char*>("O"),&my_width,&my_height,my_orig.get_data_ptr(),&my_width,my_diagonal.get_data_ptr(),my_square.get_data_ptr(),&my_width,my_orig.get_data_ptr(),&my_width,WORK,&lwork,&info);
    } else {
      dgesvd_(const_cast<char*>("O"),const_cast<char*>("S"),&my_width,&my_height,my_orig.get_data_ptr(),&my_width,my_diagonal.get_data_ptr(),my_orig.get_data_ptr(),&my_width,my_square.get_data_ptr(),&my_height,WORK,&lwork,&info);
    }
    delete [] WORK;
    if(info!=0){
      std::cerr << "error - info was " << info << std::endl;
    }
  }

  bool is_vertical(){return (my_orig.num_rows() >= my_orig.num_cols());}

  Matrix<-1,-1,RowMajor>& get_U(){if(is_vertical()){return my_orig;} else {return my_square;}}
  Vector<-1>& get_diagonal(){return my_diagonal;}
  Matrix<-1,-1,RowMajor>& get_VT(){if(is_vertical()){return my_square;} else {return my_orig;}}


  template <class Accessor>
    Matrix<> backsub(const DynamicMatrix<Accessor>& rhs, const double condition=condition_no){
    Vector<> inv_diag(my_min_dim);
    get_inv_diag(inv_diag,condition);
    return (get_VT().T() * diagmult(inv_diag, (get_U().T() * rhs)));
  }

  template <int R, int C, class Accessor>
    Matrix<R,C> backsub(const FixedMatrix<R, C, Accessor>& rhs, const double condition=condition_no){
    Vector<> inv_diag(my_min_dim);
    get_inv_diag(inv_diag,condition);
    return (get_VT().T() * diagmult(inv_diag, (get_U().T() * rhs)));
  }

  template <class Accessor>
    Vector<> backsub(const DynamicVector<Accessor>& v, const double condition=condition_no){
    Vector<> inv_diag(my_min_dim);
    get_inv_diag(inv_diag,condition);
    return (get_VT().T() * diagmult(inv_diag, (get_U().T() * v)));
  }

  template <int Size, class Accessor>
    Vector<Size> backsub(const FixedVector<Size, Accessor>& v, const double condition=condition_no){
    Vector<> inv_diag(my_min_dim);
    get_inv_diag(inv_diag,condition);
    return (get_VT().T() * diagmult(inv_diag, (get_U().T() * v)));
  }

  Matrix<> get_pinv(const double condition = condition_no){
    Vector<> inv_diag(my_min_dim);
    get_inv_diag(inv_diag,condition);
    return diagmult(get_VT().T(),inv_diag) * get_U().T();
  }


  void get_inv_diag(Vector<>& inv_diag, const double condition){
    for(int i=0; i<my_min_dim; i++){
      if(my_diagonal[i] * condition <= my_diagonal[0]){
	inv_diag[i]=0;
      } else {
	inv_diag[i]=1.0/my_diagonal[i];
      }
    }
  }



private:
  Matrix<-1,-1,RowMajor> my_orig;  // matrix with the original shape
  int my_height;
  int my_width;
  int my_min_dim;
  Vector<-1> my_diagonal;
  Matrix<-1,-1,RowMajor> my_square;   // square matrix (U or V' depending on the shape of my_orig)
};




#ifndef TOON_NO_NAMESPACE
}
#endif 








#endif
