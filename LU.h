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
#ifndef __LU_H
#define __LU_H

#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

template <int Size=-1>
class LU {
  public:

  template<class Accessor>
  LU(const FixedMatrix<Size,Size,Accessor>& m){
    compute(m);
  }

  template<class Accessor>
  void compute(const FixedMatrix<Size,Size,Accessor>& m){
    my_lu=m;
    int lda = Size;
    int M = Size;
    int N = Size;
    dgetrf_(&M,&N,my_lu.get_data_ptr(),&lda,my_IPIV,&my_info);
    if(my_info < 0){
      std::cerr << "error in LU, INFO was " << my_info << std::endl;
    }
  }

  template <int NRHS, class Accessor>
  Matrix<Size,NRHS,RowMajor> backsub(const FixedMatrix<Size,NRHS,Accessor>& rhs){
    Matrix<Size,NRHS,RowMajor> result(rhs);
    int M=NRHS;
    int N=Size;
    double alpha=1;
    int lda=Size;
    int ldb=NRHS;
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);

    // now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      for(int j=0; j<NRHS; j++){
	double temp = result[i][j];
	result[i][j] = result[swaprow][j];
	result[swaprow][j] = temp;
      }
    }
    return result;
  }

  template <class Accessor>
  Matrix<-1,-1,RowMajor> backsub(const DynamicMatrix<Accessor>& rhs){
    Matrix<-1,-1,RowMajor> result(rhs);
    assert(result.num_rows == my_lu.num_rows());
    int M=result.num_cols();
    int N=result.num_rows();
    double alpha=1;
    int lda=result.num_rows();
    int ldb=result.num_cols();
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);

    // now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      for(int j=0; j<result.num_cols(); j++){
	double temp = result(i,j);
	result(i,j) = result(swaprow,j);
	result(swaprow,j) = temp;
      }
    }
    return result;
  }


  template <class Accessor>
  Vector<Size> backsub(const FixedVector<Size,Accessor>& rhs){
    Vector<Size> result(rhs);
    int M=1;
    int N=Size;
    double alpha=1;
    int lda=Size;
    int ldb=1;
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      double temp = result[i];
      result[i]=result[swaprow];
      result[swaprow]=temp;
    }
    return result;
  }

  template <class Accessor>
  Vector<> backsub(const DynamicVector<Accessor>& rhs){
    assert(rhs.size()==Size);
    Vector<> result(rhs);
    int M=1;
    int N=Size;
    double alpha=1;
    int lda=Size;
    int ldb=1;
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      double temp = result[i];
      result[i]=result[swaprow];
      result[swaprow]=temp;
    }
    return result;
  }


  Matrix<Size,Size,RowMajor> get_inverse(){
    Matrix<Size,Size,RowMajor> Inverse = my_lu;
    int N = Size;
    int lda=Size;
    int lwork=-1;
    double size;
    dgetri_(&N, Inverse.get_data_ptr(), &lda, my_IPIV, &size, &lwork, &my_info);
    lwork=int(size);
    double* WORK = new double[lwork];
    dgetri_(&N, Inverse.get_data_ptr(), &lda, my_IPIV, WORK, &lwork, &my_info);
    delete [] WORK;
    return Inverse;
  }

  Matrix<Size,Size,RowMajor>& get_lu(){return my_lu;}
  const Matrix<Size,Size,RowMajor>& get_lu()const {return my_lu;}

  inline int get_sign() const {
    int result=1;
    for(int i=0; i<Size-1; i++){
      if(my_IPIV[i] > i+1){
	result=-result;
      }
    }
    return result;
  }

  inline double determinant() const {
    double result = get_sign();
    for (int i=0; i<Size; i++){
      result*=my_lu(i,i);
    }
    return result;
  }

  int get_info() const { return my_info; }
 private:
  Matrix<Size,Size,RowMajor> my_lu;
  int my_info;
  int my_IPIV[Size];
};
  

template <>
class LU<> {
  public:

  LU(int size):my_lu(size,size){
    my_IPIV = new int[size];
  }

  template<class Accessor>
  LU(const DynamicMatrix<Accessor>& m) :
    my_lu(m.num_rows(),m.num_cols())
  {
    my_IPIV = new int[m.num_rows()];
    assert(m.num_rows() == m.num_cols());
    compute(m);
  }

  ~LU(){delete[] my_IPIV;}


  template<class Accessor>
  void compute(const DynamicMatrix<Accessor>& m){
    my_lu=m;
    int lda = my_lu.num_cols();
    int M = lda;
    int N = lda;
    dgetrf_(&M,&N,my_lu.get_data_ptr(),&lda,my_IPIV,&my_info);
    if(my_info < 0){
      std::cerr << "error in LU, INFO was " << my_info << std::endl;
    }
  }


  template <int Rows, int Cols, class Accessor>
  Matrix<Rows,Cols,RowMajor> backsub(const FixedMatrix<Rows,Cols,Accessor>& rhs){
    assert(my_lu.num_rows() == rhs.num_rows());
    Matrix<Rows,Cols,RowMajor> result(rhs);
    int M=result.num_cols();
    int N=result.num_rows();
    double alpha=1;
    int lda=my_lu.num_rows();
    int ldb=result.num_cols();
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);

    // now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      for(int j=0; j<result.num_cols(); j++){
	double temp = result(i,j);
	result(i,j) = result(swaprow,j);
	result(swaprow,j) = temp;
      }
    }
    return result;
  }


  template <class Accessor>
  Matrix<-1,-1,RowMajor> backsub(const DynamicMatrix<Accessor>& rhs){
    assert(my_lu.num_rows() == rhs.num_rows());
    Matrix<-1,-1,RowMajor> result(rhs);
    int M=result.num_cols();
    int N=result.num_rows();
    double alpha=1;
    int lda=my_lu.num_rows();
    int ldb=result.num_cols();
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);

    // now do the row swapping (lapack dlaswp.f only shuffles fortran rows = Rowmajor cols)
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      for(int j=0; j<result.num_cols(); j++){
	double temp = result(i,j);
	result(i,j) = result(swaprow,j);
	result(swaprow,j) = temp;
      }
    }
    return result;
  }

  template <int Size, class Accessor>
  Vector<Size> backsub(const FixedVector<Size,Accessor>& rhs){
    assert(rhs.size() == my_lu.num_rows());
    Vector<Size> result(rhs);
    int M=1;
    int N=result.size();
    double alpha=1;
    int lda=result.size();
    int ldb=1;
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      double temp = result[i];
      result[i]=result[swaprow];
      result[swaprow]=temp;
    }
    return result;
  }

  
  template <class Accessor>
  Vector<> backsub(const DynamicVector<Accessor>& rhs){
    assert(rhs.size() == my_lu.num_rows());
    Vector<> result(rhs);
    int M=1;
    int N=result.size();
    double alpha=1;
    int lda=result.size();
    int ldb=1;
    dtrsm_("R","U","N","N",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    dtrsm_("R","L","N","U",&M,&N,&alpha,my_lu.get_data_ptr(),&lda,result.get_data_ptr(),&ldb);
    for(int i=N-1; i>=0; i--){
      const int swaprow = my_IPIV[i]-1; // fortran arrays start at 1
      double temp = result[i];
      result[i]=result[swaprow];
      result[swaprow]=temp;
    }
    return result;
  }

  Matrix<-1,-1,RowMajor> get_inverse(){
    Matrix<-1,-1,RowMajor> Inverse(my_lu);
    int N = my_lu.num_rows();
    int lda=my_lu.num_rows();
    int lwork=-1;
    double size;
    dgetri_(&N, Inverse.get_data_ptr(), &lda, my_IPIV, &size, &lwork, &my_info);
    lwork=int(size);
    double* WORK = new double[lwork];
    dgetri_(&N, Inverse.get_data_ptr(), &lda, my_IPIV, WORK, &lwork, &my_info);
    delete [] WORK;
    return Inverse;
  }

  Matrix<-1,-1,RowMajor>& get_lu(){return my_lu;}
  const Matrix<-1,-1,RowMajor>& get_lu()const {return my_lu;}

  inline int get_sign() const {
    int result=1;
    for(int i=0; i<my_lu.num_rows()-1; i++){
      if(my_IPIV[i] > i+1){
	result=-result;
      }
    }
    return result;
  }

  inline double determinant() const {
    double result = get_sign();
    for (int i=0; i<my_lu.num_rows(); i++){
      result*=my_lu(i,i);
    }
    return result;
  }

  int get_info() const { return my_info; }


 private:
  Matrix<-1,-1,RowMajor> my_lu;
  int my_info;
  int* my_IPIV;
};



#ifndef TOON_NO_NAMESPACE
}
#endif 





#endif
