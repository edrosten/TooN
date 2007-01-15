
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
#ifndef __NUMHELPERS_H
#define __NUMHELPERS_H

#include <TooN/TooN.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

//////////// CONTENTS OF THIS FILE ///////////

// normalizations (note US spelling)
template <class Accessor> inline void  normalize(VectorBase<Accessor>& v);
template <class Accessor> inline void  normalize_last(VectorBase<Accessor>& v);
template <class Accessor> inline void  normalize_but_last(VectorBase<Accessor>& v);

// Project
template <int Size, class Accessor> Vector<Size-1>  project(const FixedVector<Size,Accessor>& v);
template  <class Accessor> Vector<>                 project(const DynamicVector<Accessor>& v);

// Unproject
template <int Size, class Accessor> Vector<Size+1>  unproject(const FixedVector<Size,Accessor>& v);
template  <class Accessor> Vector<>                 unproject(const DynamicVector<Accessor>& v);


// as_vector
template<int Size> inline FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(double* data);
template<int Size> inline const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(const double* data);

// set a matrix to (a multiple of) the Identity
template <class Accessor> void  Identity(MatrixBase<Accessor>& M, const double factor=1);
 template <int N> Matrix<N> Identity(const double factor=1) {
     Matrix<N> M;
     Identity(M,factor);
     return M;
 }

// symmetrize a matrix
template <class Accessor> void Symmetrize(MatrixBase<Accessor>& m);

// transpose a matrix
template <class Accessor> void Transpose(MatrixBase<Accessor>& m);

// set a matrix to Zero
template <class Accessor> void Zero(MatrixBase<Accessor>&m);

// set a vector to zero
template <class Accessor> inline void Zero(VectorBase<Accessor>& v);

 template <class T, int N> struct ZeroBlock {
     static const T data[N];
 };

 template <class T, int N> const T ZeroBlock<T,N>::data[N] = {0};
 
 template <int M, int N> inline const Matrix<M,N>& zeros() { return *reinterpret_cast<const Matrix<M,N>*>(ZeroBlock<double,M*N>::data); }
 template <int N> inline const Vector<N>& zeros() { return *reinterpret_cast<const Vector<N>*>(ZeroBlock<double,N>::data); }
 //////////////////////////////////////////////



// normalizations (note US spelling)
template <class Accessor>
inline void normalize(VectorBase<Accessor>& v){
    double sumSq = v[0]*v[0];
    for (int i=1; i<v.size(); i++)
	sumSq += v[i]*v[i];
    double factor = 1.0/sqrt(sumSq);
    for (int i=0; i<v.size(); i++)
	v[i] *= factor;
}

template <class Accessor>
inline void normalize_last(VectorBase<Accessor>& v){
  const double scalefac = 1/v[v.size()-1];
  for(int i=0; i<v.size(); i++){
    v[i]*=scalefac;
  }
}

template <class Accessor>
inline void normalize_but_last(VectorBase<Accessor>& v){
  double sumsq=0;
  for(int i=0; i<v.size()-1; i++){
    sumsq += v[i] * v[i];
  }
  const double scalefac = 1/sqrt(sumsq);
  for(int i=0; i<v.size(); i++){
    v[i]*=scalefac;
  }
}


// Project
template <int Size, class Accessor>
struct FixedVProject {
  inline static void eval(Vector<Size-1>& ret, const FixedVector<Size,Accessor>& v){
    const double scalefac = 1/v[Size-1];
    for(int i=0; i<Size-1; i++){
      ret[i]=v[i]*scalefac;
    }
  }
};

  template <int Size, class Accessor> inline 
Vector<Size-1> project(const FixedVector<Size,Accessor>& v){
  return Vector<Size-1>(v,Operator<FixedVProject<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVProject : public VSizer{
  inline static void eval(Vector<>& ret, const DynamicVector<Accessor>& v){
    const int size = v.size();
    set_size(ret,size-1);
    const double scalefac = 1/v[size-1];
    for(int i=0; i<size-1; i++){
      ret[i]=v[i]*scalefac;
    }
  }
};

template  <class Accessor>
Vector<> project(const DynamicVector<Accessor>& v){
  return Vector<>(v,Operator<DynamicVProject<Accessor> >());
}

// Unproject
template <int Size, class Accessor>
struct FixedVUnproject {
  inline static void eval(Vector<Size+1>& ret, const FixedVector<Size,Accessor>& v){
    ret.template slice<0,Size>() = v;
    ret[Size]=1;
  }
};

template <int Size, class Accessor>
Vector<Size+1> unproject(const FixedVector<Size,Accessor>& v){
  return Vector<Size+1>(v,Operator<FixedVUnproject<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVUnproject : public VSizer{
  inline static void eval(Vector<>& ret, const DynamicVector<Accessor>& v){
    const int size = v.size();
    set_size(ret,size+1);
    v.copy_into(ret.get_data_ptr());
    ret[size]=1;
  }
};

template  <class Accessor>
Vector<> unproject(const DynamicVector<Accessor>& v){
  return Vector<>(v,Operator<DynamicVUnproject<Accessor> >());
}


// as_vector<Size>(double*) to convert a pointer to
// an array of data into a Vector<Size>
template<int Size>
inline FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(double* data){
  return reinterpret_cast<FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&>(*data);
}

template<int Size>
inline const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(const double* data){
  return reinterpret_cast<const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&>(*data);
}


// set a matrix to (a multiple of) the Identity
template <class Accessor>
  void Identity(MatrixBase<Accessor>& M, const double factor){
  assert(M.num_rows() == M.num_cols());
  for(int i=0; i<M.num_rows(); i++){
    for(int j=0; j<M.num_cols(); j++){
      M(i,j)=0;
    }
    M(i,i)=factor;
  }
};


// symmetrize a matrix
template <class Accessor>
void Symmetrize(MatrixBase<Accessor>& m){
  assert(m.num_rows()==m.num_cols());
  for(int r=0; r<m.num_rows()-1; r++){
    for(int c=r; c<m.num_rows(); c++){
	m(c,r) = m(r,c);// = 0.5*(m(r,c)+m(c,r));
    }
  }
}

// Transpose a matrix
template<class Accessor>
void Transpose(MatrixBase<Accessor>& m){
  assert(m.num_rows() == m.num_cols());
  for(int r=0; r<m.num_rows()-1; r++){
    for(int c=r; c<m.num_rows(); c++){
      double temp = m(r,c);
      m(r,c) = m(c,r);
      m(c,r) = temp;
    }
  }
}

// set a Matrix to zero
template <class Accessor> void Zero(MatrixBase<Accessor>&m){
  for(int r=0; r<m.num_rows(); r++){
    for(int c=0; c<m.num_cols(); c++){
      m(r,c) = 0;
    }
  }
}

// set a vector to zero
template <class Accessor> inline void Zero(VectorBase<Accessor>& v){
  for(int i=0; i<v.size(); i++){
    v[i]=0;
  }
}

 namespace util {
     template <class F, int R, int N, class A1, class A2, class A3> inline void transformCovariance(const FixedMatrix<R,N,A1>& A, const FixedMatrix<N,N,A2>& B, FixedMatrix<R,R,A3>& M)
     {
	 for (int i=0; i<R; ++i) {
	     const Vector<N> ABi = B * A[i];
	     F::eval(M[i][i],ABi * A[i]);
	     for (int j=i+1; j<R; ++j) {
		 const double v = ABi * A[j];
		 F::eval(M[i][j], v);
		 F::eval(M[j][i], v);
	     }
	 }
     }

     template <class F, int R, int N, class A1, class A2, class A3> inline void transformCovarianceUpper(const FixedMatrix<R,N,A1>& A, const FixedMatrix<N,N,A2>& B, FixedMatrix<R,R,A3>& M)
     {
	 for (int i=0; i<R; ++i) {
	     const Vector<N> ABi = B * A[i];
	     for (int j=i; j<R; ++j)
		 F::eval(M[i][j], ABi * A[j]);
	 }
     }

     template <class A1, class A2, class MatM> inline void transformCovariance(const DynamicMatrix<A1>& A, const DynamicMatrix<A2>& B, MatM& M)
     {
	 const int R = A.num_rows();
	 const int N = A.num_cols();	 
	 assert(M.num_rows() == R && 
		M.num_cols() == R && 
		B.num_rows() == N && 
		B.num_cols() == N);
	 for (int i=0; i<R; ++i) {
	     const Vector<> ABi = B * A[i];
	     M[i][i] = ABi * A[i];
	     for (int j=i+1; j<R; ++j)
		 M[j][i] = M[i][j] = ABi * A[j];
	 }
     }
 }
 
 template <class A1, class A2> Matrix<> inline transformCovariance(const DynamicMatrix<A1>& A, const DynamicMatrix<A2>& B)
 {
     Matrix<> M(A.num_rows(), A.num_rows());
     util::transformCovariance(A,B,M);
     return M;
 }

 template <int R, int N, class Accessor1, class Accessor2> inline Matrix<R> transformCovariance(const FixedMatrix<R,N,Accessor1>& A, const FixedMatrix<N,N,Accessor2>& B)
 {
     Matrix<R> M;
     util::transformCovariance<util::Assign>(A,B,M);
     return M;
 }

 template <class MatA, class MatB, class MatM> inline void transformCovariance(const MatA& A, const MatB& B, MatM& M)
 {
     util::transformCovariance<util::Assign>(A,B,M);
 }

 template <class F, class MatA, class MatB, class MatM> void transformCovariance(const MatA& A, const MatB& B, MatM& M)
 {
     util::transformCovariance<F>(A,B,M);
 }

 template <int M, int N, int C, class A, class B, class Mat> inline void add_product(const FixedMatrix<M,N,A>& Ma, const FixedMatrix<N,C,B>& Mb, Mat& Mc)
 {
     util::matrix_multiply<util::PlusEquals,M,N,C>(Ma,Mb,Mc);
 }

 template <int M, int N, class A, class B, class Vec> inline void add_product(const FixedMatrix<M,N,A>& m, const FixedVector<N,B>& v, Vec& r)
 {
     util::matrix_multiply<util::PlusEquals,M,N,1>(m,v.as_col(),r.as_col());
 }
 
 template <class A, class B, class C> inline void add_product(const DynamicMatrix<A>& Ma, const DynamicMatrix<B>& Mb, DynamicMatrix<C>& r)
 {
     const int M=Ma.num_rows();
     const int N=Ma.num_cols();
     const int R=Mb.num_cols();
     for (int i=0; i<M; ++i)
	 for (int j=0; j<R; ++j) {
	     double sum = 0;
	     for (int k=0; k<N; ++k)
		 sum += Ma(i,k)*Mb(k,j);
	     r(i,j) += sum;
	 }
 }

 template <int M, int N, int C, class A1, class A2, class Mat> void matrix_multiply(const FixedMatrix<M,N,A1>& A, const FixedMatrix<N,C,A2>& B, Mat& R)
 {
     util::matrix_multiply<M,N,C>(A,B,R);
 }

#ifndef TOON_NO_NAMESPACE
}
#endif 

#endif
