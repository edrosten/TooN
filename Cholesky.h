// -*- c++ -*-

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


#ifndef __CHOLESKY_H
#define __CHOLESKY_H

#include <iostream>

#include <TooN/lapack.h>

#include <TooN/TooN.h>
#include <TooN/util.h>
#include <limits>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 
    namespace util {
	template <int Size, class A1, class A2, class A3, class A4> inline
	void cholesky_backsub(const FixedMatrix<Size,Size,A1>& L, const FixedVector<Size,A2>& invdiag, const FixedVector<Size,A3>& v, FixedVector<Size,A4>& x) 
	{
	    Vector<Size> t;
	    // forward substitution
	    for (int i=0; i<Size; i++) {
		double b = v[i];
		for (int j=0; j<i;j++)
		    b -= L[i][j]*t[j];
		t[i] = b*invdiag[i];
	    }
	    // back substitution
	    for (int i=Size-1; i >=0; i--) {
		double b = t[i];
		for (int j=i+1; j<Size;j++)
		    b -= L[j][i]*x[j];
		x[i] = b*invdiag[i];
	    }
	}
	
	template <int Size, class A1, class A2, class A3> inline int cholesky_compute(const FixedMatrix<Size,Size,A1>& M, FixedMatrix<Size,Size,A2>& L, FixedVector<Size,A3>& invdiag) 
	{
	    int rank = Size;
	    const double eps = std::numeric_limits<double>::min();
	    for(int i=0;i<Size;i++) {
		double a=M[i][i];
		for(int k=0;k<i;k++) 
		    a-=L[i][k]*L[i][k];
		if (a < eps) {
		    --rank;
		    L[i][i] = invdiag[i] = 0;
		} else {
		    L[i][i]=sqrt(a);
		    invdiag[i] = 1.0/L[i][i];
		}
		for(int j=i+1;j<Size;j++) {
		    L[i][j] = 0;
		    a=M[i][j];
		    for(int k=0;k<i;k++) 
			a-=L[i][k]*L[j][k];
		    L[j][i]=a*invdiag[i];
		}
	    }
	    return rank;
	}

	template <int S, class A1, class A2, class A3> inline void cholesky_inverse(const FixedMatrix<S,S,A1>& L, const FixedVector<S,A2>& invdiag, FixedMatrix<S,S,A3>& I)
	{
	    for (int col = 0; col<S; col++) {
		Vector<S> t,x;
		for (int i=col; i<S; i++) {
		    double psum = 1;
		    for (int j=col; j<i; j++)
			psum -= L[i][j]*t[j];		    
		    t[i] = invdiag[i] * psum;
		}
		for (int i=S-1; i>col; i--) {
		    double psum = t[i];
		    for (int j=i+1; j<S; j++)
			psum -= L[j][i]*x[j];
		    I[i][col] = I[col][i] = x[i] = invdiag[i] * psum;
		}
		double psum = t[col];
		for (int j=col+1; j<S; j++)
		    psum -= L[j][col]*x[j];		
		I[col][col] = invdiag[col]*psum;
	    }
	}
	
    }

    template <int Size=-1>
    class Cholesky {
    public:

	template<class Accessor>
	Cholesky(const FixedMatrix<Size,Size,Accessor>& m){
	    compute(m);
	}
    
	template<class Accessor>
	void compute(const FixedMatrix<Size,Size,Accessor>& m){
	    rank = util::cholesky_compute(m,L,invdiag);
	}
	int get_rank() const { return rank; }

	double get_determinant() const {
	    double det = L[0][0];
	    for (int i=1; i<Size; i++)
		det *= L[i][i];
	    return det*det;
	}

	const Matrix<Size>& get_L() const {
	    return L;
	}

	template <class Accessor> inline 
	Vector<Size> backsub(const FixedVector<Size,Accessor>& v) const
	{
	    Vector<Size> x;
	    util::cholesky_backsub(L, invdiag, v, x);
	    return x;
	}
      
	template <class Accessor, int Cols> inline Matrix<Size,Cols> inverse_times(const FixedMatrix<Size,Cols,Accessor>& m)
	{
	    Matrix<Cols, Size> result;
	    for (int i=0; i<Cols; i++)
		result[i] = backsub(m.T()[i]);
	    return result.T();
	}

	template <class A> void get_inverse(FixedMatrix<Size,Size,A>& M) const {
	    util::cholesky_inverse(L, invdiag, M);
	}

	Matrix<Size,Size> get_inverse() const {
	    Matrix<Size,Size> M;
	    get_inverse(M);
	    return M;
	}

    private:
	Matrix<Size,Size,RowMajor> L;
	Vector<Size> invdiag;
	int rank;
    };
  
    template <>
    class Cholesky<-1> {
    public:

	template<class Accessor>
	Cholesky(const DynamicMatrix<Accessor>& m) : L(m.num_rows(), m.num_cols()), invdiag(m.num_rows()) {
	    assert(m.num_rows() == m.num_cols());
	    compute(m);
	}

	template<class Accessor>
	void compute(const DynamicMatrix<Accessor>& m){
	    int Size = m.num_rows();
	    rank = Size;
	    for(int i=0;i<Size;i++) {
		double a=m[i][i];
		for(int k=0;k<i;k++) 
		    a-=L[i][k]*L[i][k];
		if (a < std::numeric_limits<double>::epsilon()) {
		    --rank;
		    L[i][i] = 0;
		} else
		    L[i][i]=sqrt(a);
		invdiag[i] = 1.0/L[i][i];
		for(int j=i+1;j<Size;j++) {
		    a=m[i][j];
		    for(int k=0;k<i;k++) 
			a-=L[i][k]*L[j][k];
		    L[j][i]=a*invdiag[i];
		}
	    }
	}
	int get_rank() const { return rank; }
	template <class Accessor> inline 
	Vector<> backsub(const DynamicVector<Accessor>& v) const
	{
	    assert(v.size() == L.num_rows());
	    return backsub(v.get_data_ptr());
	}

	template <int N, class Accessor> inline
	Vector<> backsub(const FixedVector<N, Accessor>& v) const
	{
	    assert(N == L.num_rows());
	    return backsub(v.get_data_ptr());
	}
    
	const Matrix<>& get_L() const {
	    return L;
	}

	double get_determinant() const {
	    double det = L[0][0];
	    for (int i=1; i<L.num_rows(); i++)
		det *= L[i][i];
	    return det*det;
	}

	Matrix<> get_inverse() const {
	    int Size = L.num_rows();
	    Matrix<> M(Size, Size);
	    // compute inverse(L)
	    for (int r=0; r<Size; r++) {
		M[r][r] = invdiag[r];
		for (int i=r+1; i<Size; i++) {
		    double b = -L[i][r]*M[r][r];
		    for (int j=r+1; j<i;j++)
			b -= L[i][j]*M[r][j];
		    M[r][i] = b*invdiag[i];
		}      
	    }
	    // multiply (inverse(L))^T * inverse(L)
	    for (int i=0;i<Size; i++) {
		double s = 0;
		for (int k=i; k<Size; k++)
		    s += M[i][k]*M[i][k];
		M[i][i] = s;
		for (int j=i+1; j<Size; j++) {
		    double s = 0;
		    for (int k=j; k<Size; k++)
			s += M[j][k]*M[i][k];
		    M[i][j] = M[j][i] = s;
		}
	    }
	    return M;
	}

    private:
	Vector<> backsub(const double* v) const
	{
	    int Size = L.num_rows();
	    Vector<> t(Size);
	    // forward substitution
	    for (int i=0; i<Size; i++) {
		double b = v[i];
		for (int j=0; j<i;j++)
		    b -= L[i][j]*t[j];
		t[i] = b*invdiag[i];
	    }
	    // back substitution
	    Vector<> x(Size);
	    for (int i=Size-1; i >=0; i--) {
		double b = t[i];
		for (int j=i+1; j<Size;j++)
		    b -= L[j][i]*x[j];
		x[i] = b*invdiag[i];
	    }
	    return x;
	}
	Matrix<> L;
	Vector<> invdiag;
	int rank;
    };

#ifndef TOON_NO_NAMESPACE
}
#endif 





#endif
