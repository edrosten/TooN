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


#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <TooN/TooN.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 



/// Cholesky decomposition of a symmetric matrix.
/// Only the lower half of the matrix is considered
/// This uses the non-sqrt version of the decomposition
/// giving symmetric M = L*D*L.T() where the diagonal of L contains ones
template <int Size, class Precision=double>
class Cholesky {
public:
	Cholesky(){}

	template<class P2, class B2>
	Cholesky(const Matrix<Size, Size, P2, B2>& m)
		: my_cholesky(m) {
		compute(m);
	}
	
	// for Size=Dynamic
	Cholesky(int size) : my_cholesky(size,size) {}


	template<class P2, class B2> void compute(const Matrix<Size, Size, P2, B2>& m){
		SizeMismatch<Size,Size>::test(m.num_rows(), m.num_cols());
		SizeMismatch<Size,Size>::test(m.num_rows(), my_cholesky.num_rows());
		my_cholesky=m;
		int size=my_cholesky.num_rows();
		for(int col=0; col<size; col++){
			Precision inv_diag;
			for(int row=col; row < size; row++){
				// correct for the parts of cholesky already computed
				Precision val = my_cholesky(row,col);
				for(int col2=0; col2<col; col2++){
					// val-=my_cholesky(col,col2)*my_cholesky(row,col2)*my_cholesky(col2,col2);
					val-=my_cholesky(col2,col)*my_cholesky(row,col2);
				}
				if(row==col){
					// this is the diagonal element so don't divide
					my_cholesky(row,col)=val;
					inv_diag=1/val;
				} else {
					// cache the value without division in the upper half
					my_cholesky(col,row)=val;
					// divide my the diagonal element for all others
					my_cholesky(row,col)=val*inv_diag;
				}
			}
		}
	}

	template<int Size2, class P2, class B2>
	Vector<Size, Precision> backsub (const Vector<Size2, P2, B2>& v) {
		int size=my_cholesky.num_rows();
		SizeMismatch<Size,Size2>::test(size, v.size());

		// first backsub through L
		Vector<Size, Precision> y(size);
		for(int i=0; i<size; i++){
			Precision val = v[i];
			for(int j=0; j<i; j++){
				val -= my_cholesky(i,j)*y[j];
			}
			y[i]=val;
		}
		
		// backsub through diagonal
		for(int i=0; i<size; i++){
			y[i]/=my_cholesky(i,i);
		}

		// backsub through L.T()
		Vector<Size,Precision> result(size);
		for(int i=size-1; i>=0; i--){
			Precision val = y[i];
			for(int j=i+1; j<size; j++){
				val -= my_cholesky(j,i)*result[j];
			}
			result[i]=val;
		}

		return result;
	}


	template<int Size2, int C2, class P2, class B2>
	Matrix<Size, C2, Precision> backsub (const Matrix<Size2, C2, P2, B2>& m) {
		int size=my_cholesky.num_rows();
		SizeMismatch<Size,Size2>::test(size, m.num_rows());

		// first backsub through L
		Matrix<Size, C2, Precision> y(size, m.num_cols());
		for(int i=0; i<size; i++){
			Vector<C2, Precision> val = m[i];
			for(int j=0; j<i; j++){
				val -= my_cholesky(i,j)*y[j];
			}
			y[i]=val;
		}
		
		// backsub through diagonal
		for(int i=0; i<size; i++){
			y[i]*=(1/my_cholesky(i,i));
		}

		// backsub through L.T()
		Matrix<Size,C2,Precision> result(size, m.num_cols());
		for(int i=size-1; i>=0; i--){
			Vector<C2,Precision> val = y[i];
			for(int j=i+1; j<size; j++){
				val -= my_cholesky(j,i)*result[j];
			}
			result[i]=val;
		}
		return result;
	}


	// easy way to get inverse - could be made more efficient
	Matrix<Size,Size,Precision> get_inverse(){
		Matrix<Size,Size,Precision>I(Identity(my_cholesky.num_rows()));
		return backsub(I);
	}

	Precision determinant(){
		Precision answer=my_cholesky(0,0);
		for(int i=1; i<my_cholesky.num_rows(); i++){
			answer*=my_cholesky(i,i);
		}
		return answer;
	}

private:
	Matrix<Size,Size,Precision> my_cholesky;
};


#ifndef TOON_NO_NAMESPACE
}
#endif 

#endif
