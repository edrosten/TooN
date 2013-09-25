// -*- c++ -*-

// Copyright (C) 2008,2009 Ethan Eade, Tom Drummond (twd20@cam.ac.uk)
// and Ed Rosten (er258@cam.ac.uk)

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


#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <utility>
#include <cmath>
#include <TooN/TooN.h>

namespace TooN {
	///@ingroup gEquations
	///Return the solution for \f$Ax = b\f$, given \f$A\f$ and \f$b\f$
	///@param A \f$A\f$
	///@param b \f$b\f$
    template<int N, typename Precision>
	inline Vector<N, Precision> gaussian_elimination (Matrix<N,N,Precision> A, Vector<N, Precision> b) {
		using std::swap;
		using std::abs;

		int size=b.size();

		for (int i=0; i<size; ++i) {
			int argmax = i;
			Precision maxval = abs(A[i][i]);
			
			for (int ii=i+1; ii<size; ++ii) {
				double v =  abs(A[ii][i]);
				if (v > maxval) {
					maxval = v;
					argmax = ii;
				}
			}
			Precision pivot = A[argmax][i];
			//assert(abs(pivot) > 1e-16);
			Precision inv_pivot = static_cast<Precision>(1)/pivot;
			if (argmax != i) {
				for (int j=i; j<size; ++j)
					swap(A[i][j], A[argmax][j]);
				swap(b[i], b[argmax]);
			}
			//A[i][i] = 1;
			for (int j=i+1; j<size; ++j)
				A[i][j] *= inv_pivot;
			b[i] *= inv_pivot;
			
			for (int u=i+1; u<size; ++u) {
				double factor = A[u][i];
				//A[u][i] = 0;
				for (int j=i+1; j<size; ++j)
					A[u][j] -= factor * A[i][j];
				b[u] -= factor * b[i];
			}
		}
		
		Vector<N,Precision> x(size);
		for (int i=size-1; i>=0; --i) {
			x[i] = b[i];
			for (int j=i+1; j<size; ++j)
				x[i] -= A[i][j] * x[j];
		}
		return x;
    }
	
	namespace Internal
	{
		template<int i, int j, int k> struct Size3
		{
			static const int s = !IsStatic<i>::is?i: (!IsStatic<j>::is?j:k);
		};

	};

	///@ingroup gEquations
	///Return the solution for \f$Ax = b\f$, given \f$A\f$ and \f$b\f$
	///@param A \f$A\f$
	///@param b \f$b\f$
    template<int R1, int C1, int R2, int C2, typename Precision>
	inline Matrix<Internal::Size3<R1, C1, R2>::s, C2, Precision> gaussian_elimination (Matrix<R1,C1,Precision> A, Matrix<R2, C2, Precision> b) {
		using std::swap;
		using std::abs;
		SizeMismatch<R1, C1>::test(A.num_rows(), A.num_cols());
		SizeMismatch<R1, R2>::test(A.num_rows(), b.num_rows());

		int size=A.num_rows();

		for (int i=0; i<size; ++i) {
			int argmax = i;
			Precision maxval = abs(A[i][i]);
			
			for (int ii=i+1; ii<size; ++ii) {
				Precision v =  abs(A[ii][i]);
				if (v > maxval) {
					maxval = v;
					argmax = ii;
				}
			}
			Precision pivot = A[argmax][i];
			//assert(abs(pivot) > 1e-16);
			Precision inv_pivot = static_cast<Precision>(1)/pivot;
			if (argmax != i) {
				for (int j=i; j<size; ++j)
					swap(A[i][j], A[argmax][j]);

				for(int j=0; j < b.num_cols(); j++)
					swap(b[i][j], b[argmax][j]);
			}
			//A[i][i] = 1;
			for (int j=i+1; j<size; ++j)
				A[i][j] *= inv_pivot;
			b[i] *= inv_pivot;
			
			for (int u=i+1; u<size; ++u) {
				Precision factor = A[u][i];
				//A[u][i] = 0;
				for (int j=i+1; j<size; ++j)
					A[u][j] -= factor * A[i][j];
				b[u] -= factor * b[i];
			}
		}
		
		Matrix<Internal::Size3<R1, C1, R2>::s,C2,Precision> x(b.num_rows(), b.num_cols());
		for (int i=size-1; i>=0; --i) {
			for(int k=0; k <b.num_cols(); k++)
			{
				x[i][k] = b[i][k];
				for (int j=i+1; j<size; ++j)
					x[i][k] -= A[i][j] * x[j][k];
			}
		}
		return x;
    }
}
#endif
