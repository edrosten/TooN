//-*- c++ -*-

// Copyright (C) 2012, Edward Rosten (ed@edwardrosten.com)

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


#ifndef TOON_INC_QR_H
#define TOON_INC_QR_H
#include <TooN/TooN.h>
#include <cmath>

namespace TooN
{
/**
Performs %QR decomposition.

@warning this will only work if the number of columns is greater than 
the number of rows!

The QR decomposition operates on a matrix A. 
In general:
\f[
A = QR
\f]

@ingroup gDecomps
*/
template<int Rows=Dynamic, int Cols=Rows, typename Precision=double> class QR
{
	
	private:
		static const int square_Size = (Rows>=0 && Cols>=0)?(Rows<Cols?Rows:Cols):Dynamic;

	public:	
		/// Construct the %QR decomposition of a matrix. This initialises the class, and
		/// performs the decomposition immediately.
		/// @param m The matrix to decompose
		template<int R, int C, class P, class B> 
		QR(const Matrix<R,C,P,B>& m_)
		:m(m_), Q(Identity(square_size()))
		{
			//pivot is set to all zeros, which means all columns are free columns
			//and can take part in column pivoting.

			compute();
		}
		
		///Return R
		const Matrix<Rows, Cols, Precision>& get_R()
		{
			return m;
		}
		
		///Return Q
		const Matrix<square_Size, square_Size, Precision, ColMajor>& get_Q()
		{
			return Q;
		}	

	private:

		template<class B1, class B2>		
		void pre_multiply_by_householder(Matrix<Dynamic, Dynamic, Precision, B1>& A, const Vector<Dynamic, Precision, B2>& v, Precision s)
		{
			//The Householder matrix is P = 1 - v*v/s
			
			//PA = (1 - v*v^T/s) A = A - v*v^T A 
			//                                               
			//          [ v1 ]  [ v1   v2  v2  v4 ]  [    |     |     ]
			//   = A -  [ v2 ]                       [ a1 |  a2 | ... ]   
			//          [ v3 ]                       [    |     |     ]
			//          [ v4 ]                       [    |     |     ]

			//          [ v1 ]  [ va1   va2  va3 va4 ]
			//   = A -  [ v2 ]                      
			//          [ v3 ]                     
			//          [ v4 ]                    

			//          [ v1 (v a1) ]
			//   = A -  [ v2 (v a2) ]
			//          [ v3 (v a3) ]
			//          [ v4 (v a4) ]

			for(int col=0; col < A.num_cols(); col++)
			{
				Precision va = v * A.T()[col] / s;

				for(int row=0; row < A.num_rows(); row++)
					A[row][col] -= v[row] * va;
			}
		}

		void compute()
		{
			
			//QR decomposition makes use of Householder reflections.
			// A = QR, 

			//Q1 A = Q1 QR
			// Q2 Q1 A = Q2 Q1 R
			// ... 
			// Q^-1 A  = R
			//
			// Pick a sequence of Qn which makes R upper triangular.
			//
			// The sequence Qn ... Q1 is equal to Q^-1
			//
			// Qi has the form of a Householder reflection

			for(int n=0; n < square_size()-1; n++)	
			{
				using std::sqrt;

				int sz = square_size() - n;
				int nc = m.num_cols() - n;

				//Compute the reflection on a submatrix
				//such that it never breaks the triangular 
				//properties of the matrix being created.

				Matrix<Dynamic, Dynamic, Precision, typename Matrix<Rows,Cols,Precision>::SliceBase> s = m.slice(n, n, sz, nc);

				//Now perform a householder reduction on the first column

				//Define x to be the first column
				//auto x = s.T()[0];


				//The reflection vector is u = x - |x| * [ 1 0 ... 0] * sgn(x_0)
				//
				//let a = |x| * sgn(x_0])

				Precision  nx2 = norm_sq(s.T()[0]);

				Precision a = sqrt(nx2) * (s.T()[0][0] < 0?-1:1);

				//We now want the vector u = x - ae
				//
				// Multipling (I - 2 hat(u) * hat(u)^T) * s will zero out the first column of s except
				// for the first element. The matrix P = is orthogonal, a bit like Qn 
				//
				//Since x is no longer needed, we can just modify the first element
				s.T()[0][0] -= a;
				//auto& u = x;

				//We want H = norm_sq(u)/2 =  a (a - x_0) = a * -u_0
				Precision H = -a * s.T()[0][0];


				//Note that we're working on a reduced sized matrix here.
				//
				//The actual householder matrix, P,  is:
				//                                                                  
				//  [ 1  |   0         ]                                               
				//  [___1|_____________]                                                
				//  [    |[ ^   ^T  ]  ]                                                
				//  [  0 |[ u * u   ]  ]                                                        
				//  [    |[         ]  ]                                                

				// We want m <- P * m
				// and Q <- P * Q
				//
				//
				// Since m is going towards upper triangular and so has lots of zeros
				// we can compute it by performing the multiplication in just the
				// lower block:
				//

				//  [ 1  |   0         ]  [ m1 | m2     ]     [ m1 |  m2  ]
				//  [___1|_____________]  [____|_______ ]     [____|______]               
				//  [    |[  ^  ^T  ]  ]  [    |        ]   = [    |      ]               
				//  [  0 |[I-u* u   ]  ]  [ 0  |  m3    ]     [  0 | .... ]                       
				//  [    |[         ]  ]  [    |        ]     [    |      ]                
				
				// Q does not have the column of zeros, so the multiplication has to be performed on
				// the whole width

				//This can be done in place except that the first column of s holds u
				pre_multiply_by_householder(s.slice(0, 1, sz, nc-1).ref(), s.T()[0], H);

				//Also update the Q matrix
				pre_multiply_by_householder(Q.slice(n,0,sz,square_size()).ref(), s.T()[0], H);

				//Now do the first column multiplication...
				//Which makes it upper triangular.

				s[0][0] = a;
				s.T()[0].slice(1, sz-1) = Zeros;

			}
				
			//Note above that we've build up Q^-1, so we need to transpose Q now 
			//to invert it

			using std::swap;
			for(int r=0; r < Q.num_rows(); r++)
				for(int c=r+1; c < Q.num_cols(); c++)
					swap(Q[r][c], Q[c][r]);

		}

		Matrix<Rows, Cols, Precision> m;
		Matrix<square_Size, square_Size, Precision, ColMajor> Q;
		

		int square_size()
		{
			return std::min(m.num_rows(), m.num_cols());	
		}
};






}

#endif
