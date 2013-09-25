// -*- c++ -*-

//    Copyright (C) 2009 Ed Rosten (er258@cam.ac.uk)

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


#ifndef TOON_INC_GAUSS_JORDAN_H
#define TOON_INC_GAUSS_JORDAN_H

#include <utility>
#include <cmath>
#include <TooN/TooN.h>

namespace TooN
{
/// Perform Gauss-Jordan reduction on m
///
/// If m is of the form \f$[A | I ]\f$, then after reduction, m
/// will be \f$[ I | A^{-1}]\f$. There is no restriction on the input, 
/// in that the matrix augmenting A does not need to be I, or square.
/// The reduction is performed using elementary row operations and 
/// partial pivoting.
///
/// @param m The matrix to be reduced.
/// @ingroup gDecomps
template<int R, int C, class Precision, class Base> void gauss_jordan(Matrix<R, C, Precision, Base>& m)
{
	using std::swap;

	//Loop over columns to reduce.
	for(int col=0; col < m.num_rows(); col++)
	{
		//Reduce the current column to a single element


		//Search down the current column in the lower triangle for the largest
		//absolute element (pivot).  Then swap the pivot row, so that the pivot
		//element is on the diagonal. The benchmarks show that it is actually
		//faster to swap whole rows than it is to access the rows via indirection 
		//and swap the indirection element. This holds for both pointer indirection
		//and using a permutation vector over rows.
		{
		  using std::abs;
			int pivotpos = col;
			double pivotval = abs(m[pivotpos][col]);
			for(int p=col+1; p <m.num_rows(); p++)
			  if(abs(m[p][col]) > pivotval)
				{
					pivotpos = p;
					pivotval = abs(m[pivotpos][col]);
				}
			
			if(col != pivotpos)
				swap(m[col].ref(), m[pivotpos].ref());
		}

		//Reduce the current column in every row to zero, excluding elements on
		//the leading diagonal.
		for(int row = 0; row < m.num_rows(); row++)
		{
			if(row != col)
			{
				double multiple = m[row][col] / m[col][col];
		
				//We could eliminate some of the computations in the augmented
				//matrix, if the augmented half is the identity. In general, it
				//is not. 

				//Subtract the pivot row from all other rows, to make 
				//column col zero.
				m[row][col] = 0;
				for(int c=col+1; c < m.num_cols(); c++)
					m[row][c] = m[row][c] - m[col][c] * multiple;
			}
		}
	}
	
	//Final pass to make diagonal elements one. Performing this in a final
	//pass allows us to avoid any significant computations on the left-hand
	//square matrix, since it is diagonal, and ends up as the identity.
	for(int row=0;row < m.num_rows(); row++)
	{
		double mul = 1/m[row][row];

		m[row][row] = 1;

		for(int col=m.num_rows(); col < m.num_cols(); col++)
			m[row][col] *= mul;
	}
}

}
#endif
