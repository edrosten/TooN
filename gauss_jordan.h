#ifndef TOON_INC_GAUSS_JORDAN_H
#define TOON_INC_GAUSS_JORDAN_H

#include <utility>
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
template<int R, int C, class X> void gauss_jordan(FixedMatrix<R, C, X>& m)
{
	using std::swap;

	//Loop over columns to reduce.
	for(int col=0; col < R; col++)
	{
		//Reduce the current column to a single element


		//Search down the current column in the lower triangle for the largest
		//absolute element (pivot).  Then swap the pivot row, so that the pivot
		//element is on the diagonal. The benchmarks show that it is actually
		//faster to swap whole rows than it is to access the rows via indirection 
		//and swap the indirection element. This holds for both pointer indirection
		//and using a permutation vector over rows.
		{
			int pivotpos = col;
			double pivotval = abs(m[pivotpos][col]);
			for(int p=col+1; p <R; p++)
				if(abs(m[p][col]) > pivotval)
				{
					pivotpos = p;
					pivotval = abs(m[pivotpos][col]);
				}

			swap(m[col], m[pivotpos]);
		}

		//Reduce the current column in every row to zero, excluding elements on
		//the leading diagonal.
		for(int row = 0; row < R; row++)
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
				for(int c=col+1; c < C; c++)
					m[row][c] = m[row][c] - m[col][c] * multiple;
			}
		}
	}
	
	//Final pass to make diagonal elements one. Performing this in a final
	//pass allows us to avoid any significant computations on the left-hand
	//square matrix, since it is diagonal, and ends up as the identity.
	for(int row=0;row < R; row++)
	{
		double mul = 1/m[row][row];

		m[row][row] = 1;

		for(int col=R; col < C; col++)
			m[row][col] *= mul;
	}
}

}
#endif
