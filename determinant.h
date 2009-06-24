#ifndef TOON_INCLUDE_DETERMINANT_H
#define TOON_INCLUDE_DETERMINANT_H
#ifdef TOON_DETERMINANT_LAPACK
	#include <TooN/LU.h>
#endif

namespace TooN
{
	namespace Internal
	{
		template<int R, int C> struct Square
		{
		};

		template<int R> struct Square<R, R>
		{
			static const int Size = R;
		};

		template<int R> struct Square<R, Dynamic>
		{
			static const int Size = R;
		};
		template<int C> struct Square<Dynamic, C>
		{
			static const int Size = C;
		};
		template<> struct Square<Dynamic, Dynamic>
		{
			static const int Size = Dynamic;
		};
	};


	/** Compute the determinant using Gaussian elimination.
		@param A The matrix to find the determinant of.
		@returns determinant.
		@ingroup gLinAlg
	*/
	template<int R, int C, typename Precision, typename Base>
	Precision determinant_gaussian_elimination(const Matrix<R, C, Precision, Base>& A_)
	{
		Matrix<Internal::Square<R,C>::Size, Internal::Square<R,C>::Size,Precision> A = A_;
		TooN::SizeMismatch<R, C>::test(A.num_rows(), A.num_cols());
		using std::swap;

		int size=A.num_rows();
		
		//If row operations of the form row_a += alpha * row_b
		//then the determinant is unaffected. However, if a row
		//is scaled, then the determinant is scaled by the same 
		//amount. The total scaling is held in divisor.
		Precision determinant=1;
		Precision divisor=1;

		for (int i=0; i<size; ++i) {

			//Find the pivot
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
			
			//Swap the current row with the pivot row if necessary.
			if (argmax != i) {
				for (int j=i; j<size; ++j)
					swap(A[i][j], A[argmax][j]);
			}

			determinant *= A[i][i];

			if(determinant == 0)
				return 0;
			
			for (int u=i+1; u<size; ++u) {
				//Multiply out the usual 1/pivot term
				//to avoid division.
				double factor = A[u][i];
				//A[u][i] = 0;
				divisor*=pivot;
				for (int j=i+1; j<size; ++j)
					A[u][j] = A[u][j]*pivot - factor * A[i][j];
			}
		}

		if(size %2 == 0) 
			return -determinant/divisor;
		else
			return determinant/divisor;
	}
	
	/** Compute the determinant using TooN::LU.
		@param A The matrix to find the determinant of.
		@returns determinant.
		@ingroup gLinAlg
	*/
	template<int R, int C, class P, class B>
	P determinant_LU(const Matrix<R, C, P, B>& A)
	{
		TooN::SizeMismatch<R, C>::test(A.num_rows(), A.num_cols());
		LU<Internal::Square<R,C>::Size, P> lu(A);
		return lu.determinant();
	}

	/** 
		Compute the determinant of a matrix using an appropriate method. The
		obvious method is used for 2x2, otherwise
		determinant_gaussian_elimination() or determinant_LU() is used depending
		on the value of \c TOON_DETERMINANT_LAPACK.  See also \ref sConfigLapack.
		@param A The matrix to find the determinant of.
		@returns determinant.
		@ingroup gLinAlg
	*/
	template<int R, int C, class P, class B>
	P determinant(const Matrix<R, C, P, B>& A)
	{
		TooN::SizeMismatch<R, C>::test(A.num_rows(), A.num_cols());
		if(A.num_rows() == 2)
			return A[0][0]*A[1][1] - A[1][0]*A[0][1];
		#if defined TOON_DETERMINANT_LAPACK && TOON_DETERMINANT_LAPACK != -1
			else if(A.num_rows() >= TOON_DETERMINANT_LAPACK)
				return determinant_LU(A);
		#endif
		else
			return determinant_gaussian_elimination(A);
	}
}

#endif
