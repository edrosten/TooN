#ifndef TOON_INCLUDE_QR_LAPACK_H
#define TOON_INCLUDE_QR_LAPACK_H


#include <TooN/TooN.h>
#include <TooN/lapack.h>
#include <utility>

namespace TooN{

/**
Performs %QR decomposition.

***WARNING*** this will only work if the number of columns is greater than 
the number of rows!

@ingroup gDecomps
*/
template<int Rows=Dynamic, int Cols=Rows, class Precision=double>
class QR_Lapack{

	private:
		static const int square_Size = (Rows>=0 && Cols>=0)?(Rows<Cols?Rows:Cols):Dynamic;

	public:	
		/// Construct the %QR decomposition of a matrix. This initialises the class, and
		/// performs the decomposition immediately.
		template<int R, int C, class P, class B> 
		QR_Lapack(const Matrix<R,C,P,B>& m)
		:copy(m),tau(square_size()), Q(square_size(), square_size())
		{
			compute();
		}
		
		///Return L
		const Matrix<Rows, Cols, Precision, ColMajor>& get_R()
		{
			return copy;
		}
		
		///Return Q
		const Matrix<square_Size, square_Size, Precision, ColMajor>& get_Q()
		{
			return Q;
		}	


	private:

		void compute()
		{	
			int M = copy.num_rows();
			int N = copy.num_cols();
			
			int LWORK=-1;
			int INFO;
			int lda = M;

			Precision size;
			
			//Compute the working space
			geqrf_(&M, &N, copy.get_data_ptr(), &lda, tau.get_data_ptr(), &size, &LWORK, &INFO);

			LWORK = (int) size;

			Precision* work = new Precision[LWORK];

			geqrf_(&M, &N, copy.get_data_ptr(), &lda, tau.get_data_ptr(), work, &LWORK, &INFO);


			if(INFO < 0)
				std::cerr << "error in QR, INFO was " << INFO << std::endl;

			//The upper "triangle+" of copy is R
			//The lower right and tau contain enough information to reconstruct Q
			
			//LAPACK provides a handy function to do the reconstruction
			Q = copy.template slice<0,0,square_Size, square_Size>(0,0,square_size(), square_size());
			
			int K = square_size();
			M=K;
			N=K;
			lda = K;
			orgqr_(&M, &N, &K, Q.get_data_ptr(), &lda, tau.get_data_ptr(), work, &LWORK, &INFO);

			if(INFO < 0)
				std::cerr << "error in QR, INFO was " << INFO << std::endl;

			delete [] work;
			
			//Now zero out the lower triangle
			for(int r=1; r < square_size(); r++)
				for(int c=0; c<r; c++)
					copy[r][c] = 0;

		}

		Matrix<Rows, Cols, Precision, ColMajor> copy;
		Matrix<square_Size, square_Size, Precision, ColMajor> Q;
		Vector<square_Size, Precision> tau;


		int square_size()
		{
			return std::min(copy.num_rows(), copy.num_cols());	
		}


};

}


#endif
