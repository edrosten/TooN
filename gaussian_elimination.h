#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <utility>
#include <TooN/TooN.h>

namespace TooN {
    template<int N, typename Precision>
	inline Vector<N, Precision> gaussian_elimination (Matrix<N,N,Precision> A, Vector<N, Precision> b) {
		using std::swap;

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
			static const int s=(i!= -1)?i:(j!=-1?j:k);
		};

	};

    template<int R1, int C1, int R2, int C2, typename Precision>
	inline Matrix<Internal::Size3<R1, C1, R2>::s, C2, Precision> gaussian_elimination (Matrix<R1,C1,Precision> A, Matrix<R2, C2, Precision> b) {
		using std::swap;
		SizeMismatch<R1, C1>::test(A.num_rows(), A.num_cols());
		SizeMismatch<R1, R2>::test(A.num_rows(), b.num_rows());

		int size=A.num_rows();

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

				for(int j=0; j < b.num_cols(); j++)
					swap(b[i][j], b[argmax][j]);
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
