#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <utility>
#include <TooN/TooN.h>

namespace TooN {
    template<int N, typename Precision>
	inline Vector<N, Precision>::type>  gaussian_elimination (Matrix<N,N,Precision> A, Vector<N, Precision> b) {
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
		
		Vector<N,Precision>::type> x(size);
		for (int i=size-1; i>=0; --i) {
			x[i] = b[i];
			for (int j=i+1; j<size; ++j)
				x[i] -= A[i][j] * x[j];
		}
		return x;
    }
}
#endif
