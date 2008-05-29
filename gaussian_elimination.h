#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <utility>
#include <TooN/TooN.h>

namespace TooN {
    template <class R, int N, class B> inline
    R gaussian_elimination(Matrix<N> A, B b)
    {
	using std::swap;

	for (int i=0; i<N; ++i) {
	    int argmax = i;
	    double maxval = abs(A[i][i]);
	
	    for (int ii=i+1; ii<N; ++ii) {
		double v =  abs(A[ii][i]);
		if (v > maxval) {
		    maxval = v;
		    argmax = ii;
		}
	    }
	    double pivot = A[argmax][i];
	    //assert(abs(pivot) > 1e-16);
	    double inv_pivot = 1.0/pivot;
	    if (argmax != i) {
		for (int j=i; j<N; ++j)
		    swap(A[i][j], A[argmax][j]);
		swap(b[i], b[argmax]);
	    }
	    //A[i][i] = 1;
	    for (int j=i+1; j<N; ++j)
		A[i][j] *= inv_pivot;
	    b[i] *= inv_pivot;

	    for (int u=i+1; u<N; ++u) {
		double factor = A[u][i];
		//A[u][i] = 0;
		for (int j=i+1; j<N; ++j)
		    A[u][j] -= factor * A[i][j];
		b[u] -= factor * b[i];
	    }
	}
    
	R x;
	for (int i=N-1; i>=0; --i) {
	    x[i] = b[i];
	    for (int j=i+1; j<N; ++j)
		x[i] -= A[i][j] * x[j];
	}
	return x;
    }

    template <int N, class A1, class A2>
    inline Vector<N> gaussian_elimination(const FixedMatrix<N,N,A1>& A, const FixedVector<N,A2>& b) {
	return gaussian_elimination<Vector<N>, N, Vector<N> >(Matrix<N>(A), Vector<N>(b));
    }

    template <int N, class A1, class A2, int M>
    inline Matrix<N> gaussian_elimination(const FixedMatrix<N,N,A1>& A, const FixedMatrix<N,M,A2>& b) {
	return gaussian_elimination<Matrix<N,M>,N, Matrix<N,M> >(Matrix<N>(A), Matrix<N,M>(b));
    }

}
#endif
