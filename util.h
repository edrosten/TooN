#ifndef UTIL_H
#define UTIL_H

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 
    namespace util {

	template <bool Cond> struct Assert;
	template <> struct Assert<true> {};

	template <int B, int E, bool Valid=(B<=E)> struct Dot { 
	    template <class V1, class V2> static inline double eval(const V1&, const V2&) { return 0; }
	};

	template <int B, int E> struct Dot<B,E,true> {
	    template <class V1, class V2> static inline double eval(const V1& v1, const V2& v2) { return v1[B]* v2[B] + Dot<B+1,E>::eval(v1,v2); }
	};
	
	template <int N> struct Dot<N,N,true> {
	    template <class V1, class V2> static inline double eval(const V1& v1, const V2& v2) { return v1[N]*v2[N];  }
	};	


	template <int B, int E> struct AddV {
	    template <class V1, class V2> static inline void eval(V1& v1, const V2& v2) { v1[B]+=v2[B]; AddV<B+1,E>::eval(v1,v2); }
	};
	
	template <int N> struct AddV<N,N> {
	    template <class V1, class V2> static inline void eval(V1& v1, const V2& v2) { v1[N] += v2[N];  }
	};	
	
#if 0
	template <int A, int N, int B, int Row, int Col=0> struct MatrixProductRow {
	    template <class M1, class M2, class V> static inline void eval(const M1& a, const M2& b, V& v) { 		
		v[Col] = Dot<0,N-1>::eval(a[Row], b[Col]);
		MatrixProductRow<A,N,B,Row,Col+1>::eval(a,b,v);
	    }
	};

	template <int A, int N, int B, int Row> struct MatrixProductRow<A,N,B,Row,B> {
	    template <class M1, class M2, class V> static inline void eval(const M1& a, const M2& b, V& v) {
		v[B] = Dot<0,N-1>::eval(a[Row], b[B]);
	    }
	};
	
	template <int A, int N, int B, int Row=0> struct MatrixProduct {
	    template <class M1, class M2, class M3> static inline void eval(const M1& a, const M2& b, M3& m) {
		MatrixProductRow<A,N,B,Row>::eval(a,b,m[Row]);
		MatrixProduct<A,N,B,Row+1>::eval(a,b,m);
	    }
	};
	
	template <int A, int N, int B> struct MatrixProduct<A,N,B,A> {
	    template <class M1, class M2, class M3> static inline void eval(const M1& a, const M2& b, M3& m) {
		MatrixProductRow<A,N,B,A>::eval(a,b,m[A]);
	    }
	};

	template <int A, int N, int B, class M1, class M2, class M3> inline void matrix_multiply(const M1& a, const M2& b, M3& m) {
	    MatrixProduct<A-1,N,B-1>::eval(a,b.T(),m);
	}
#endif

	template <int B, int Col=0> struct MatrixProductRow {
	    template <class F, class M1, class M2, class V> static inline void eval(const M1& a, const M2& b, V& v, int row) { 		
		F::eval(v[Col], a[row] * b[Col]);
		MatrixProductRow<B,Col+1>::template eval<F>(a,b,v, row);
	    }
	};

	template <int B> struct MatrixProductRow<B,B> {
	    template <class F, class M1, class M2, class V> static inline void eval(const M1& a, const M2& b, V& v, int row) {
		F::eval(v[B], a[row] * b[B]);
	    }
	};
	struct Assign { template <class LHS, class RHS> static inline void eval(LHS& lhs, const RHS& rhs) { lhs = rhs; } };
	struct PlusEquals { template <class LHS, class RHS> static inline void eval(LHS& lhs, const RHS& rhs) { lhs += rhs; } };
	struct MinusEquals { template <class LHS, class RHS> static inline void eval(LHS& lhs, const RHS& rhs) { lhs -= rhs; } };

	template <class F, int A, int N, int B, class M1, class M2, class M3> inline void matrix_multiply(const M1& a, const M2& b, M3& m) {
	    for (int i=0; i<m.num_rows(); i++)
		MatrixProductRow<B-1>::template eval<F>(a,b.T(),m[i],i);
	}
	    	
	template <int A, int N, int B, class M1, class M2, class M3> inline void matrix_multiply(const M1& a, const M2& b, M3& m) {
	    matrix_multiply<Assign,A,N,B,M1,M2,M3>(a,b,m);
	}

    }



#ifndef TOON_NO_NAMESPACE
}
#endif 





#endif 
