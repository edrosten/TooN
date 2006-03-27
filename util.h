#ifndef UTIL_H
#define UTIL_H

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

    namespace util {
	template <int B, int E> struct Dot {
	    template <class V1, class V2> static inline double eval(const V1& v1, V2& v2) { return v1[B]* v2[B] + Dot<B+1,E>::eval(v1,v2); }
	};
	
	template <int N> struct Dot<N,N> {
	    template <class V1, class V2> static inline double eval(const V1& v1, V2& v2) { return v1[N]*v2[N];  }
	};	
    }



#ifndef TOON_NO_NAMESPACE
}
#endif 





#endif 
