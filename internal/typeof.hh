#ifdef TOON_TYPEOF_DECLTYPE
	#define TOON_TYPEOF(X) decltype((X))
#elif defined TOON_TYPEOF_TYPEOF
	#define TOON_TYPEOF(X) typeof((X))
#elif defined TOON_TYPEOF___TYPEOF__
	#define TOON_TYPEOF(X) __typeof__((X))
#elif defined TOON_TYPEOF_BOOST
    #include <boost/typeof/typeof.hpp>
	#define TOON_TYPEOF(X) BOOST_TYPEOF((X))
#else
	#include <complex>
	namespace TooN{
		namespace Internal{
			#include <TooN/internal/builtin_typeof.h>
		}
	}
	#define TOON_TYPEOF(X) typename Internal::DeEnumerate<sizeof Internal::enumerate(X)>::type
#endif
