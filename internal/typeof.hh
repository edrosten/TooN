// -*- c++ -*-

// Copyright (C) 2009 Ed Rosten (er258@cam.ac.uk)

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



#ifdef TOON_TYPEOF_DECLTYPE
	#define TOON_TYPEOF(X) decltype((X))
#elif defined TOON_TYPEOF_TYPEOF
	#define TOON_TYPEOF(X) typeof((X))
#elif defined TOON_TYPEOF___TYPEOF__
	#define TOON_TYPEOF(X) __typeof__((X))
#elif defined TOON_TYPEOF_BOOST
    #include <boost/typeof/typeof.hpp>
	#define TOON_TYPEOF(X) BOOST_TYPEOF((X))
#elif (__cplusplus > 199711L ) && ! defined TOON_TYPEOF_BUILTIN
    #define TOON_TYPEOF(X) decltype((X))
#elif defined __GNUC__ && ! defined TOON_TYPEOF_BUILTIN
    #define TOON_TYPEOF(X) typeof((X))
#else
	#include <complex>
	namespace TooN{
		namespace Internal{
			#include <TooN/internal/builtin_typeof.h>
		}
	}
	#define TOON_TYPEOF(X) typename Internal::DeEnumerate<sizeof Internal::enumerate(X)>::type
#endif
