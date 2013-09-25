// -*- c++ -*-

// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

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

namespace TooN {

// class to generate compile time error
// general case which doesn't exist
template<int Size1, int Size2>
struct SizeMismatch_;

// special cases which do exist
template<int Size>
struct SizeMismatch_<Size,Size>{
  static inline void test(int, int){}
};

template<int Size>
struct SizeMismatch_<Dynamic,Size>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

template<int Size>
struct SizeMismatch_<Size,Dynamic>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

template <>
struct SizeMismatch_<Dynamic,Dynamic>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

#if 0
namespace Internal
{
	struct BadSize;
}
#endif

#ifdef TOON_TEST_INTERNALS
template<int Size1, int Size2>
struct SizeMismatch_
{
	static inline void test(int, int)
	{
		throw Internal::StaticSizeMismatch();
	}
};
#endif

template<int Size1, int Size2>
struct SizeMismatch
{
	static inline void test(int s1, int s2)
	{
		SizeMismatch_< (Size1 == Dynamic || Size1 == Resizable)?Dynamic:Size1,
		               (Size2 == Dynamic || Size2 == Resizable)?Dynamic:Size2 >::test(s1, s2);
	}
};

}
