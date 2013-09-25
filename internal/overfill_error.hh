//Copyright (C) Edward Rosten 2009, 2010

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

namespace TooN{
namespace Internal{

template<bool b> struct overfill;
template<> struct overfill<0>{};

template<int N, int Size> struct CheckOverFill
{
	static void check(int)
	{
		#ifdef TOON_TEST_INTERNALS
			if(N >= Size)
				throw StaticVectorOverfill();
		#else
			Internal::overfill<(N>=Size)> overfilled_vector;
		#endif
	};
};

template<int N> struct CheckOverFill<N, -1>
{
	static void check(int s)
	{
		#ifdef TOON_TEST_INTERNALS
			if(N >= s)
				throw VectorOverfill();
		#elif !defined TOON_NDEBUG_FILL
			if(N >= s)
			{
				std::cerr << "TooN overfilled vector" << std::endl;
				std::abort();
			}
		#endif
	};
};


template<int N, int R, int C, bool IsDynamic=(R==-1||C==-1)> struct CheckMOverFill
{
	static void check(int)
	{
		#ifdef TOON_TEST_INTERNALS
			if(N >= R*C)
				throw StaticMatrixOverfill();
		#else
			Internal::overfill<(N>=R*C)>();
		#endif
	}
};

template<int N, int R, int C> struct CheckMOverFill<N, R, C, 1>
{
	static void check(int s)
	{
		#ifdef TOON_TEST_INTERNALS
			if(N >= s)
				throw StaticMatrixOverfill();
		#else
			if(N >= s)
			{
				std::cerr << "TooN overfilled matrix" << std::endl;
				std::abort();
			}
		#endif
	}
};

}}
