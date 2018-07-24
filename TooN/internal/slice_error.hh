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

namespace TooN {

namespace Internal
{
	template<bool StaticBad> 
	struct BadSlice;
	
	///@internal
	///@brief A static slice is OK.
	///This class is used after it has been determined that a slice is OK.
	///It does nothing except provide a callable function. By contrast, 
	///if the slice is not OK, then the class is not specified and the function
	///is therefore not callable, and a compile error results.
	///@ingroup gInternal
	template<> 
	struct BadSlice<0>{
		static void check(){} ///<This function does nothing: it merely exists.
	};
	
	///@internal
	///@brief Check if a slice is OK.
	///This class is used to see if a slice is OK. It provides a
	///callable function which checks the run-time slice information.
	///If the compile time information is bad, then it will not compile
	///at all. Otherwise, the sizes are checked at run-time. The check 
	///will be optimized away if the sizes are known statically.
	///@ingroup gInternal
	template<int Size, int Start, int Length> 
	struct CheckSlice
	{
		
		///@internal
		///@brief choose a number statically or dynamically.
		template<int Num> struct N
		{
			static int n(int num)
			{
				return (Num==Dynamic||Num==Resizable)?num:Num;
			}
		};

		///@internal 
		///@brief Check the slice.
		///This is full static checking, which is stricter than 
		///mixed chacking. For instance, none of the slice parameters. This
		///should be used in addition to the other check function.
		///are allowed to be -1 (Dynamic).
		static void check()
		{
			//Sanity check all basic static sizes
			BadSlice<!(Size== Dynamic || Size==Resizable || Size > 0)>::check();
			BadSlice<!(Start >= 0)>::check();
			BadSlice<!(Length >= 0)>::check();
			BadSlice<(Size != Dynamic && Size != Resizable && (Start + Length > Size))>::check();
		}	

		///@internal 
		///@brief Check the slice.
		///The policy is that static sized where present are used.
		///However, for extra debugging one can test to see if the
		///static and dynamic sizes are mismatched.
		///any sense whatsoever.
		///@param size Vector size
		///@param start Start position of the slice
		///@param length Length of the slice.
		static void check(int size, int start, int length)
		{
			//Sanity check all basic static sizes
			BadSlice<!(Size   == Dynamic || Size==Resizable || Size > 0)>::check();
			BadSlice<!(Start  == Dynamic || Start >= 0)>::check();
			BadSlice<!(Length == Dynamic || Length >= 0)>::check();
			
			//We can make sure Length <= Size, even if Start is unknown
			BadSlice<(Size!=Dynamic && Size != Resizable &&  Length != Dynamic && Length > Size)>::check();
			
			//We can make sure Start < Size even if Length is unknown
			BadSlice<(Start != Dynamic && Size != Dynamic && Size != Resizable && Start > Size)>::check();

			BadSlice<(Size != Dynamic && Size != Resizable && Start != Dynamic && Length != Dynamic && Start + Length > Size)>::check();
			#ifndef TOON_NDEBUG_MISMATCH
				if(Start != Dynamic && Size != Resizable && Start != start)
				{
					std::cerr << "TooN slice: mismatch between static and dynamic start.\n";
					std::abort();
				}
				if(Length != Dynamic && Size != Resizable && Length != length)
				{
					std::cerr << "TooN slice: mismatch between static and dynamic length.\n";
					std::abort();
				}
				if(Size != Dynamic && Size != Resizable && Size != size)
				{
					std::cerr << "TooN slice: mismatch between static and dynamic size.\n";
					std::abort();
				}
			#endif
			if( N<Start>::n(start) + N<Length>::n(length) > N<Size>::n(size) || 
			   N<Start>::n(start) < 0 ||
			   N<Length>::n(length) < 0)
			{
				#ifdef TOON_TEST_INTERNALS
					throw Internal::SliceError();
				#elif !defined TOON_NDEBUG_SLICE
					std::cerr << "TooN slice out of range" << std::endl;
					std::abort();
				#endif
			}
		}
	};	

	#ifdef TOON_TEST_INTERNALS
		template<bool StaticBad> 
		struct BadSlice{
			static void check(){
				throw Internal::StaticSliceError();
			}
		};
	#endif
}

}
