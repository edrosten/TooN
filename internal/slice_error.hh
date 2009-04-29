// -*- c++ -*-

// Copyright (C) 2009 Ed Rosten (er258@cam.ac.uk)
//
// This file is part of the TooN Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

namespace TooN {

namespace Internal
{
	template<bool StaticBad> 
	struct BadSlice;

	template<> 
	struct BadSlice<0>{
		static void check(){}
	};

	template<int Size, int Start, int Length> 
	struct CheckStaticSlice
	{
		static void check(int/*size*/)
		{
			BadSlice<(Start < 0) || (Length < 1) || (Start+Length>Size)>::check();
		}
	};	
	
	template<int Start, int Length> struct CheckStaticSlice<-1, Start, Length>{
		static void check(int size)
		{
			BadSlice<(Start < 0) || (Length < 1)>::check();
			if(Start + Length > size)
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

	struct CheckDynamicSlice{
		static void check(int size, int start, int length){
			if(start < 0 || start + length > size)
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
