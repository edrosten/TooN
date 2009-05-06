//-*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk), Gerhard Reitmayr (gr281@cam.ac.uk)
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


#ifndef TOON_INCLUDE_TOON_H
#define TOON_INCLUDE_TOON_H
#include <iostream>
#include <cstdlib>
#include <limits>
#include <new>
#include <TooN/internal/config.hh>
#include <TooN/internal/typeof.hh>

#ifdef TOON_INITIALIZE_RANDOM
#include <ctime>
#endif

namespace TooN {

#ifdef TOON_TEST_INTERNALS
	namespace Internal
	{
		struct BadIndex{};
		struct SliceError{};
		struct StaticSliceError{};
		struct SizeMismatch{};
		struct StaticSizeMismatch{};

	}
#endif
	
	//Is the number a field? ie, *, -, *, / defined.
	//Specialize this to make TooN work properly with new types
	using std::numeric_limits;
	template<class C> struct IsField
	{
		static const int value = numeric_limits<C>::is_specialized;
	};
	
	namespace Internal
	{
		static const unsigned int max_bytes_on_stack=1000;
		struct Slicing{};
		template<int RowStride, int ColStride> struct Slice;
	}

	template<int Size, class Precision, class Base> struct Vector;
	template<int Rows, int Cols, class Precision, class Base> struct Matrix;
	template<int Size, class Precision, class Base> struct DiagonalMatrix;
	template<typename T> struct Operator;
	
	static const int Dynamic = -1;

	typedef double DefaultPrecision;
}

#include <TooN/internal/allocator.hh>

#include <TooN/internal/size_mismatch.hh>
#include <TooN/internal/slice_error.hh>
#include <TooN/internal/debug.hh>

#include <TooN/internal/vbase.hh>
#include <TooN/internal/vector.hh>
	
#include <TooN/internal/mbase.hh>
#include <TooN/internal/matrix.hh>
#include <TooN/internal/reference.hh>

#include <TooN/internal/make_vector.hh>
#include <TooN/internal/operators.hh>
	
#include <TooN/internal/objects.h>

#include <TooN/internal/diagmatrix.h>

#endif
