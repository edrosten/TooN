//-*- c++ -*-
#ifndef TOON_INCLUDE_TOON_H
#define TOON_INCLUDE_TOON_H
#include <iostream>
#include <cstdlib>
#include <limits>
#include <new>
#include <TooN/internal/config.hh>
#include <TooN/internal/typeof.hh>

namespace TooN
{
	

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

		struct SpecifySize{};
		struct SpecifyStride{};
		struct SpecifyRows{};
		struct SpecifyCols{};

		struct NoIgnore{};
		
		//Some helper classes to label integer data with a specific meaning,
		//for creating slices of external data.
		struct Spec_____{                 };
		struct Spec____C{int           cs;};
		struct Spec___R_{int       rs    ;};
		struct Spec___RC{int       rs, cs;};
		struct Spec__C__{int    c        ;};
		struct Spec__C_C{int    c    , cs;};
		struct Spec__CR_{int    c, rs    ;};
		struct Spec__CRC{int    c, rs, cs;};
		struct Spec_R___{int r           ;};
		struct Spec_R__C{int r,        cs;};
		struct Spec_R_R_{int r,    rs    ;};
		struct Spec_R_RC{int r,    rs, cs;};
		struct Spec_RC__{int r, c        ;};
		struct Spec_RC_C{int r, c,     cs;};
		struct Spec_RCR_{int r, c, rs    ;};
		struct Spec_RCRC{int r, c, rs, cs;};
	}

	template<int Size, class Precision, class Base> struct Vector;
	template<int Rows, int Cols, class Precision, class Base> struct Matrix;
	template<typename T> class Operator;
	template<typename T> class SliceSpec: public T{};

	static const int Dynamic = -1;

	#include <TooN/internal/allocator.hh>

	#include <TooN/internal/size_mismatch.hh>
	#include <TooN/internal/slice_error.hh>
	#include <TooN/internal/debug.hh>

	#include <TooN/internal/vbase.hh>
	#include <TooN/internal/vector.hh>
	
	#include <TooN/internal/mbase.hh>
	#include <TooN/internal/matrix.hh>
	#include <TooN/internal/foreign_matrix.hh>

	#include <TooN/internal/make_vector.hh>
	#include <TooN/internal/operators.hh>
}
#endif
