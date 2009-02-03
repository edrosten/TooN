//-*- c++ -*-
#ifndef TOON_INCLUDE_TOON_H
#define TOON_INCLUDE_TOON_H
#include <iostream>
#include <cstdlib>
namespace TooN
{
	
	static const unsigned int max_bytes_on_stack=1000;
	struct Slicing{};


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
	
	#include <TooN/internal/allocator.hh>
	#include <TooN/internal/size_mismatch.hh>
	#include <TooN/internal/slice_error.hh>
	#include <TooN/internal/debug.hh>
	#include <TooN/internal/vbase.hh>
	#include <TooN/internal/vector.hh>

	#include <TooN/internal/mbase.hh>
	#include <TooN/internal/matrix.hh>

	#include <TooN/internal/make_vector.hh>
	#include <TooN/operators.h>
}
#endif
