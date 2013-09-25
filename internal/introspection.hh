// -*- c++ -*-

// Copyright (C) 2012
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

// Allocators always copy data on copy construction.
//
// When a Vector/Matrix is constructed from a different, but compatible type
// copying is done at a much higher level: above the level that knows how the
// data is laid out in memory.
//
// At this level, copy construction is required since it is only known here 
// whether data or a reference to data should be copied.


///Pretty generic SFINAE introspection generator
///@ingroup gInternal
namespace TooN{
namespace Internal{
//Fake function to pretend to return an instance of a type
template<class C>
const C& get();

//Two typedefs guaranteed to have different sizes
typedef char OneSized[1];
typedef char TwoSized[2];


//Class to give us a size 2 return value or fail on 
//substituting S
template<int S> 
struct SFINAE_dummy
{
        typedef TwoSized Type;
};


#define TOON_CREATE_THING_DETECTOR(Y, X, Z) \
template<class C>\
struct Has_##X##_##Z\
{\
        static OneSized& detect_##X##_##Z(...);\
\
        template<class S>\
        static typename SFINAE_dummy<sizeof(Y)>::Type& detect_##X##_##Z(const S&);\
\
        static const bool Has = sizeof(detect_##X##_##Z(get<C>())) == 2;\
}

#define TOON_CREATE_MEMBER_DETECTOR(X)  TOON_CREATE_THING_DETECTOR(S::X, X, Member)
#define TOON_CREATE_METHOD_DETECTOR(X)  TOON_CREATE_THING_DETECTOR(&S::X, X, Method)
#define TOON_CREATE_TYPEDEF_DETECTOR(X) TOON_CREATE_THING_DETECTOR(typename S::X, X, Typedef)

}}
