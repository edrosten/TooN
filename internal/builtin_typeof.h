//Copyright (C) Edward Rosten 2009, Gerhard Reitmayr 2011

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

template<int N> struct Enumerate{char i[N];};
Enumerate<0> enumerate(const unsigned char&);
Enumerate<1> enumerate(const char&);
Enumerate<2> enumerate(const int&);
Enumerate<3> enumerate(const unsigned int&);
Enumerate<4> enumerate(const float&);
Enumerate<5> enumerate(const double&);
Enumerate<6> enumerate(const std::complex<float>&);
Enumerate<7> enumerate(const std::complex<double>&);
template<int N> struct DeEnumerate{};
template<> struct DeEnumerate<0>{typedef unsigned char type;};
template<> struct DeEnumerate<1>{typedef char type;};
template<> struct DeEnumerate<2>{typedef int type;};
template<> struct DeEnumerate<3>{typedef unsigned int type;};
template<> struct DeEnumerate<4>{typedef float type;};
template<> struct DeEnumerate<5>{typedef double type;};
template<> struct DeEnumerate<6>{typedef std::complex<float> type;};
template<> struct DeEnumerate<7>{typedef std::complex<double> type;};

#ifdef _FADBAD_H
Enumerate<8> enumerate(const ::fadbad::F<float>&);
Enumerate<9> enumerate(const ::fadbad::F<double>&);
template<> struct DeEnumerate<8>{typedef ::fadbad::F<float> type;};
template<> struct DeEnumerate<9>{typedef ::fadbad::F<double> type;};
#endif
