
/*                       
	 Copyright (C) 2005 Tom Drummond

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __TOON_H
#define __TOON_H

#ifdef WIN32 // to get M_PI, etc.
#define _USE_MATH_DEFINES   
#endif

#include <string.h>  // for memcpy
#include <cmath>    // for sqrt
#include <iomanip>
#include <cassert>
#include <iostream> // for input and output of vectors and matrices

#include <TooN/util.h>
#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

static const int General=-1;


// actual values are needed in transpose operation
// which does nothing except return a ref to
// a matrix with the opposite layout
//enum{RowMajor=0,ColMajor=1};


// just placeholders
struct RowMajor{
  static bool is_rowmajor(){return true;}
};
struct ColMajor{
  static bool is_rowmajor(){return false;}
};

struct NUMERICS {
  // enum{Stack,Heap}; // allocation zone for fixed size stuff
  enum{Owner,Reference}; // ownership of variable size stuff

  // static const int DefaultLayout = RowMajor;

  typedef RowMajor DefaultLayout;

  // multiplication policy
  enum{BlasMult,CPPMult};

  // maximum no of mults in a matrix/vector * matrix/vector
  // product where we use inline c++
  // above this number we switch to using BLAS
  static const int MaxCPPMultCount = 100;

  // maximum no of doubles in an object before
  // we put it on the heap instead of the stack
  static const int MaxStackSize=100;  // ie 10x10 matrix
};


// forward declarations of all needed classes


/////////////// Memory Management //////////////////

// Fixed Size Memory Access
template<int Size>
class Stack;

template<int Size>
class Heap;

//////////// Vectors ///////////////////

template <int Size, class AllocZone>
class FixedVAccessor;

template <int Size, int Skip>
class SkipAccessor;

class DynamicVAccessor;
class RefSkipAccessor;

template <class Accessor>
class VectorBase;

template <int Size, class Accessor>
class FixedVector;

template <class Accessor>
class DynamicVector;

//class RefVector;
//class ConstRefVector;
//class RefSkipVector;
//class ConstRefSkipVector;

template <int Size=General>
class Vector;


/////////// Matrices ////////////

template <int Rows, int Cols, class Layout, class AllocZone>
class FixedMAccessor;

template<int Rows, int Cols, int Skip, class Layout>
class SkipMAccessor;

template<class Layout>
class DynamicMAccessor;

template <class Layout>
class RefSkipMAccessor;

template <class Accessor>
class MatrixBase;

template<int Rows, int Cols, class Accessor>
class FixedMatrix;

template<class Accessor>
class DynamicMatrix;

template<class Layout>
class RefMatrix;
template<class Layout>
class ConstRefMatrix;

template<class Layout>
class RefSkipMatrix;
template<class Layout>
class ConstRefSkipMatrix;

template <int Rows=General,
          int Cols=Rows,
          class Layout = typename NUMERICS::DefaultLayout >
class Matrix; 


// never actually use one of these
// they're just used as dummy arguments
// to destinguish between templated constructors
// most compilers should be able to optimise them out
template <class T>
class Operator {};  

template <int Size, int I>
class ZoneHandler;

template<int Size>
class ZoneHandler<Size,0> {
public:
  typedef Stack<Size> get_zone;
};

template<int Size>
class ZoneHandler<Size,1> {
public:
  typedef Heap<Size> get_zone;
};

template <int Size>
class SizeTraits : public ZoneHandler<Size,(Size > NUMERICS::MaxStackSize ? 1 : 0)> {
};


#ifdef TOON_DEBUG
  #define TOON_ASSERT(X,Y) if(!(X)) throw Y()
  #define TOON_THROW
  #include <TooN/accessorexceptions.hh>
#else
  #define TOON_ASSERT(X,Y)
  #define TOON_THROW throw()
#endif

#include <TooN/vmagic.hh>

#include <TooN/membase.hh>
#include <TooN/vbase.hh>
#include <TooN/vaccessor.hh>
#include <TooN/mbase.hh>
#include <TooN/maccessor.hh>
#include <TooN/vclasses.hh>
#include <TooN/mclasses.hh>
#include <TooN/blasoperators.hh>
#include <TooN/linoperators.hh>

 namespace util {
#include <TooN/generated.h>
 }

#ifndef TOON_NO_NAMESPACE
}
#endif

#ifdef TOON_USING_NAMESPACE
	using namespace TooN;
#endif

#endif

