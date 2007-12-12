
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
#ifndef __VBASE_HH
#define __VBASE_HH

// VectorBase //

template <class Accessor>
class VectorBase : public Accessor {
public:
  const double* get_data_ptr()const{return Accessor::my_values;}
  double* get_data_ptr(){return Accessor::my_values;}

};

// operator ostream& <<
template <class Accessor>
std::ostream& operator << (std::ostream& os, const VectorBase<Accessor>& v){
  for (int i=0; i<v.size(); i++){
    os << std::setw(12) << v[i] << "  ";
  }
  return os;
}

// operator istream& >>
template <class Accessor>
std::istream& operator >> (std::istream& is, VectorBase<Accessor>& v){
  for (int i=0; i<v.size(); i++){
    is >>  v[i];
  }
  return is;
}

template <class Accessor1, class Accessor2> struct VectorCopy;

////////////////////////////////////////////////////////
//                                                    //
// Fixed and Dynamic Vector classes                   //
// all the arithmetic and assignment                  //
// operations are applied to these                    //
//                                                    //
////////////////////////////////////////////////////////

template <int Size, class Accessor>
class FixedVector : public VectorBase<Accessor> {
public:
  // assignment from correct sized FixedVector
  template<class Accessor2>
  inline FixedVector& operator=(const FixedVector<Size,Accessor2>& from){
    VectorCopy<Accessor, Accessor2>::eval(*this,from);
    return *this;
  }

  // copy assignment
  inline FixedVector& operator=(const FixedVector& from){
    VectorCopy<Accessor,Accessor>::eval(*this,from);
    return *this;
  }

  // assignment from any DynamicVector
  template<class Accessor2>
    inline FixedVector& operator=(const DynamicVector<Accessor2>& from){
    assert(from.size() == Size);
    VectorCopy<Accessor, Accessor2>::eval(*this, from);
    return *this;
  }

  // assignment from a double - uses vector magic
  VectorMagic::VectorFiller<1,Size, FixedVector<Size, Accessor>, VectorMagic::CommaStyle> operator=(double t) {
    (*this)[0] = t;
    return VectorMagic::VectorFiller<1,Size, FixedVector<Size, Accessor>, VectorMagic::CommaStyle>(*this);
  }

  // insertion operators - uses vector magic
  VectorMagic::VectorFiller<1,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle> operator<<(double t) {
    (*this)[0] = t;
    return VectorMagic::VectorFiller<1,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle>(*this);
  }

  template <int N> VectorMagic::VectorFiller<N,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle> operator<<(const VectorMagic::ComponentPlaceHolder<N>& t) {
    return VectorMagic::VectorFiller<N,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle>(*this);
  }

  template <int N> VectorMagic::VectorFiller<N,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle> operator<<(const Vector<N>& t) {
    (*this).template slice<0,N>() = t;
    return VectorMagic::VectorFiller<N,Size, FixedVector<Size,Accessor>, VectorMagic::InsertionStyle>(*this);
  }

    static void dummy() {}

  template<class A, int I> FixedVector<Size, Accessor>& operator=(const VectorMagic::VectorCreator<A,I>& v)
  {
    v.assign(*this);
    return *this;
  }

};

template <class Accessor>
class DynamicVector : public VectorBase<Accessor>{
  typedef VectorBase<Accessor> parent;

public:
  // assignment from any VectorBase
  template<class Accessor2>
  DynamicVector& operator=(const VectorBase<Accessor2>& from){
    assert(parent::my_size == from.size());
    VectorCopy<Accessor,Accessor2>::eval(*this,from);
    return *this;
  }

  // repeated for explicit copy assignment
  DynamicVector& operator=(const DynamicVector& from){
    assert(parent::my_size == from.size());
    VectorCopy<Accessor,Accessor>::eval(*this,from);
    return *this;
  }
    operator DynamicVector& () { return *this; }
    operator const DynamicVector& () const { return *this; }

  template<class A, int I> DynamicVector<Accessor> & operator=(const VectorMagic::VectorCreator<A,I>& v)
  {
    v.assign(*this);
    return *this;
  }

};


// Special kinds of DynamicVector only constructed internally
// e.g. from DynamicMAccessor<>::operator[]

template <class V> struct NonConst : public V {
  inline operator const V&() const { return *this; }
  inline operator V&() { return *this; }
  template <class T> inline NonConst& operator=(const T& t) { V::operator=(t);  return *this; }
};

#include <TooN/vaccessor.hh>


// Data copying //

template <class Accessor1, class Accessor2>
struct VectorCopy {
  inline static void eval(VectorBase<Accessor1>& to, const VectorBase<Accessor2>& from){
    for(int i=0; i<from.size(); i++){
      to[i]=from[i];
    }
  }
};

template <int Size, class Zone1, class Zone2>
struct VectorCopy<FixedVAccessor<Size,Zone1>,FixedVAccessor<Size,Zone2> > {
  inline static void eval(VectorBase<FixedVAccessor<Size,Zone1> >& to,
			  const VectorBase<FixedVAccessor<Size,Zone2> >& from){
    memcpy(to.get_data_ptr(),from.get_data_ptr(),Size*sizeof(double));
  }
};

template<>
struct VectorCopy<DynamicVAccessor,DynamicVAccessor>{
  inline static void eval(VectorBase<DynamicVAccessor>& to,
			  const VectorBase<DynamicVAccessor>& from){
    memcpy(to.get_data_ptr(),from.get_data_ptr(),from.size()*sizeof(double));
  }
};


#endif
