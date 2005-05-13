//     Copyright (C) 2002 Tom Drummond (twd20@eng.cam.ac.uk)

//     This library is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation; either
//     version 2.1 of the License, or (at your option) any later version.

//     This library is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.

//     You should have received a copy of the GNU Lesser General Public
//     License along with this library; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef __VCLASSES_HH
#define __VCLASSES_HH

#include <assert.h>

///////////////////////////
//                       //
// Actual Vector classes //
//                       //
///////////////////////////

namespace VectorMagic {
  template <int N=1>
  struct ComponentPlaceHolder {
  };

  struct InsertionStyle {};
  struct CommaStyle {};

  template <int Index, int Limit, class Vec, class Style> struct VectorFiller;
  template <int Index, int Limit, class Vec> struct VectorFiller<Index,Limit,Vec,CommaStyle> {
    Vec& v;
    bool final_initializer_but_Vector_not_filled;
    inline VectorFiller(Vec& vec) : v(vec), final_initializer_but_Vector_not_filled(true) {}
    template <class T> inline VectorFiller<Index+1,Limit,Vec,CommaStyle> operator,(const T& t) {
      v[Index] = t;
      final_initializer_but_Vector_not_filled = false;
      return VectorFiller<Index+1,Limit,Vec,CommaStyle>(v);
    }
    template <int N> inline VectorFiller<Index+N,Limit,Vec,CommaStyle> operator,(const ComponentPlaceHolder<N>& ph) {
      final_initializer_but_Vector_not_filled = false;
      return (VectorFiller<Index+1,Limit,Vec,CommaStyle>(v), ComponentPlaceHolder<N-1>());
    }
    inline VectorFiller<Index+1,Limit,Vec,CommaStyle> operator,(const ComponentPlaceHolder<1>& ph) {
      final_initializer_but_Vector_not_filled = false;
      return VectorFiller<Index+1,Limit,Vec,CommaStyle>(v);
    }
    inline ~VectorFiller() {
      assert(!final_initializer_but_Vector_not_filled);
    }
    inline operator Vec () const { return v; }
  };
  
  template <int Index, int Limit, class Vec> struct VectorFiller<Index,Limit,Vec,InsertionStyle> {
    Vec& v;
    inline VectorFiller(Vec& vec) : v(vec){}
    template <class T> inline VectorFiller<Index+1,Limit,Vec,InsertionStyle> operator<<(const T& t) {
      v[Index] = t;
      return VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v);
    }
    template <int N> inline VectorFiller<Index+N,Limit,Vec,InsertionStyle> operator<<(const ComponentPlaceHolder<N>& ph) {
      return (VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v) << ComponentPlaceHolder<N-1>());
    }
    inline VectorFiller<Index+1,Limit,Vec,InsertionStyle> operator<<(const ComponentPlaceHolder<1>& ph) {
      return VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v);
    }
    inline operator Vec () const { return v; }
  };

  template <int Index, class Vec> struct VectorFiller<Index, Index, Vec, CommaStyle> {
    Vec& v;
    inline VectorFiller(Vec& vec) : v(vec) {}
    template <class T> inline void operator,(const T& t) const { too_many_elements_given(); }
    inline operator Vec () const { return v; }
  };

  template <int Index, class Vec> struct VectorFiller<Index, Index, Vec, InsertionStyle> {
    Vec& v;
    inline VectorFiller(Vec& vec) : v(vec) {}
    template <class T> inline void operator,(const T& t) const { too_many_elements_given(); }
    inline operator Vec () const { return v; }
  };

  template <class Left, int Size, class Val> struct VectorCreator
  {
    const Left& left;
    const Val& val;
    VectorCreator(const Left& l, const Val& v) : left(l), val(v) { }
    template <class T> VectorCreator<VectorCreator<Left,Size,Val>, Size+1, T> operator,(const T& t) const {
      return VectorCreator<VectorCreator<Left,Size,Val>, Size+1, T>(*this, t);
    }
    template <class V> void assign(V& v) const { v[Size-1] = val; left.assign(v); }
    operator Vector<Size> () const {
      Vector<Size> v;
      assign(v);
      return v;
    }
  };

  struct BaseVectorCreator
  {
    template <class T> inline VectorCreator<BaseVectorCreator, 1, T> operator,(const T& t) const {
      return VectorCreator<BaseVectorCreator, 1, T>(*this, t);
    }
    template <class V> inline void assign(V& v) const {}
  };
}

static VectorMagic::BaseVectorCreator make_Vector;

namespace VectorMagic 
{
  inline void dummy_make_Vector_user() { int i; make_Vector.assign(i); }
}

template <int N> VectorMagic::ComponentPlaceHolder<N> no_change() { return VectorMagic::ComponentPlaceHolder<N>(); }
inline VectorMagic::ComponentPlaceHolder<1> no_change() { return VectorMagic::ComponentPlaceHolder<1>(); }
  

template <int Size>
class Vector : public FixedVector<Size, FixedVAccessor<Size,typename SizeTraits<Size>::get_zone> > {
 public:
  // default constructor does nothing
  inline Vector(){}

  // constructor from c-style array
  inline Vector(const double from[Size]){
    memcpy(this->my_values,from,Size*sizeof(double));
  }

  // constructor from 1-ary operator
  template <class T, class Op>
    inline Vector(const T& arg, const Operator<Op>&){Op::eval(*this,arg);}

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
    inline Vector(const LHS& lhs, const RHS& rhs, const Operator<Op>&){Op::eval(*this,lhs,rhs);}

  // constructor from correct sized FixedVector
  template<class Accessor>
  inline Vector(const FixedVector<Size,Accessor>& from){
    FixedVector<Size, FixedVAccessor<Size,typename SizeTraits<Size>::get_zone> >::operator=(from);
  }

  // constructor from any DynamicVector
  template<class Accessor>
    inline Vector(const DynamicVector<Accessor>& from){
    assert(from.size() == Size);
    FixedVector<Size, FixedVAccessor<Size,typename SizeTraits<Size>::get_zone> >::operator=(from);
  }
  
  template <class Accessor> inline Vector<Size>& operator=(const FixedVector<Size,Accessor>& fv) {
    *this = Vector<Size>(fv);
    return *this;
  }

  template <class Accessor> inline Vector<Size>& operator=(const DynamicVector<Accessor>& dv) {
    *this = Vector<Size>(dv);
    return *this;
  }

  template <int N> VectorMagic::VectorFiller<N,Size, Vector<Size>,VectorMagic::CommaStyle> operator=(const VectorMagic::ComponentPlaceHolder<N>& t) {
    return VectorMagic::VectorFiller<N,Size, Vector<Size>,VectorMagic::CommaStyle>(*this);
  }

  template <class T> VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::CommaStyle> operator=(const T& t) {
    (*this)[0] = t;
    return VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::CommaStyle>(*this);
  }

  template <int N> VectorMagic::VectorFiller<N,Size, Vector<Size>, VectorMagic::InsertionStyle> operator<<(const VectorMagic::ComponentPlaceHolder<N>& t) {
    return VectorMagic::VectorFiller<N,Size, Vector<Size>, VectorMagic::InsertionStyle>(*this);
  }

  template <class T> VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::InsertionStyle> operator<<(const T& t) {
    (*this)[0] = t;
    return VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::InsertionStyle>(*this);
  }
  
};


template <>
class Vector<> : public DynamicVector<DynamicVAccessor> {
 public:
  Vector(int Size) {
    this->my_size=Size; this->my_values = new double[Size];
  }

  inline Vector(int Size, double* from){
    this->my_size=Size;
    this->my_values = new double[Size];
    memcpy(this->my_values,from,Size*sizeof(double));
  }
  
  inline ~Vector(){
    delete[] this->my_values;
  }

  // constructor from 1-ary operator
  template <class T, class Op>
  Vector(const T& arg, const Operator<Op>&){
    Op::eval(*this,arg);
  }

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  Vector(const LHS& lhs, const RHS& rhs, const Operator<Op>&){
    Op::eval(*this,lhs,rhs);
  }

  // constructor from any VectorBase
  template<class Accessor>
  Vector(const VectorBase<Accessor>& from){
    this->my_size = from.size();
    this->my_values = new double[this->my_size];
    DynamicVector<DynamicVAccessor>::operator=(from);
  }
  
  // trap copy constructor here
  Vector(const Vector& from){
    this->my_size = from.my_size;
    this->my_values = new double[this->my_size];
    memcpy(this->my_values,from.my_values,this->my_size*sizeof(double));
  }
};

struct VSizer{
  static inline void set_size(Vector<>& v, int size){
    v.my_size=size;
    v.my_values=new double[size];
  }
};


#endif
