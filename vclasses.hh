
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
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __VCLASSES_HH
#define __VCLASSES_HH

///////////////////////////
//                       //
// Actual Vector classes //
//                       //
///////////////////////////

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

  // vector magic assignment operators (uses comma style) - insertion style defined in vbase.hh
  VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::CommaStyle> operator=(double t) {
    (*this)[0] = t;
    return VectorMagic::VectorFiller<1,Size, Vector<Size>, VectorMagic::CommaStyle>(*this);
  }

  template <int N> VectorMagic::VectorFiller<N,Size, Vector<Size>,VectorMagic::CommaStyle> operator=(const VectorMagic::ComponentPlaceHolder<N>& t) {
    return VectorMagic::VectorFiller<N,Size, Vector<Size>,VectorMagic::CommaStyle>(*this);
  }

};





template <>
class Vector<> : public DynamicVector<DynamicVAccessor> {
 public:
  Vector(){
  	this->my_size = 0;
	this->my_values = NULL;
  }

  inline void resize(int new_size)
  {
  	if(this->my_size != new_size)
	{
		delete[] this->my_values;
		this->my_size = new_size;
		this->my_values = new double[this->my_size];
	}
  }

    template <class Accessor> inline void assign(const VectorBase<Accessor>& v) {
	resize(v.size());
	DynamicVector<DynamicVAccessor>::operator=(v);
    }
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
