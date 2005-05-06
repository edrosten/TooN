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


#ifndef __VBASE_HH
#define __VBASE_HH

// VectorBase //
#include <assert.h>

template <class Accessor>
struct VectorBase : public Accessor {
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


////////////////////////////////////////////////////////
//                                                    //
// Fixed and Dynamic Vector classes                   //
// all the arithmetic and assignment                  //
// operations are applied to these                    //
//                                                    //
////////////////////////////////////////////////////////




template <int Size, class Accessor>
struct FixedVector : public VectorBase<Accessor> {
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
};

template <class Accessor>
struct DynamicVector : public VectorBase<Accessor>{
  typedef VectorBase<Accessor> parent;
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


};


// Special kinds of DynamicVector only constructed internally
// e.g. from DynamicMAccessor<>::operator[]
struct RefVector : public DynamicVector<DynamicVAccessor>{
  RefVector(int size, double* ptr){
    my_size=size;
    my_values=ptr;
  }
};

class DynamicSkipAccessor{
 public:

  //CHECK THIS
  template<int Start, int Length>
  inline DynamicVector<DynamicSkipAccessor> slice() 
  {
    //DynamicSkipAccessors do not own memory, so destruction of one will not free it 
	
	DynamicVector<DynamicSkipAccessor> r;

	r.my_size = Length;
	r.my_skip = my_skip;
	r.my_values = my_values + my_skip * Start;
	
	return r;
  }

  template<int Start, int Length>
  inline const DynamicVector<DynamicSkipAccessor> slice() const 
  {
    //DynamicSkipAccessors do not own memory, so destruction of one will not free it 
	
	DynamicVector<DynamicSkipAccessor> r;
	r.my_size = Length;
	r.my_skip = my_skip;
	r.my_values = my_values + my_skip * Start;
	
	return r;
  }

  inline const double& operator[] (int i)const 
  {
    return my_values[i*my_skip];
  }

  inline double& operator[] (int i) 
  {
    return my_values[i*my_skip];	
  }
  
  inline int size() const 
  {
    return my_size;
  }

  inline RefSkipMatrix<ColMajor> as_row(); // implemented in linoperators.hh
  inline RefSkipMatrix<RowMajor> as_col(); //


 protected:
  int my_size;
  int my_skip;
  double* my_values;
};

// e.g. from SkipMAccessor<>::operator[]
struct RefSkipVector : public DynamicVector<DynamicSkipAccessor> {
  RefSkipVector(int size, int skip, double* ptr){
    my_size=size;
    my_skip=skip;
    my_values=ptr;
  }
};

#endif
