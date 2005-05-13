//-*- c++ -*-

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


#ifndef __NUMHELPERS_H
#define __NUMHELPERS_H

#include <TooN/toon.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif 

//////////// CONTENTS OF THIS FILE ///////////

// normalizations (note US spelling)
template <class Accessor> inline void  normalize(VectorBase<Accessor>& v);
template <class Accessor> inline void  normalize_last(VectorBase<Accessor>& v);
template <class Accessor> inline void  normalize_but_last(VectorBase<Accessor>& v);

// Project
template <int Size, class Accessor> Vector<Size-1>  project(const FixedVector<Size,Accessor>& v);
template  <class Accessor> Vector<>                 project(const DynamicVector<Accessor>& v);

// Unproject
template <int Size, class Accessor> Vector<Size+1>  unproject(const FixedVector<Size,Accessor>& v);
template  <class Accessor> Vector<>                 unproject(const DynamicVector<Accessor>& v);


// as_vector
template<int Size> inline FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(double* data);
template<int Size> inline const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(const double* data);

// set a matrix to (a multiple of) the Identity
template <int Size, class Accessor> void  Identity(FixedMatrix<Size,Size,Accessor>&m, const double factor=1);

// symmetrize a matrix
template <int Size, class Accessor> void Symmetrize(FixedMatrix<Size,Size,Accessor>& m);

// transpose a matrix
template <int Size, class Accessor> void Transpose(FixedMatrix<Size,Size,Accessor>& m);

// set a vector to zero
template <class Accessor> inline void Zero(VectorBase<Accessor>& v);

//////////////////////////////////////////////



// normalizations (note US spelling)
template <class Accessor>
inline void normalize(VectorBase<Accessor>& v){
  double sumsq=0;
  for(int i=0; i<v.size(); i++){
    sumsq += v[i] * v[i];
  }
  const double scalefac = 1/sqrt(sumsq);
  for(int i=0; i<v.size(); i++){
    v[i]*=scalefac;
  }
}

template <class Accessor>
inline void normalize_last(VectorBase<Accessor>& v){
  const double scalefac = 1/v[v.size()-1];
  for(int i=0; i<v.size(); i++){
    v[i]*=scalefac;
  }
}

template <class Accessor>
inline void normalize_but_last(VectorBase<Accessor>& v){
  double sumsq=0;
  for(int i=0; i<v.size()-1; i++){
    sumsq += v[i] * v[i];
  }
  const double scalefac = 1/sqrt(sumsq);
  for(int i=0; i<v.size(); i++){
    v[i]*=scalefac;
  }
}


// Project
template <int Size, class Accessor>
struct FixedVProject {
  inline static void eval(Vector<Size-1>& ret, const FixedVector<Size,Accessor>& v){
    const double scalefac = 1/v[Size-1];
    for(int i=0; i<Size-1; i++){
      ret[i]=v[i]*scalefac;
    }
  }
};

template <int Size, class Accessor>
Vector<Size-1> project(const FixedVector<Size,Accessor>& v){
  return Vector<Size-1>(v,Operator<FixedVProject<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVProject : public VSizer{
  inline static void eval(Vector<>& ret, const DynamicVector<Accessor>& v){
    const int size = v.size();
    set_size(ret,size-1);
    const double scalefac = 1/v[size-1];
    for(int i=0; i<size-1; i++){
      ret[i]=v[i]*scalefac;
    }
  }
};

template  <class Accessor>
Vector<> project(const DynamicVector<Accessor>& v){
  return Vector<>(v,Operator<DynamicVProject<Accessor> >());
}

// Unproject
template <int Size, class Accessor>
struct FixedVUnproject {
  inline static void eval(Vector<Size+1>& ret, const FixedVector<Size,Accessor>& v){
    ret.template slice<0,Size>() = v;
    ret[Size]=1;
  }
};

template <int Size, class Accessor>
Vector<Size+1> unproject(const FixedVector<Size,Accessor>& v){
  return Vector<Size+1>(v,Operator<FixedVUnproject<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVUnproject : public VSizer{
  inline static void eval(Vector<>& ret, const DynamicVector<Accessor>& v){
    const int size = v.size();
    set_size(ret,size+1);
    v.copy_into(ret.get_data_ptr());
    ret[size]=1;
  }
};

template  <class Accessor>
Vector<> unproject(const DynamicVector<Accessor>& v){
  return Vector<>(v,Operator<DynamicVUnproject<Accessor> >());
}


// as_vector<Size>(double*) to convert a pointer to
// an array of data into a Vector<Size>
template<int Size>
inline FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(double* data){
  return reinterpret_cast<FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&>(*data);
}

template<int Size>
inline const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&  as_vector(const double* data){
  return reinterpret_cast<const FixedVector<Size,FixedVAccessor<Size,Stack<Size> > >&>(*data);
}


// set a matrix to (a multiple of) the Identity
template <int Size, class Accessor>
void Identity(FixedMatrix<Size,Size,Accessor>&m, const double factor){
  for(int i=0; i<Size; i++){
    for(int j=0; j<Size; j++){
      m(i,j)=0;
    }
    m(i,i)=factor;
  }
};


// symmetrize a matrix
template <int Size, class Accessor>
void Symmetrize(FixedMatrix<Size,Size,Accessor>& m){
  for(int r=0; r<Size-1; r++){
    for(int c=r; c<Size; c++){
      m(c,r) = m(r,c) = 0.5*(m(r,c)+m(c,r));
    }
  }
}

// Transpose a matrix
template<int Size, class Accessor>
void Transpose(FixedMatrix<Size,Size,Accessor>& m){
  for(int r=0; r<Size-1; r++){
    for(int c=r; c<Size; c++){
      double temp = m(r,c);
      m(r,c) = m(c,r);
      m(c,r) = temp;
    }
  }
}

// set a vector to zero
template <class Accessor> inline void Zero(VectorBase<Accessor>& v){
  for(int i=0; i<v.size(); i++){
    v[i]=0;
  }
}

#ifndef TOON_NO_NAMESPACE
}
#endif 

#endif
