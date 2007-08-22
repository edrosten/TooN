/*
     Copyright (C) 2007 The Authors

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

/* This code mostly made by copying from so2.h !! */

#ifndef __SO2_H
#define __SO2_H

#include <TooN/TooN.h>
#include <TooN/LU.h>
#include <TooN/helpers.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

    class SE2;

class SO2 {
public:
  friend std::istream& operator>>(std::istream& is, SO2& rhs);
  friend std::istream& operator>>(std::istream& is, class SE2& rhs);
  friend class SE2;
  inline SO2();
  inline SO2(const Matrix<2>& rhs);

  inline SO2& operator=(const Matrix<2>& rhs);
  template <class Accessor> static inline void coerce(FixedMatrix<2,2,Accessor>& M);

  inline static SO2 exp(const double);

  inline double ln() const;

  inline double operator[](int i){return my_matrix[i/2][i%2];}

  inline SO2 inverse() const;

  inline SO2& operator *=(const SO2& rhs){
    my_matrix=my_matrix*rhs.my_matrix;
    return *this;
  }

  inline SO2 operator *(const SO2& rhs) const {
      return SO2(*this,rhs);
  }

  inline const Matrix<2>& get_matrix()const {return my_matrix;}
  // returns i-th generator times pos
  inline static Vector<2> generator(Vector<2> pos);

 private:
  struct Invert {};
  inline SO2(const SO2& so2, const Invert&) : my_matrix(so2.my_matrix.T()) {}
  inline SO2(const SO2& a, const SO2& b) : my_matrix(a.my_matrix*b.my_matrix) {}
      
  Matrix<2> my_matrix;
};

inline std::ostream& operator<< (std::ostream& os, const SO2& rhs){
  return os << rhs.get_matrix();
}

inline std::istream& operator>>(std::istream& is, SO2& rhs){
  return is >> rhs.my_matrix;
}

template<class Accessor> inline
Vector<2> operator*(const SO2& lhs, const FixedVector<2,Accessor>& rhs){
  return lhs.get_matrix() * rhs;
}

template<class Accessor> inline
Vector<2> operator*(const SO2& lhs, const DynamicVector<Accessor>& rhs){
  assert(rhs.size() == 2);
  return lhs.get_matrix() * rhs;
}

template<class Accessor> inline
Vector<2> operator*(const DynamicVector<Accessor>& lhs, const SO2& rhs){
  assert(lhs.size() == 2);
  return lhs * rhs.get_matrix();
}

template<class Accessor> inline
Vector<2> operator*(const FixedVector<2,Accessor>& lhs, const SO2& rhs){
  return lhs * rhs.get_matrix();
}

template <int RHS, class Accessor> inline
Matrix<2,RHS> operator*(const SO2& lhs, const FixedMatrix<2,RHS,Accessor>& rhs){
  return lhs.get_matrix() * rhs;
}

template <int LHS, class Accessor> inline
Matrix<LHS,2> operator*(const FixedMatrix<LHS,2,Accessor>& lhs, const SO2& rhs){
  return lhs * rhs.get_matrix();
}

namespace SO2static
{
  static double identity[4]={1,0,0,1};
}

inline SO2::SO2() :
  my_matrix(SO2static::identity)
{}

inline SO2::SO2(const Matrix<2>& rhs){
    *this = rhs;
}

inline SO2& SO2::operator=(const Matrix<2>& rhs){
  my_matrix = rhs;
  coerce(my_matrix);
  return *this;
}

template <class Accessor> inline void SO2::coerce(FixedMatrix<2,2, Accessor>& M){
  normalize(M[0]);
  M[1] -= M[0] * (M[0]*M[1]);
  normalize(M[1]);
}

inline SO2 SO2::exp(const double d){
  Matrix<2> m;
  m[0][0] = m[1][1] = cos(d);
  m[1][0] = sin(d);
  m[0][1] = -m[1][0];
  return SO2(m);
}

inline double SO2::ln() const{
  return atan2(my_matrix[1][0], my_matrix[0][0]);
}

inline SO2 SO2::inverse() const{
    return SO2(*this, Invert());
}

inline Vector<2> SO2::generator(Vector<2> pos){
  Vector<2> result;
  result[0] = pos[1];
  result[1] = -pos[0];
  return result;
}
 
#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
