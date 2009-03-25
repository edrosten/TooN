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
#include <TooN/helpers.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

template<class Precision> class SO2;
template<class Precision> class SE2;

template<class Precision> inline std::ostream & operator<<(std::ostream &, const SO2<Precision> & );
template<class Precision> inline std::istream & operator>>(std::istream &, SO2<Precision> & );

template<class Precision = double>
class SO2 {
public:
  friend std::istream& operator>><Precision>(std::istream& is, SO2& rhs);
  friend std::ostream& operator<< <Precision>(std::ostream& is, const SO2& rhs);
  friend class SE2<Precision>;
  inline SO2(); 
  SO2(const Matrix<2,2,Precision>& rhs) { 
      *this = rhs;
  }
  SO2(const Precision l) {
    *this = exp(l);
  }
  
  template <int R, int C, class P, class A> inline SO2& operator=(const Matrix<R,C,P,A>& rhs);
  template <int R, int C, class P, class A> static inline void coerce(Matrix<R,C,P,A>& M);

  inline static SO2 exp(const Precision);

  Precision ln() const {
      return atan2(my_matrix[1][0], my_matrix[0][0]);
  }

  SO2 inverse() const {
    return SO2(*this, Invert());
  }

  SO2& operator *=(const SO2& rhs){
    my_matrix=my_matrix*rhs.my_matrix;
    return *this;
  }

  SO2 operator *(const SO2& rhs) const {
      return SO2(*this,rhs);
  }

  const Matrix<2,2,Precision>& get_matrix() const {return my_matrix;}
  
  /// returns generator matrix
  static Matrix<2,2,Precision> generator() {return my_generator;}

 private:
  struct Invert {};
  inline SO2(const SO2& so2, const Invert&) : my_matrix(so2.my_matrix.T()) {}
  inline SO2(const SO2& a, const SO2& b) : my_matrix(a.my_matrix*b.my_matrix) {}
      
  Matrix<2,2,Precision> my_matrix;
  static Matrix<2,2,Precision> my_generator;
};

template <class Precision>
inline std::ostream& operator<< (std::ostream& os, const SO2<Precision> & rhs){
  return os << rhs.get_matrix();
}

template <class Precision>
inline std::istream& operator>>(std::istream& is, SO2<Precision>& rhs){
  return is >> rhs.my_matrix;
}

template<int D, class P1, class PV, class Accessor> inline
Vector<2, typename Internal::MultiplyType<P1, PV>::type> operator*(const SO2<P1> & lhs, const Vector<D, PV, Accessor> & rhs){
  return lhs.get_matrix() * rhs;
}

template<int D, class P1, class PV, class Accessor>  inline
Vector<2, typename Internal::MultiplyType<P1,PV>::type> operator*(const Vector<D, PV, Accessor>& lhs, const SO2<P1> & rhs){
  return lhs * rhs.get_matrix();
}

template <int R, int C, class P1, class P2, class Accessor> inline
Matrix<2,C,typename Internal::MultiplyType<P1,P2>::type> operator*(const SO2<P1> & lhs, const Matrix<R,C,P2,Accessor>& rhs){
  return lhs.get_matrix() * rhs;
}

template <int R, int C, class P1, class P2, class Accessor> inline
Matrix<R,2,typename Internal::MultiplyType<P1,P2>::type> operator*(const Matrix<R,C,P1,Accessor>& lhs, const SO2<P2>& rhs){
  return lhs * rhs.get_matrix();
}

namespace SO2static
{
  static double identity[4] = {1,0,0,1};
  static double generator[4] = {0,-1,1,0};
}

template <class Precision> 
inline SO2<Precision>::SO2() //:
  //my_matrix(SO2static::identity,2,2, RowMajor::Layout<2,2,Precision>())
{}

template <class P1>
template <int R, int C, class P2, class A>
inline SO2<P1> & SO2<P1>::operator=(const Matrix<R,C,P2,A>& rhs){
  my_matrix = rhs;
  coerce(my_matrix);
  return *this;
}

template <class Precision>
template <int R, int C, class P, class A> 
inline void SO2<Precision>::coerce(Matrix<R,C,P,A>& M){
  SizeMismatch<2,R>::test(2, M.num_rows());
  SizeMismatch<2,C>::test(2, M.num_cols());
  normalize(M[0]);
  M[1] -= M[0] * (M[0]*M[1]);
  normalize(M[1]);
}

template <class Precision>
inline SO2<Precision> SO2<Precision>::exp(const Precision d){
  SO2<Precision> result;
  result.my_matrix[0][0] = result.my_matrix[1][1] = cos(d);
  result.my_matrix[1][0] = sin(d);
  result.my_matrix[0][1] = -result.my_matrix[1][0];
  return result;
}

template <class Precision> Matrix<2,2, Precision> SO2<Precision>::my_generator; //(SO2static::generator, 2, 2, RowMajor::Layout<2,2,Precision>());

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
