
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
#ifndef __SO3_H
#define __SO3_H

#include <TooN/TooN.h>
#include <TooN/helpers.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

class SO3 {
public:
  friend std::istream& operator>>(std::istream& is, SO3& rhs);
  friend std::istream& operator>>(std::istream& is, class SE3& rhs);
  inline SO3();
  inline SO3(const Matrix<3>& rhs);

  inline SO3& operator=(const SO3& rhs);
  inline SO3& operator=(const Matrix<3>& rhs);
  inline void coerce();

  static inline SO3 exp(const double* vect);

  template<class Accessor>
  inline static SO3 exp(const FixedVector<3,Accessor>& vect) {return exp(vect.get_data_ptr());}

  inline Vector<3> ln() const;

  inline double operator[](int i){return my_matrix[i/3][i%3];}

  inline SO3 inverse() const;

  inline SO3& operator *=(const SO3& rhs){
    my_matrix=my_matrix*rhs.my_matrix;
    return *this;
  }

  inline SO3 operator *(const SO3& rhs) const{
    return SO3(*this)*=rhs;
  }

  inline const Matrix<3>& get_matrix()const {return my_matrix;}

  // returns i-th generator times pos
  inline static Vector<3> generator_field(int i, Vector<3> pos);

  // adjoint transformation on the Lie algebra
  inline Vector<3> adjoint(Vector<3> vect) const ;

 private:
  Matrix<3> my_matrix;
};

inline std::ostream& operator<< (std::ostream& os, const SO3& rhs){
  return os << rhs.get_matrix();
}

inline std::istream& operator>>(std::istream& is, SO3& rhs){
  return is >> rhs.my_matrix;
}




template<class Accessor> inline
Vector<3> operator*(const SO3& lhs, const FixedVector<3,Accessor>& rhs){
  return lhs.get_matrix() * rhs;
}

template<class Accessor> inline
Vector<3> operator*(const SO3& lhs, const DynamicVector<Accessor>& rhs){
  //FIXME: check size
  return lhs.get_matrix() * rhs;
}




template<class Accessor> inline
Vector<3> operator*(const DynamicVector<Accessor>& lhs, const SO3& rhs){
  //FIXME: check size
  return lhs * rhs.get_matrix();
}

template<class Accessor> inline
Vector<3> operator*(const FixedVector<3,Accessor>& lhs, const SO3& rhs){
  return lhs * rhs.get_matrix();
}

template <int RHS, class Accessor> inline
Matrix<3,RHS> operator*(const SO3& lhs, const FixedMatrix<3,RHS,Accessor>& rhs){
  return lhs.get_matrix() * rhs;
}

template <int LHS, class Accessor> inline
Matrix<LHS,3> operator*(const FixedMatrix<LHS,3,Accessor>& lhs, const SO3& rhs){
  return lhs * rhs.get_matrix();
}

namespace SO3static
{
  static double identity[9]={1,0,0,0,1,0,0,0,1};
}

inline SO3::SO3() :
  my_matrix(SO3static::identity)
{}

inline SO3::SO3(const Matrix<3>& rhs){
    *this = rhs;
}

inline SO3& SO3::operator=(const SO3& rhs){
  my_matrix = rhs.my_matrix;
  coerce();
  return *this;
}

inline SO3& SO3::operator=(const Matrix<3>& rhs){
  my_matrix = rhs;
  coerce();
  return *this;
}

inline void SO3::coerce(){
  normalize(my_matrix[0]);
  my_matrix[1] -= my_matrix[0] * (my_matrix[0]*my_matrix[1]);
  normalize(my_matrix[1]);
  my_matrix[2] -= my_matrix[0] * (my_matrix[0]*my_matrix[2]);
  my_matrix[2] -= my_matrix[1] * (my_matrix[1]*my_matrix[2]);
  normalize(my_matrix[2]);
}


inline SO3 SO3::exp(const double* vect){
  SO3 result;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      result.my_matrix[i][j] = 0;
    }
  }

  double theta = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);

  double stot = 1;
  double shtot = 0.5;

  double ct=cos(theta);

  if(theta > 0.001) {
    stot = sin(theta)/theta;
    shtot = sin(theta/2)/theta;
  }

  result.my_matrix[0][1] -= stot*vect[2];
  result.my_matrix[0][2] += stot*vect[1];
  result.my_matrix[1][0] += stot*vect[2];
  result.my_matrix[1][2] -= stot*vect[0];
  result.my_matrix[2][0] -= stot*vect[1];
  result.my_matrix[2][1] += stot*vect[0];

  result.my_matrix[0][0] += ct;
  result.my_matrix[1][1] += ct;
  result.my_matrix[2][2] += ct;

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      result.my_matrix[i][j] += 2*shtot*shtot*vect[i]*vect[j];
    }
  }
  return result;
}


inline Vector<3> SO3::ln() const{
  Vector<3> result;

  double trace = my_matrix[0][0] + my_matrix[1][1] + my_matrix[2][2];

  result[0] = (my_matrix[2][1]-my_matrix[1][2])/2;
  result[1] = (my_matrix[0][2]-my_matrix[2][0])/2;
  result[2] = (my_matrix[1][0]-my_matrix[0][1])/2;

  double sin_angle_abs = sqrt(result*result);

  // FIX PAS 10/5/02 added the min since was occasionally giving a fraction over 1.0
  if(sin_angle_abs > 1.0)
     sin_angle_abs = 1.0;
  // END FIX

  double tost=1;
  if (sin_angle_abs > 0.00001){
    double angle = asin(sin_angle_abs);
    if(trace < 1) {
      angle = M_PI - angle;
    }

    tost = angle / sin_angle_abs;
  }


  result *=tost;

  return result;
}

inline SO3 SO3::inverse() const{
  SO3 result;
  result.my_matrix[0][0] = my_matrix[0][0];
  result.my_matrix[0][1] = my_matrix[1][0];
  result.my_matrix[0][2] = my_matrix[2][0];
  result.my_matrix[1][0] = my_matrix[0][1];
  result.my_matrix[1][1] = my_matrix[1][1];
  result.my_matrix[1][2] = my_matrix[2][1];
  result.my_matrix[2][0] = my_matrix[0][2];
  result.my_matrix[2][1] = my_matrix[1][2];
  result.my_matrix[2][2] = my_matrix[2][2];
  return result;
}


// SO3& SO3::operator *=(const SO3& rhs){
//   *this = *this * rhs;
//   return *this;
// }

// SO3 SO3::operator *(const SO3& rhs)const{
//   SO3 result;
//   result.my_matrix = my_matrix * rhs.my_matrix;
//   return result;
// }

inline Vector<3> SO3::generator_field(int i, Vector<3> pos){
  Vector<3> result;
  result[i]=0;
  result[(i+1)%3] = - pos[(i+2)%3];
  result[(i+2)%3] = pos[(i+1)%3];
  return result;
}

inline Vector<3> SO3::adjoint(Vector<3> vect)const {
  return my_matrix * vect;
}

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
