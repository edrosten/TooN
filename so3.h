
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
  template <class Accessor> static inline void coerce(FixedMatrix<3,3,Accessor>& M);

  static inline SO3 exp(const double* vect);

  template<class Accessor> inline static SO3 exp(const FixedVector<3,Accessor>& vect);
  template <class Accessor> inline static double SO3::exp_with_half(const FixedVector<3,Accessor>& vect, SO3& first, SO3& second);

  inline Vector<3> ln() const;

  inline double operator[](int i){return my_matrix[i/3][i%3];}

  inline SO3 inverse() const;

  inline SO3& operator *=(const SO3& rhs){
    my_matrix=my_matrix*rhs.my_matrix;
    return *this;
  }

  inline SO3 operator *(const SO3& rhs) const {
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
  assert(rhs.size() == 3);
  return lhs.get_matrix() * rhs;
}




template<class Accessor> inline
Vector<3> operator*(const DynamicVector<Accessor>& lhs, const SO3& rhs){
  assert(lhs.size() == 3);
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
  return *this;
}

inline SO3& SO3::operator=(const Matrix<3>& rhs){
  my_matrix = rhs;
  coerce(my_matrix);
  return *this;
}

 template <class Accessor> inline void SO3::coerce(FixedMatrix<3,3, Accessor>& M){
     normalize(M[0]);
     M[1] -= M[0] * (M[0]*M[1]);
     normalize(M[1]);
     M[2] -= M[0] * (M[0]*M[2]);
     M[2] -= M[1] * (M[1]*M[2]);
     normalize(M[2]);
}


template <class Accessor> inline SO3 SO3::exp(const FixedVector<3,Accessor>& vect){
  SO3 result;
  
  double theta = sqrt(vect*vect);
  
  double stot = 1;
  double shtot = 0.5;

  double ct=cos(theta);

  if(theta > 0.001) {
      stot = sin(theta)/theta;
      shtot = sin(theta*0.5)/theta;
  }

  result.my_matrix[0][0] = result.my_matrix[1][1] = result.my_matrix[2][2] = ct;
  result.my_matrix[0][1] = -(result.my_matrix[1][0] = stot*vect[2]);
  result.my_matrix[2][0] = -(result.my_matrix[0][2] = stot*vect[1]);
  result.my_matrix[1][2] = -(result.my_matrix[2][1] = stot*vect[0]);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      result.my_matrix[i][j] += 2*shtot*shtot*vect[i]*vect[j];
    }
  }
  return result;
}

 template <class Accessor> inline double SO3::exp_with_half(const FixedVector<3,Accessor>& vect, SO3& first, SO3& second){
     double theta = sqrt(vect*vect);

     double stot = 1;
     double shtot = 0.5;
     double sqtoth = 0.5;

     double cth = cos(theta*0.5);
     double ct = 2 * cth*cth - 1;
     
     if(theta > 0.001) {	 
	 double sth = sin(theta*0.5);
	 double st = 2*cth*sth;
	 double thetaInv = 1.0/theta;

	 stot = st*thetaInv;
	 shtot = sth*thetaInv;
	 if (theta > 0.0005)
	     sqtoth = 2*sin(theta*0.25)*thetaInv;
     }
     
     first.my_matrix[0][0] = first.my_matrix[1][1] = first.my_matrix[2][2] = ct;
     first.my_matrix[0][1] = -(first.my_matrix[1][0] = stot*vect[2]);
     first.my_matrix[2][0] = -(first.my_matrix[0][2] = stot*vect[1]);
     first.my_matrix[1][2] = -(first.my_matrix[2][1] = stot*vect[0]);

     second.my_matrix[0][0] = second.my_matrix[1][1] = second.my_matrix[2][2] = cth;
     second.my_matrix[0][1] = -(second.my_matrix[1][0] = shtot*vect[2]);
     second.my_matrix[2][0] = -(second.my_matrix[0][2] = shtot*vect[1]);
     second.my_matrix[1][2] = -(second.my_matrix[2][1] = shtot*vect[0]);
     
     for(int i=0; i<3; i++){
	 for(int j=0; j<3; j++){
	     first.my_matrix[i][j] += 2*shtot*shtot*vect[i]*vect[j];
	 }
     }
     
     for(int i=0; i<3; i++){
	 for(int j=0; j<3; j++){
	     second.my_matrix[i][j] += 0.5*sqtoth*sqtoth*vect[i]*vect[j];
	 }
     }
     return theta;
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
  } else {
    if(trace < 1){
      if(my_matrix[0][0] > 0){
        result[0] = M_PI;
      } else if(my_matrix[1][1] > 0){
        result[1] = M_PI;
      } else if(my_matrix[2][2] > 0){
        result[2] = M_PI;
      }
    }
  }


  result *=tost;

  return result;
}

inline SO3 SO3::inverse() const{
  SO3 result;
  result.my_matrix = my_matrix.T();
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
