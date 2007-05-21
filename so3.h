
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
#include <TooN/LU.h>
#include <TooN/helpers.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

    class SE3;

class SO3 {
public:
  friend std::istream& operator>>(std::istream& is, SO3& rhs);
  friend std::istream& operator>>(std::istream& is, class SE3& rhs);
  friend class SE3;
  inline SO3();
  inline SO3(const Matrix<3>& rhs);

  inline SO3& operator=(const Matrix<3>& rhs);
  template <class Accessor> static inline void coerce(FixedMatrix<3,3,Accessor>& M);

  template<class Accessor> inline static SO3 exp(const FixedVector<3,Accessor>& vect);
  template <class Accessor> inline static double exp_with_half(const FixedVector<3,Accessor>& vect, SO3& first, SO3& second, double& shtot);

  inline Vector<3> ln() const;

  inline double operator[](int i){return my_matrix[i/3][i%3];}

  inline SO3 inverse() const;

  inline SO3& operator *=(const SO3& rhs){
    my_matrix=my_matrix*rhs.my_matrix;
    return *this;
  }

  inline SO3 operator *(const SO3& rhs) const {
      return SO3(*this,rhs);
  }

  inline const Matrix<3>& get_matrix()const {return my_matrix;}
  // returns i-th generator times pos
  inline static Vector<3> generator_field(int i, Vector<3> pos);

  // adjoint transformation on the Lie algebra
  inline Vector<3> adjoint(Vector<3> vect) const ;

 private:
  struct Invert {};
  inline SO3(const SO3& so3, const Invert&) : my_matrix(so3.my_matrix.T()) {}
  inline SO3(const SO3& a, const SO3& b) : my_matrix(a.my_matrix*b.my_matrix) {}
      
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

template <class Accessor> inline double SO3::exp_with_half(const FixedVector<3,Accessor>& vect, SO3& first, SO3& second, double& shtot){
     const double thetasq = vect*vect;
     const double theta = sqrt(thetasq);
     double stot, sqtoth, cth, ct;
     if (thetasq < 1e-6) {
	 const double sixth_thetasq = thetasq/6.0;
	 stot = 1 - sixth_thetasq;
	 shtot = 0.5 - sixth_thetasq*0.125;
	 sqtoth = 0.5 - sixth_thetasq*0.03125;
	 cth = 1.0 - thetasq*0.125;
	 ct = 1.0 - thetasq*0.5;
     } else {
	 const double thetainv = 1.0/theta;
	 const double ht = theta*0.5;
	 const double sth = sin(ht);
	 cth = cos(ht);
	 ct = 2*cth*cth - 1;
	 shtot = sth*thetainv;
	 stot = 2*cth*shtot;
	 sqtoth = 2*sqrt(0.5*(1 - cth))*thetainv;
     }

     {
	 const double a = shtot*vect[0];
	 const double b = shtot*vect[1];
	 const double c = shtot*vect[2];
	 first.my_matrix[0][0] = ct + 2*a*a;
	 first.my_matrix[1][1] = ct + 2*b*b;
	 first.my_matrix[2][2] = ct + 2*c*c;

	 const double f = stot*vect[2];
	 const double cr1 = 2*a*b;
	 first.my_matrix[1][0] = f + cr1;
	 first.my_matrix[0][1] = cr1 - f;
	 
	 const double e = stot*vect[1];
	 const double cr2 = 2*a*c;
	 first.my_matrix[0][2] = e + cr2;
	 first.my_matrix[2][0] = cr2 - e;
	 	 
	 const double d = stot*vect[0];
	 const double cr3 = 2*b*c;
	 first.my_matrix[2][1] = d + cr3;
	 first.my_matrix[1][2] = cr3 - d;
     }

     {
	 const double a = sqtoth*vect[0];
	 const double b = sqtoth*vect[1];
	 const double c = sqtoth*vect[2];
	 second.my_matrix[0][0] = cth + 0.5*a*a;
	 second.my_matrix[1][1] = cth + 0.5*b*b;
	 second.my_matrix[2][2] = cth + 0.5*c*c;

	 const double f = shtot*vect[2];
	 const double cr1 = 0.5*a*b;
	 second.my_matrix[1][0] = f + cr1;
	 second.my_matrix[0][1] = cr1 - f;
	 
	 const double e = shtot*vect[1];
	 const double cr2 = 0.5*a*c;
	 second.my_matrix[0][2] = e + cr2;
	 second.my_matrix[2][0] = cr2 - e;
	 	 
	 const double d = shtot*vect[0];
	 const double cr3 = 0.5*b*c;
	 second.my_matrix[2][1] = d + cr3;
	 second.my_matrix[1][2] = cr3 - d;
     }
     return theta;
}

inline Vector<3> SO3::ln() const{
  Vector<3> result;

  double trace = my_matrix[0][0] + my_matrix[1][1] + my_matrix[2][2];
 
  // 2 * cos(theta) - 1 == trace

  result[0] = (my_matrix[2][1]-my_matrix[1][2])/2;
  result[1] = (my_matrix[0][2]-my_matrix[2][0])/2;
  result[2] = (my_matrix[1][0]-my_matrix[0][1])/2;

  if (trace > -0.95) {
      double sin_angle_abs = sqrt(result*result);
      if (sin_angle_abs > 0.00001) {
	  double angle;
	  if(sin_angle_abs > 1.0) {
	      sin_angle_abs = 1.0;
	      angle = M_PI_2;
	  } else {
	      angle = asin(sin_angle_abs);
	      if (trace < 1)
		  angle = M_PI - angle;
	  }
	  result *= angle / sin_angle_abs;
      }
      return result;
  } else {
      TooN::Matrix<3> A = my_matrix;
      A[0][0] -= 1.0;
      A[1][1] -= 1.0;
      A[2][2] -= 1.0;
      TooN::LU<3> lu(A);     
      const TooN::Matrix<3,3,TooN::ColMajor>& u = lu.get_lu().T();
      const double u0 = fabs(u[0][0]);
      const double u1 = fabs(u[1][1]);
      const double u2 = fabs(u[2][2]);
      int row;      
      if (u0 < u1)
	  row = u0 < u2 ? 0 : 2;
      else
	  row = u1 < u2 ? 1 : 2;      
      //std::cerr << u << std::endl;
      //std::cerr << "row = " << row << std::endl;
      TooN::Vector<3> r;
      const double factor = acos(0.5*(trace-1));
      switch (row) {
      case 0: r[0] = factor; r[1] = 0; r[2] = 0; break;
      case 1: r[0] = u[0][1]*u[2][2]; r[1] = -u[0][0]*u[2][2]; r[2] = 0; r *= factor/sqrt(r*r); break;
      case 2: r[0] = u[0][1]*u[1][2] - u[0][2]*u[1][1]; r[1] = -u[0][0]*u[1][2]; r[2] = u[0][0]*u[1][1];  r *= factor/sqrt(r*r); break;
      }
      if (result * r < 0)
	  return -r;
      else
	  return r;
  }
}

inline SO3 SO3::inverse() const{
    return SO3(*this, Invert());
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
