
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

  inline double operator[](int i) const {return my_matrix[i/3][i%3];}

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

template <class A1, class A2>
inline void rodrigues_so3_exp(const TooN::FixedVector<3,A1>& w, const double A, const double B, TooN::FixedMatrix<3,3,A2>& R)
{
    {
	const double wx2 = w[0]*w[0];
	const double wy2 = w[1]*w[1];
	const double wz2 = w[2]*w[2];

	R[0][0] = 1.0 - B*(wy2 + wz2);
	R[1][1] = 1.0 - B*(wx2 + wz2);
	R[2][2] = 1.0 - B*(wx2 + wy2);
    }
    {
	const double a = A*w[2];
	const double b = B*(w[0]*w[1]);
	R[0][1] = b - a;
	R[1][0] = b + a;
    }
    {
	const double a = A*w[1];
	const double b = B*(w[0]*w[2]);
	R[0][2] = b + a;
	R[2][0] = b - a;
    }
    {
	const double a = A*w[0];
	const double b = B*(w[1]*w[2]);
	R[1][2] = b - a;
	R[2][1] = b + a;
    }
}

template <class Accessor>
inline SO3 SO3::exp(const FixedVector<3,Accessor>& w){
    static const double one_6th = 1.0/6.0;
    static const double one_20th = 1.0/20.0;

    SO3 result;

    const double theta_sq = w*w;
    const double theta = sqrt(theta_sq);
    double A, B;

    if (theta_sq < 1e-8) {
	A = 1.0 - one_6th * theta_sq;
	B = 0.5;
    } else {
	if (theta_sq < 1e-6) {
	    B = 0.5 - 0.25 * one_6th * theta_sq;
	    A = 1.0 - theta_sq * one_6th*(1.0 - one_20th * theta_sq);
	} else {
	    const double inv_theta = 1.0/theta;
	    A = sin(theta) * inv_theta;
	    B = (1 - cos(theta)) * (inv_theta * inv_theta);
	}
    }
    rodrigues_so3_exp(w, A, B, result.my_matrix);
    return result;
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
