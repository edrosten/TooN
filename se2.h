
/*                       
     Copyright (C) 2005 The Authors

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

/* This code mostly made by copying from se3.h !! */

#ifndef __SE2_H
#define __SE2_H

#include <TooN/so2.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

class SE2 {
  friend SE2 operator*(const SO2& lhs, const SE2& rhs);
  friend std::istream& operator>> (std::istream& is, SE2& rhs);
  
 public:
  inline SE2();
  template <class A> inline SE2(const SO2& R, const FixedVector<2,A>& T) : my_rotation(R), my_translation(T) {}
      

  inline SO2& get_rotation(){return my_rotation;}
  inline const SO2& get_rotation() const {return my_rotation;}
  inline Vector<2>& get_translation() {return my_translation;}
  inline const Vector<2>& get_translation() const {return my_translation;}

  static inline SE2 exp(const Vector<3>& vect);
  static inline Vector<3> ln(const SE2& se2);
  inline Vector<3> ln() const { return SE2::ln(*this); }

  inline SE2 inverse() const;

  inline SE2& operator *=(const SE2& rhs);
  inline SE2 operator *(const SE2& rhs) const { return SE2(my_rotation*rhs.my_rotation, my_translation + my_rotation*rhs.my_translation); }
  inline SE2& left_multiply_by(const SE2& left);

  static inline Vector<3> generator_field(int i, Vector<3> pos);


  template<class Accessor>
  inline void adjoint(FixedVector<3,Accessor>& vect)const;

  template <class Accessor>
  inline void adjoint(FixedMatrix<3,3,Accessor>& M)const;

private:
  SO2 my_rotation;
  Vector<2> my_translation;
};


// left multiply an SE2 by an SO2 
inline SE2 operator*(const SO2& lhs, const SE2& rhs);


// transfers a vector in the Lie algebra
// from one coord frame to another
// so that exp(adjoint(vect)) = (*this) * exp(vect) * (this->inverse())
template<class Accessor>
  inline void SE2::adjoint(FixedVector<3,Accessor>& vect)const{
  vect.template slice<0,2>() = my_rotation * vect.template slice<0,2>();
  vect[0] += vect[2] * my_translation[1];
  vect[1] -= vect[2] * my_translation[0];
 }
 
template <class Accessor>
inline void SE2::adjoint(FixedMatrix<3,3,Accessor>& M)const{
  for(int i=0; i<3; i++){
    adjoint(M.T()[i]);
  }
  for(int i=0; i<3; i++){
    adjoint(M[i]);
  }
}


// operator ostream& <<
inline std::ostream& operator <<(std::ostream& os, const SE2& rhs){
  for(int i=0; i<2; i++){
    os << rhs.get_rotation().get_matrix()[i] << rhs.get_translation()[i] << std::endl;
  }
  return os;
}

// operator istream& >>
inline std::istream& operator>>(std::istream& is, SE2& rhs){
  for(int i=0; i<2; i++){
    is >> rhs.get_rotation().my_matrix[i] >> rhs.get_translation()[i];
  }
  return is;
}


//////////////////
// operator *   //
// SE2 * Vector //
//////////////////

template<class VectorType>
struct SE2VMult {
  inline static void eval(Vector<3>& ret, const SE2& lhs, const VectorType& rhs){
    ret.template slice<0,2>()=lhs.get_rotation()*rhs.template slice<0,2>();
    ret.template slice<0,2>()+=lhs.get_translation() * rhs[2];
    ret[2] = rhs[2];
  }
};

template<class Accessor> inline
Vector<3> operator*(const SE2& lhs, const FixedVector<3,Accessor>& rhs){
  return Vector<3>(lhs,rhs,Operator<SE2VMult<FixedVector<3, Accessor> > >());
}

template<class Accessor> inline
Vector<3> operator*(const SE2& lhs, const DynamicVector<Accessor>& rhs){
	//FIXME: size checking
  return Vector<3>(lhs,rhs,Operator<SE2VMult<DynamicVector<Accessor> > >());
}

template <class Accessor> inline
Vector<2> operator*(const SE2& lhs, const FixedVector<2,Accessor>& rhs){
  return lhs.get_translation() + lhs.get_rotation() * rhs;
}


//////////////////
// operator *   //
// Vector * SE2 //
//////////////////

template<class Accessor>
struct VSE2Mult {
  inline static void eval(Vector<4>& ret, const FixedVector<4,Accessor>& lhs, const SE2& rhs){
    ret.template slice<0,2>() = lhs.template slice<0,2>() * rhs.get_rotation();
    ret[2] = lhs[2];
    ret[2] += lhs.template slice<0,2>() * rhs.get_translation();
  }
};

template<class Accessor> inline
Vector<4> operator*(const FixedVector<4,Accessor>& lhs, const SE2& rhs){
  return Vector<4>(lhs,rhs,Operator<VSE2Mult<Accessor> >());
}



//////////////////
// operator *   //
// SE2 * Matrix //
//////////////////

template <int RHS, class Accessor>
struct SE2MMult {
  inline static void eval(Matrix<4,RHS>& ret, const SE2& lhs, const FixedMatrix<4,RHS,Accessor>& rhs){
    for(int i=0; i<RHS; i++){
      ret.T()[i].template slice<0,2>() = lhs.get_rotation() * rhs.T()[i].template slice<0,2>();
      ret.T()[i].template slice<0,2>() += lhs.get_translation() * rhs(2,i);
      ret(2,i) = rhs(2,i);
    }
  }
};


template <int RHS, class Accessor> inline 
Matrix<4,RHS> operator*(const SE2& lhs, const FixedMatrix<4,RHS,Accessor>& rhs){
  return Matrix<4,RHS>(lhs,rhs,Operator<SE2MMult<RHS,Accessor> >());
}


//////////////////
// operator *   //
// Matrix * SE2 //
//////////////////

template <int LHS, class Accessor> 
struct MSE2Mult {
  inline static void eval(Matrix<LHS,4>& ret, const FixedMatrix<LHS,4,Accessor>& lhs, const SE2& rhs){
    for(int i=0; i<LHS; i++){
      ret[i].template slice<0,2>() = lhs[i].template slice<0,2>() * rhs.get_rotation();
      ret(i,2) = rhs.get_translation() * lhs[i].template slice<0,2>();
      ret(i,2) += lhs(i,2);
    }
  }
};


template <int LHS, class Accessor> inline 
Matrix<LHS,4> operator*(const FixedMatrix<LHS,4,Accessor>& lhs, const SE2& rhs){
  return Matrix<LHS,4>(lhs,rhs,Operator<MSE2Mult<LHS,Accessor> >());
}


namespace SE2static 
{
  static double zero[2]={0,0};
}

inline SE2::SE2() :
  my_translation(SE2static::zero)
{}


/* inline SE2 SE2::exp(const Vector<3>& mu){ */
/*   SE2 result; */
/*   double theta = mu[2]; */
/*   result.get_rotation() = SO2::exp(theta); */
/*   Matrix<2> m2; */
/*   m2[0][0] = m2[1][1] = result.get_rotation().get_matrix()[1][0]; */
/*   m2[0][1] = result.get_rotation().get_matrix()[0][0] - 1.0; */
/*   m2[1][0] = - m2[0][1]; */
/*   if(theta != 0.0)  */
/*     result.get_translation() = m2 * mu.slice<0,2>() / fabs(theta); */
/*   else */
/*     result.get_translation() = mu.slice<0,2>(); */
/*   return result; */
/* } */

inline SE2 SE2::exp(const Vector<3>& mu)
{
  static const double one_6th = 1.0/6.0;
  static const double one_20th = 1.0/20.0;
  
  SE2 result;
  
  const double theta = mu[2];
  const double theta_sq = theta * theta;
  
  Vector<2> cross;
  cross[0] = -theta * mu[1];
  cross[1] = theta * mu[0];
  
  result.my_rotation = SO2::exp(theta);
  
  if (theta_sq < 1e-8) 
    {
      result.get_translation() = mu.slice<0,2>() + 0.5 * cross;
    } 
  else 
    {
      double A, B, C;
      if (theta_sq < 1e-6) 
	{
	  C = one_6th*(1.0 - one_20th * theta_sq);
	  A = 1.0 - theta_sq * C;
	  B = 0.5 - 0.25 * one_6th * theta_sq;
	} 
      else
	{
	  const double inv_theta = (1.0/theta);
	  double sine = result.my_rotation.get_matrix()[1][0];
	  double cosine = result.my_rotation.get_matrix()[0][0];
	  A = sine * inv_theta;
	  B = (1 - cosine) * (inv_theta * inv_theta);
	  C = (1 - A) * (inv_theta * inv_theta);
	}
      result.get_translation() = (1.0 - C * theta_sq) * mu.slice<0,2>() + B * cross;
    }
  return result;
}
 

inline Vector<3> SE2::ln(const SE2& se2) {
  double theta = se2.my_rotation.ln();
  double shtot = 0.5;
  
  if(fabs(theta) > 0.00001) {
    shtot = sin(theta/2)/theta;
  }
    
  // now do the rotation
  SO2 halfrotator = SO2::exp(-0.5 * theta );
  Vector<2> rottrans = halfrotator * se2.my_translation;
  rottrans /= (2 * shtot);
  
  Vector<3> result;
  result.slice<0,2>()=rottrans;
  result[2] = theta;
  return result;
}

inline SE2 SE2::inverse() const {
  const SO2& rinv = my_rotation.inverse();
    return SE2(rinv, -(rinv*my_translation));
}

inline SE2& SE2::left_multiply_by(const SE2& left) {
    my_translation = left.my_translation + left.get_rotation() * my_translation;
    my_rotation = left.my_rotation * my_rotation;
    return *this;
}

inline Vector<3> SE2::generator_field(int i, Vector<3> pos){
  double result_d[]={0,0,0};
  Vector<3> result(result_d);
  if(i < 2){
    result[i]=pos[2];
    return result;
  }
  result[0] = - pos[1];
  result[1] =   pos[0];
  return result;
}


inline SE2 operator*(const SO2& lhs, const SE2& rhs){
  SE2 result;
  result.my_rotation = lhs*rhs.my_rotation;
  result.my_translation = lhs*rhs.my_translation;
  return result;
}

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
