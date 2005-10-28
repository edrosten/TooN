
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
#ifndef __WLS_H
#define __WLS_H

#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <TooN/helpers.h>

#include <cassert>
#include <cmath>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

/// Performs weighted least squares computation.
/// @param Size The number of dimensions in the system
/// @ingroup gEquations
template <int Size=-1>
class WLS {
public:
	/// Default constructor
  WLS(){clear();}
  /// Construct using a given regularisation prior
  WLS(double prior){clear(prior);}

  /// Clear all the measurements and apply a constant regularisation term. 
  /// Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
  /// @param prior The strength of the prior
  void clear(double prior=0){
    Identity(my_C_inv,prior);
    for(int i=0; i<Size; i++){
      my_vector[i]=0;
    }
  }

	/// Applies a constant regularisation term. 
	/// Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
	/// @param val The strength of the prior
  void add_prior(double val){
    for(int i=0; i<Size; i++){
      my_C_inv(i,i)+=val;
    }
  }
  
	/// Applies a regularisation term with a different strength for each parameter value. 
	/// Equates to a prior that says all the parameters are zero with \f$\sigma_i^2 = \frac{1}{\text{v}_i}\f$.
	/// @param v The vector of priors
  template<class Accessor>
  void add_prior(const FixedVector<Size,Accessor>& v){
    for(int i=0; i<Size; i++){
      my_C_inv(i,i)+=v[i];
    }
  }

	/// Applies a whole-matrix regularisation term. 
	/// This is the same as adding the \f$m\f$ to the inverse covariance matrix.
	/// @param m The inverse covariance matrix to add
  template<class Accessor>
  void add_prior(const FixedMatrix<Size,Size,Accessor>& m){
    my_C_inv+=m;
  }

	/// Add a single measurement 
	/// @param m The value of the measurement
	/// @param J The Jacobian for the measurement \f$\frac{\partial\text{m}}{\partial\text{param}_i}\f$
	/// @param weight The inverse variance of the measurement (default = 1)
  template<class Accessor>
  inline void add_df(double m, const FixedVector<Size,Accessor>& J, double weight = 1) {
    Vector<Size> Jw = J*weight;
    for(int i=0; i<Size; i++){
      for(int j=0; j<Size; j++){
	my_C_inv[i][j]+=J[i]*Jw[j];
      }
      my_vector[i]+=Jw[i]*m;
    }
  }

	/// Add multiple measurements at once (much more efficiently)
	/// @param N The number of measurements
	/// @param m The measurements to add
	/// @param J The Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
	/// @param invcov The inverse covariance of the measurement values
  template<int N, class Accessor1, class Accessor2, class Accessor3>
  inline void add_df(const FixedVector<N,Accessor1>& m,
		     const FixedMatrix<Size,N,Accessor2>& J,
		     const FixedMatrix<N,N,Accessor3>& invcov){
    my_C_inv += J * invcov * J.T();
    my_vector += J * invcov * m;
  }

  void compute(){
    my_svd.compute(my_C_inv);
    my_mu=my_svd.backsub(my_vector);
  }

  /// Combine measurements from two WLS systems
  /// @param meas The measurements to combine with
  void operator += (const WLS& meas){
    my_vector+=meas.my_vector;
    my_C_inv += meas.my_C_inv;
  }

  /// Returns the inverse covariance matrix
  Matrix<Size,Size,RowMajor>& get_C_inv() {return my_C_inv;}
  /// Returns the inverse covariance matrix
  const Matrix<Size,Size,RowMajor>& get_C_inv() const {return my_C_inv;}
  Vector<Size>& get_mu(){return my_mu;}
  const Vector<Size>& get_mu() const {return my_mu;}
  Vector<Size>& get_vector(){return my_vector;}
  const Vector<Size>& get_vector() const {return my_vector;}
  SVD<Size>& get_svd(){return my_svd;}
  const SVD<Size>& get_svd() const {return my_svd;}

private:
  Vector<Size> my_mu;
  Matrix<Size,Size,RowMajor> my_C_inv;
  Vector<Size> my_vector;
  SVD<Size,Size> my_svd;

  // comment out to allow bitwise copying
  WLS( WLS& copyof );
  int operator = ( WLS& copyof );
};







//  syg21: Dynamic WLS

/// Performs weighted least squares computation.
/// @param Size The number of dimensions in the system
/// @ingroup gMaths
template<> 
class WLS<-1> {
public:
  /// Default constructor 
   //WLS(){clear();} 
  WLS(unsigned int s) 
    : Size(s), my_mu(Size), my_vector(Size), my_C_inv(Size, Size)
    { 
      clear();
    } 
     
    void Identity(Matrix<> &M, double d)
    {
      for(int r=0; r<M.num_rows(); r++)
	{
	  for(int c=0; c<M.num_cols(); c++)
	    {
	      M[r][c] = 0.0;
	    }
	  M[r][r] = 1.0;
	}
    }

  /// Clear all the measurements and apply a constant regularisation term. 
  /// Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
  /// @param prior The strength of the prior
  void clear(double prior=0){
    Identity(my_C_inv,prior);
    for(int i=0; i<Size; i++){
      my_vector[i]=0;
    }
  }

  /// Applies a constant regularisation term. 
  /// Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
  /// @param val The strength of the prior
  void add_prior(double val){
      for(int i=0; i<Size; i++){
      my_C_inv(i,i)+=val;
    }
  }
  
  /// Applies a regularisation term with a different strength for each parameter value. 
  /// Equates to a prior that says all the parameters are zero with \f$\sigma_i^2 = \frac{1}{\text{v}_i}\f$.
  /// @param v The vector of priors
  template<int VSize, class Accessor>
  void add_prior(const FixedVector<VSize,Accessor>& v){
    assert(VSize==Size);
    for(int i=0; i<VSize; i++){
      my_C_inv(i,i)+=v[i];
    }
  }
  /// Applies a whole-matrix regularisation term. 
  /// This is the same as adding the \f$m\f$ to the inverse covariance matrix.
  /// @param m The inverse covariance matrix to add
  template<int MSize, class Accessor>
  void add_prior(const FixedMatrix<MSize,MSize,Accessor>& m){
    my_C_inv+=m;
  }

  /// Add a single measurement 
  /// @param m The value of the measurement
  /// @param J The Jacobian for the measurement \f$\frac{\partial\text{m}}{\partial\text{param}_i}\f$
  /// @param weight The inverse variance of the measurement (default = 1)
  template<int VSize, class Accessor>
  inline void add_df(double m, const FixedVector<VSize,Accessor>& J, double weight = 1) {
    assert(VSize==Size);
    Vector<VSize> Jw = J*weight;
    for(int i=0; i<VSize; i++){
      for(int j=0; j<VSize; j++){
	my_C_inv[i][j]+=J[i]*Jw[j];
      }
      my_vector[i]+=Jw[i]*m;
    }
  }

  /// Add a single measurement 
  /// @param m The value of the measurement
  /// @param J The Jacobian for the measurement \f$\frac{\partial\text{m}}{\partial\text{param}_i}\f$
  /// @param weight The inverse variance of the measurement (default = 1)
  template<class Accessor>
    inline void add_df(double m, const DynamicVector<Accessor>& J, double weight = 1) {
    Vector<> Jw(Size);
    Jw = J * weight;
    for(int i=0; i<Size; i++){
      for(int j=0; j<Size; j++){
	my_C_inv[i][j]+=J[i]*Jw[j];
      }
      my_vector[i]+=Jw[i]*m;
    }
  }


  void compute(){
    //  create SVD
    SVD<> my_svd(my_C_inv);
    //my_svd.compute(my_C_inv);
    my_mu=my_svd.backsub(my_vector);
  }

  /// Returns mu
  const Vector<>& get_mu() const {return my_mu;}
  Vector<>& get_mu() {return my_mu;}
  /// Returns the inverse covariance matrix
  const Matrix<>& get_C_inv() const {return my_C_inv;}
  Matrix<>& get_C_inv() {return my_C_inv;}
  /// Returns my_vector
  const Vector<>& get_vector() const {return my_vector;}
  Vector<>& get_vector(){return my_vector;}


 private: 
  int Size; 
  Vector<> my_mu;
  Vector<> my_vector;
  Matrix<> my_C_inv;

  // comment out to allow bitwise copying
  WLS( WLS& copyof );
  int operator = ( WLS& copyof );


  //  To return SVD, have to create a point to an SVD<> instead of creating a temporary variable in compute()
#if 0
  SVD<Size>& get_svd(){return my_svd;}
  const SVD<Size>& get_svd() const {return my_svd;}
#endif

};













#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
