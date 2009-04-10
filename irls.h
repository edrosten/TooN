// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk)
//
// This file is part of the TooN Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.


#ifndef __IRLS_H
#define __IRLS_H

#include <TooN/wls.h>
#include <cassert>
#include <cmath>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

/// Robust reweighting (type I) for IRLS.
/// A reweighting class with \f$w(x)=\frac{1}{\sigma + |x|}\f$.
/// This structure can be passed as the second template argument in IRLS.
/// @ingroup gEquations
struct RobustI {
  double sd_inlier; ///< The inlier standard deviation, \f$\sigma\f$.
  inline double reweight(double x) {return 1/(sd_inlier+fabs(x));}  ///< Returns \f$w(x)\f$.
  inline double true_scale(double x) {return reweight(x) - fabs(x)*reweight(x)*reweight(x);}  ///< Returns \f$w(x) + xw'(x)\f$.
  inline double objective(double x) {return fabs(x) + sd_inlier*log(sd_inlier*reweight(x));}  ///< Returns \f$\int xw(x)dx\f$.
};

/// Robust reweighting (type II) for IRLS.
/// A reweighting class with \f$w(x)=\frac{1}{\sigma + x^2}\f$.
/// This structure can be passed as the second template argument in IRLS.
/// @ingroup gEquations
struct RobustII {
  double sd_inlier; ///< The inlier standard deviation, \f$\sigma\f$.
  inline double reweight(double d){return 1/(sd_inlier+d*d);} ///< Returns \f$w(x)\f$.
  inline double true_scale(double d){return d - 2*d*reweight(d);} ///< Returns \f$w(x) + xw'(x)\f$.
  inline double objective(double d){return 0.5 * log(1 + d*d/sd_inlier);} ///< Returns \f$\int xw(x)dx\f$.
};

/// A reweighting class representing no reweighting in IRLS.
/// \f$w(x)=1\f$
/// This structure can be passed as the second template argument in IRLS.
/// @ingroup gEquations
struct ILinear {
  inline double reweight(double d){return 1;} ///< Returns \f$w(x)\f$.
  inline double true_scale(double d){return 1;} ///< Returns \f$w(x) + xw'(x)\f$.
  inline double objective(double d){return d*d;} ///< Returns \f$\int xw(x)dx\f$.
};


/// Performs iterative reweighted least squares.
/// @param Size the size
/// @param Reweight The reweighting functor. This structure must provide reweight(), 
/// true-scale() and objective() methods. Existing examples are  Robust I, Robust II and ILinear.
/// @ingroup gEquations
template <int Size, class Reweight>
class IRLS
  : public Reweight,
    public WLS<Size>
{
public:
  IRLS(){Identity(my_true_C_inv,0);my_residual=0;}

  inline void add_df(double d, const Vector<Size>& f) {
    double scale = Reweight::reweight(d);
    double ts = Reweight::true_scale(d);
    my_residual += Reweight::objective(d);

    WLS<Size>::add_df(d,f,scale);

    for(int i=0; i<Size; i++){
      for(int j=0; j<Size; j++){
	my_true_C_inv[i][j]+=f[i]*f[j]*ts;
      }
    }
  }

  void operator += (const IRLS& meas){
    WLS<Size>::operator+=(meas);
    my_true_C_inv += meas.my_true_C_inv;
  }


  Matrix<Size,Size,RowMajor>& get_true_C_inv() {return my_true_C_inv;}
  const Matrix<Size,Size,RowMajor>& get_true_C_inv()const {return my_true_C_inv;}

  double get_residual() {return my_residual;}

private:

  double my_residual;

  Matrix<Size,Size,RowMajor> my_true_C_inv;

  // comment out to allow bitwise copying
  IRLS( IRLS& copyof );
  int operator = ( IRLS& copyof );
};

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
