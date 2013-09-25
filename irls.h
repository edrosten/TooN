// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk)

//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND OTHER CONTRIBUTORS ``AS IS''
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR OTHER CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.


#ifndef __IRLS_H
#define __IRLS_H

#include <TooN/wls.h>
#include <cassert>
#include <cmath>

namespace TooN {

	/// Robust reweighting (type I) for IRLS.
	/// A reweighting class with \f$w(x)=\frac{1}{\sigma + |x|}\f$.
	/// This structure can be passed as the second template argument in IRLS.
	/// @ingroup gEquations
	template<typename Precision>
	struct RobustI {
		void set_sd(Precision x){ sd_inlier = x;} ///<Set the noise standard deviation.
		double sd_inlier; ///< The inlier standard deviation, \f$\sigma\f$.
		inline Precision reweight(Precision x) {return 1/(sd_inlier+fabs(x));}  ///< Returns \f$w(x)\f$.
		inline Precision true_scale(Precision x) {return reweight(x) - fabs(x)*reweight(x)*reweight(x);}  ///< Returns \f$w(x) + xw'(x)\f$.
		inline Precision objective(Precision x) {return fabs(x) + sd_inlier*::log(sd_inlier*reweight(x));}  ///< Returns \f$\int xw(x)dx\f$.
	};

	/// Robust reweighting (type II) for IRLS.
	/// A reweighting class with \f$w(x)=\frac{1}{\sigma + x^2}\f$.
	/// This structure can be passed as the second template argument in IRLS.
	/// @ingroup gEquations
	template<typename Precision>
	struct RobustII {
		void set_sd(Precision x){ sd_inlier = x*x;} ///<Set the noise standard deviation.
		Precision sd_inlier; ///< The inlier standard deviation squared, \f$\sigma\f$.
		inline Precision reweight(Precision d){return 1/(sd_inlier+d*d);} ///< Returns \f$w(x)\f$.
		inline Precision true_scale(Precision d){return d - 2*d*reweight(d);} ///< Returns \f$w(x) + xw'(x)\f$.
		inline Precision objective(Precision d){return 0.5 * ::log(1 + d*d/sd_inlier);} ///< Returns \f$\int xw(x)dx\f$.
	};

	/// A reweighting class representing no reweighting in IRLS.
	/// \f$w(x)=1\f$
	/// This structure can be passed as the second template argument in IRLS.
	/// @ingroup gEquations
	template<typename Precision>
	struct ILinear {
		void set_sd(Precision){} ///<Set the noise standard deviation (does nothing).
		inline Precision reweight(Precision d){return 1;} ///< Returns \f$w(x)\f$.
		inline Precision true_scale(Precision d){return 1;} ///< Returns \f$w(x) + xw'(x)\f$.
		inline Precision objective(Precision d){return d*d;} ///< Returns \f$\int xw(x)dx\f$.
	};
	
	///A reweighting class where the objective function tends to a 
	///fixed value, rather than infinity. Note that this is not therefore
	///a proper distribution since its integral is not finite. It is considerably
	///more efficient than RobustI and II, since log() is not used.
	/// @ingroup gEquations
	template<typename Precision>
	struct RobustIII {

		void set_sd(Precision x){ sd_inlier = x*x;} ///<Set the noise standard deviation.
		Precision sd_inlier; ///< Inlier standard deviation squared.
		/// Returns \f$w(x)\f$.
		Precision reweight(Precision x) const
		{
			double d = (1 + x*x/sd_inlier);
			return 1/(d*d);
		}	
		///< Returns \f$\int xw(x)dx\f$.
		Precision objective(Precision x) const 
		{
			return x*x / (2*(1 + x*x/sd_inlier));
		}
	};

	/// Performs iterative reweighted least squares.
	/// @param Size the size
	/// @param Reweight The reweighting functor. This structure must provide reweight(), 
	/// true-scale() and objective() methods. Existing examples are  Robust I, Robust II and ILinear.
	/// @ingroup gEquations
	template <int Size, typename Precision, template <typename Precision> class Reweight>
	class IRLS
		: public Reweight<Precision>,
		  public WLS<Size,Precision>
	{
	public:
		IRLS(int size=Size):
			WLS<Size,Precision>(size),
			my_true_C_inv(Zeros(size))
		{
			my_residual=0;
		}
		
		template<int Size2, typename Precision2, typename Base2>
		inline void add_mJ(Precision m, const Vector<Size2,Precision2,Base2>& J) {
			SizeMismatch<Size,Size2>::test(my_true_C_inv.num_rows(), J.size());

			Precision scale = Reweight<Precision>::reweight(m);
			Precision ts = Reweight<Precision>::true_scale(m);
			my_residual += Reweight<Precision>::objective(m);

			WLS<Size>::add_mJ(m,J,scale);

			Vector<Size,Precision> scaledm(m*ts);

			my_true_C_inv += scaledm.as_col() * scaledm.as_row();

		}

		void operator += (const IRLS& meas){
			WLS<Size>::operator+=(meas);
			my_true_C_inv += meas.my_true_C_inv;
		}


		Matrix<Size,Size,Precision>& get_true_C_inv() {return my_true_C_inv;}
		const Matrix<Size,Size,Precision>& get_true_C_inv()const {return my_true_C_inv;}

		Precision get_residual() {return my_residual;}

		void clear(){
			WLS<Size,Precision>::clear();
			my_residual=0;
			my_true_C_inv = Zeros;
		}

	private:

		Precision my_residual;

		Matrix<Size,Size,Precision> my_true_C_inv;

		// comment out to allow bitwise copying
		IRLS( IRLS& copyof );
		int operator = ( IRLS& copyof );
	};

}

#endif
