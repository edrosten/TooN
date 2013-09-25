// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk), Gerhard Reitmayr (gr281@cam.ac.uk)

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

/* This code mostly made by copying from se3.h !! */

#ifndef TOON_INCLUDE_SE2_H
#define TOON_INCLUDE_SE2_H

#include <TooN/so2.h>


namespace TooN {

/// Represent a two-dimensional Euclidean transformation (a rotation and a translation). 
/// This can be represented by a \f$2\times 3\f$ matrix operating on a homogeneous co-ordinate, 
/// so that a vector \f$\underline{x}\f$ is transformed to a new location \f$\underline{x}'\f$
/// by
/// \f[\begin{aligned}\underline{x}' &= E\times\underline{x}\\ \begin{bmatrix}x'\\y'\end{bmatrix} &= \begin{pmatrix}r_{11} & r_{12} & t_1\\r_{21} & r_{22} & t_2\end{pmatrix}\begin{bmatrix}x\\y\\1\end{bmatrix}\end{aligned}\f]
/// 
/// This transformation is a member of the Special Euclidean Lie group SE2. These can be parameterised with
/// three numbers (in the space of the Lie Algebra). In this class, the first two parameters are a
/// translation vector while the third is the amount of rotation in the plane as for SO2.
/// @ingroup gTransforms
template <typename Precision = DefaultPrecision>
class SE2 {
public:
	/// Default constructor. Initialises the the rotation to zero (the identity) and the translation to zero
	SE2() : my_translation(Zeros) {}
	template <class A> SE2(const SO2<Precision>& R, const Vector<2,Precision,A>& T) : my_rotation(R), my_translation(T) {}
	template <int S, class P, class A> SE2(const Vector<S, P, A> & v) { *this = exp(v); }

	/// Returns the rotation part of the transformation as a SO2
	SO2<Precision> & get_rotation(){return my_rotation;}
	/// @overload
	const SO2<Precision> & get_rotation() const {return my_rotation;}
	/// Returns the translation part of the transformation as a Vector
	Vector<2, Precision> & get_translation() {return my_translation;}
	/// @overload
	const Vector<2, Precision> & get_translation() const {return my_translation;}

	/// Exponentiate a Vector in the Lie Algebra to generate a new SE2.
	/// See the Detailed Description for details of this vector.
	/// @param vect The Vector to exponentiate
	template <int S, typename P, typename A>
	static inline SE2 exp(const Vector<S,P, A>& vect);
	
	/// Take the logarithm of the matrix, generating the corresponding vector in the Lie Algebra.
	/// See the Detailed Description for details of this vector.
	static inline Vector<3, Precision> ln(const SE2& se2);
	/// @overload
	Vector<3, Precision> ln() const { return SE2::ln(*this); }

	/// compute the inverse of the transformation
	SE2 inverse() const {
		const SO2<Precision> & rinv = my_rotation.inverse();
		return SE2(rinv, -(rinv*my_translation));
	};

	/// Right-multiply by another SE2 (concatenate the two transformations)
	/// @param rhs The multipier
	template <typename P>
	SE2<typename Internal::MultiplyType<Precision,P>::type> operator *(const SE2<P>& rhs) const { 
		return SE2<typename Internal::MultiplyType<Precision,P>::type>(my_rotation*rhs.get_rotation(), my_translation + my_rotation*rhs.get_translation()); 
	}

	/// Self right-multiply by another SE2 (concatenate the two transformations)
	/// @param rhs The multipier
	template <typename P>
	inline SE2& operator *=(const SE2<P>& rhs) { 
		*this = *this * rhs; 
		return *this; 
	}

	/// returns the generators for the Lie group. These are a set of matrices that
	/// form a basis for the vector space of the Lie algebra.
	/// - 0 is translation in x
	/// - 1 is translation in y
	/// - 2 is rotation in the plane
	static inline Matrix<3,3, Precision> generator(int i) {
		Matrix<3,3,Precision> result(Zeros);
		if(i < 2){
			result[i][2] = 1;
			return result;
		}
		result[0][1] = -1;
		result[1][0] = 1;
		return result;
	}

	/// transfers a vector in the Lie algebra, from one coord frame to another
	/// so that exp(adjoint(vect)) = (*this) * exp(vect) * (this->inverse())
	template<typename Accessor>
	Vector<3, Precision> adjoint(const Vector<3,Precision, Accessor> & vect) const {
		Vector<3, Precision> result;
		result[2] = vect[2];
		result.template slice<0,2>() = my_rotation * vect.template slice<0,2>();
		result[0] += vect[2] * my_translation[1];
		result[1] -= vect[2] * my_translation[0];
		return result;
	}

	template <typename Accessor>
	Matrix<3,3,Precision> adjoint(const Matrix<3,3,Precision,Accessor>& M) const {
		Matrix<3,3,Precision> result;
		for(int i=0; i<3; ++i)
			result.T()[i] = adjoint(M.T()[i]);
		for(int i=0; i<3; ++i)
			result[i] = adjoint(result[i]);
		return result;
	}

private:
	SO2<Precision> my_rotation;
	Vector<2, Precision> my_translation;
};

/// Write an SE2 to a stream 
/// @relates SE2
template <class Precision>
inline std::ostream& operator<<(std::ostream& os, const SE2<Precision> & rhs){
	std::streamsize fw = os.width();
	for(int i=0; i<2; i++){
		os.width(fw);
		os << rhs.get_rotation().get_matrix()[i];
		os.width(fw);
		os << rhs.get_translation()[i] << '\n';
	}
	return os;
}

/// Read an SE2 from a stream 
/// @relates SE2
template <class Precision>
inline std::istream& operator>>(std::istream& is, SE2<Precision>& rhs){
	for(int i=0; i<2; i++)
		is >> rhs.get_rotation().my_matrix[i].ref() >> rhs.get_translation()[i];
	rhs.get_rotation().coerce();
	return is;
}


//////////////////
// operator *   //
// SE2 * Vector //
//////////////////

namespace Internal {
template<int S, typename P, typename PV, typename A>
struct SE2VMult;
}

template<int S, typename P, typename PV, typename A>
struct Operator<Internal::SE2VMult<S,P,PV,A> > {
	const SE2<P> & lhs;
	const Vector<S,PV,A> & rhs;
	
	Operator(const SE2<P> & l, const Vector<S,PV,A> & r ) : lhs(l), rhs(r) {}
	
	template <int S0, typename P0, typename A0>
	void eval(Vector<S0, P0, A0> & res ) const {
		SizeMismatch<3,S>::test(3, rhs.size());
		res.template slice<0,2>()=lhs.get_rotation()*rhs.template slice<0,2>();
		res.template slice<0,2>()+=lhs.get_translation() * rhs[2];
		res[2] = rhs[2];
	}
	int size() const { return 3; }
};

/// Right-multiply with a Vector<3>
/// @relates SE2
template<int S, typename P, typename PV, typename A>
inline Vector<3, typename Internal::MultiplyType<P,PV>::type> operator*(const SE2<P> & lhs, const Vector<S,PV,A>& rhs){
	return Vector<3, typename Internal::MultiplyType<P,PV>::type>(Operator<Internal::SE2VMult<S,P,PV,A> >(lhs,rhs));
}

/// Right-multiply with a Vector<2> (special case, extended to be a homogeneous vector)
/// @relates SE2
template <typename P, typename PV, typename A>
inline Vector<2, typename Internal::MultiplyType<P,PV>::type> operator*(const SE2<P>& lhs, const Vector<2,PV,A>& rhs){
	return lhs.get_translation() + lhs.get_rotation() * rhs;
}

//////////////////
// operator *   //
// Vector * SE2 //
//////////////////

namespace Internal {
template<int S, typename P, typename PV, typename A>
struct VSE2Mult;
}

template<int S, typename P, typename PV, typename A>
struct Operator<Internal::VSE2Mult<S,P,PV,A> > {
	const Vector<S,PV,A> & lhs;
	const SE2<P> & rhs;
	
	Operator(const Vector<S,PV,A> & l, const SE2<P> & r ) : lhs(l), rhs(r) {}
	
	template <int S0, typename P0, typename A0>
	void eval(Vector<S0, P0, A0> & res ) const {
		SizeMismatch<3,S>::test(3, lhs.size());
		res.template slice<0,2>() = lhs.template slice<0,2>()*rhs.get_rotation();
		res[2] = lhs[2];
		res[2] += lhs.template slice<0,2>() * rhs.get_translation();
	}
	int size() const { return 3; }
};

/// Left-multiply with a Vector<3>
/// @relates SE2
template<int S, typename P, typename PV, typename A>
inline Vector<3, typename Internal::MultiplyType<PV,P>::type> operator*(const Vector<S,PV,A>& lhs, const SE2<P> & rhs){
	return Vector<3, typename Internal::MultiplyType<PV,P>::type>(Operator<Internal::VSE2Mult<S, P,PV,A> >(lhs,rhs));
}

//////////////////
// operator *   //
// SE2 * Matrix //
//////////////////

namespace Internal {
template <int R, int C, typename PM, typename A, typename P>
struct SE2MMult;
}

template<int R, int Cols, typename PM, typename A, typename P>
struct Operator<Internal::SE2MMult<R, Cols, PM, A, P> > {
	const SE2<P> & lhs;
	const Matrix<R,Cols,PM,A> & rhs;
	
	Operator(const SE2<P> & l, const Matrix<R,Cols,PM,A> & r ) : lhs(l), rhs(r) {}
	
	template <int R0, int C0, typename P0, typename A0>
	void eval(Matrix<R0, C0, P0, A0> & res ) const {
		SizeMismatch<3,R>::test(3, rhs.num_rows());
		for(int i=0; i<rhs.num_cols(); ++i)
			res.T()[i] = lhs * rhs.T()[i];
	}
	int num_cols() const { return rhs.num_cols(); }
	int num_rows() const { return 3; }
};

/// Right-multiply with a Matrix<3>
/// @relates SE2
template <int R, int Cols, typename PM, typename A, typename P> 
inline Matrix<3,Cols, typename Internal::MultiplyType<P,PM>::type> operator*(const SE2<P> & lhs, const Matrix<R,Cols,PM, A>& rhs){
	return Matrix<3,Cols,typename Internal::MultiplyType<P,PM>::type>(Operator<Internal::SE2MMult<R, Cols, PM, A, P> >(lhs,rhs));
}

//////////////////
// operator *   //
// Matrix * SE2 //
//////////////////

namespace Internal {
template <int Rows, int C, typename PM, typename A, typename P>
struct MSE2Mult;
}

template<int Rows, int C, typename PM, typename A, typename P>
struct Operator<Internal::MSE2Mult<Rows, C, PM, A, P> > {
	const Matrix<Rows,C,PM,A> & lhs;
	const SE2<P> & rhs;
	
	Operator( const Matrix<Rows,C,PM,A> & l, const SE2<P> & r ) : lhs(l), rhs(r) {}
	
	template <int R0, int C0, typename P0, typename A0>
	void eval(Matrix<R0, C0, P0, A0> & res ) const {
		SizeMismatch<3, C>::test(3, lhs.num_cols());
		for(int i=0; i<lhs.num_rows(); ++i)
			res[i] = lhs[i] * rhs;
	}
	int num_cols() const { return 3; }
	int num_rows() const { return lhs.num_rows(); }
};

/// Left-multiply with a Matrix<3>
/// @relates SE2
template <int Rows, int C, typename PM, typename A, typename P> 
inline Matrix<Rows,3, typename Internal::MultiplyType<PM,P>::type> operator*(const Matrix<Rows,C,PM, A>& lhs, const SE2<P> & rhs ){
	return Matrix<Rows,3,typename Internal::MultiplyType<PM,P>::type>(Operator<Internal::MSE2Mult<Rows, C, PM, A, P> >(lhs,rhs));
}

template <typename Precision>
template <int S, typename PV, typename Accessor>
inline SE2<Precision> SE2<Precision>::exp(const Vector<S, PV, Accessor>& mu)
{
	SizeMismatch<3,S>::test(3, mu.size());

	static const Precision one_6th = 1.0/6.0;
	static const Precision one_20th = 1.0/20.0;
  
	SE2<Precision> result;
  
	const Precision theta = mu[2];
	const Precision theta_sq = theta * theta;
  
	const Vector<2, Precision> cross = makeVector( -theta * mu[1], theta * mu[0]);
	result.get_rotation() = SO2<Precision>::exp(theta);

	if (theta_sq < 1e-8){
		result.get_translation() = mu.template slice<0,2>() + 0.5 * cross;
	} else {
		Precision A, B;
		if (theta_sq < 1e-6) {
			A = 1.0 - theta_sq * one_6th*(1.0 - one_20th * theta_sq);
			B = 0.5 - 0.25 * one_6th * theta_sq;
		} else {
			const Precision inv_theta = (1.0/theta);
			const Precision sine = result.my_rotation.get_matrix()[1][0];
			const Precision cosine = result.my_rotation.get_matrix()[0][0];
			A = sine * inv_theta;
			B = (1 - cosine) * (inv_theta * inv_theta);
		}
		result.get_translation() = TooN::operator*(A,mu.template slice<0,2>()) + TooN::operator*(B,cross);
	}
	return result;
}
 
template <typename Precision>
inline Vector<3, Precision> SE2<Precision>::ln(const SE2<Precision> & se2) {
	const Precision theta = se2.get_rotation().ln();

	Precision shtot = 0.5;  
	if(fabs(theta) > 0.00001)
		shtot = sin(theta/2)/theta;

	const SO2<Precision> halfrotator(theta * -0.5);
	Vector<3, Precision> result;
	result.template slice<0,2>() = (halfrotator * se2.get_translation())/(2 * shtot);
	result[2] = theta;
	return result;
}

/// Multiply a SO2 with and SE2
/// @relates SE2
/// @relates SO2
template <typename Precision>
inline SE2<Precision> operator*(const SO2<Precision> & lhs, const SE2<Precision>& rhs){
	return SE2<Precision>( lhs*rhs.get_rotation(), lhs*rhs.get_translation());
}

}
#endif
