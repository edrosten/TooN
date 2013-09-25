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

#ifndef TOON_INCLUDE_SO2_H
#define TOON_INCLUDE_SO2_H

#include <TooN/TooN.h>
#include <TooN/helpers.h>

namespace TooN {

template <typename Precision> class SO2;
template <typename Precision> class SE2;
template <typename Precision> class SIM2;

template<typename Precision> inline std::istream & operator>>(std::istream &, SO2<Precision> & );
template<typename Precision> inline std::istream & operator>>(std::istream &, SE2<Precision> & );
template<typename Precision> inline std::istream & operator>>(std::istream &, SIM2<Precision> & );

/// Class to represent a two-dimensional rotation matrix. Two-dimensional rotation
/// matrices are members of the Special Orthogonal Lie group SO2. This group can be parameterised with
/// one number (the rotation angle).
/// @ingroup gTransforms
template<typename Precision = DefaultPrecision>
class SO2 {
	friend std::istream& operator>> <Precision>(std::istream&, SO2& );
	friend std::istream& operator>> <Precision>(std::istream&, SE2<Precision>& );
	friend std::istream& operator>> <Precision>(std::istream&, SIM2<Precision>& );

public:
	/// Default constructor. Initialises the matrix to the identity (no rotation)
	SO2() : my_matrix(Identity) {} 
	
	/// Construct from a rotation matrix.
	SO2(const Matrix<2,2,Precision>& rhs) {  
		*this = rhs; 
		coerce();
	}

	/// Construct from an angle.
	SO2(const Precision l) { *this = exp(l); }
  
	/// Assigment operator from a general matrix. This also calls coerce()
	/// to make sure that the matrix is a valid rotation matrix.
	template <int R, int C, typename P, typename A> 
	SO2& operator=(const Matrix<R,C,P,A>& rhs){
		my_matrix = rhs;
		coerce();
		return *this;
	}

	/// Modifies the matrix to make sure it is a valid rotation matrix.
	void coerce(){
		my_matrix[0] = unit(my_matrix[0]);
		my_matrix[1] -= my_matrix[0] * (my_matrix[0]*my_matrix[1]);
		my_matrix[1] = unit(my_matrix[1]);
	}

	/// Exponentiate an angle in the Lie algebra to generate a new SO2.
	inline static SO2 exp(const Precision & d){
		SO2<Precision> result;
		result.my_matrix[0][0] = result.my_matrix[1][1] = cos(d);
		result.my_matrix[1][0] = sin(d);
		result.my_matrix[0][1] = -result.my_matrix[1][0];
		return result;
	}

	/// extracts the rotation angle from the SO2
	Precision ln() const { return atan2(my_matrix[1][0], my_matrix[0][0]); }

	/// Returns the inverse of this matrix (=the transpose, so this is a fast operation)
	SO2 inverse() const { return SO2(*this, Invert()); }

	/// Self right-multiply by another rotation matrix
	template <typename P>
	SO2& operator *=(const SO2<P>& rhs){
		my_matrix=my_matrix*rhs.get_matrix();
		return *this;
	}

	/// Right-multiply by another rotation matrix
	template <typename P>
	SO2<typename Internal::MultiplyType<Precision, P>::type> operator *(const SO2<P>& rhs) const { 
		return SO2<typename Internal::MultiplyType<Precision, P>::type>(*this,rhs); 
	}

	/// Returns the SO2 as a Matrix<2>
	const Matrix<2,2,Precision>& get_matrix() const {return my_matrix;}

	/// returns generator matrix
	static Matrix<2,2,Precision> generator() {
		Matrix<2,2,Precision> result;
		result[0] = makeVector(0,-1);
		result[1] = makeVector(1,0);
		return result;
	}

private:
	struct Invert {};
	inline SO2(const SO2& so2, const Invert&) : my_matrix(so2.my_matrix.T()) {}
	template <typename PA, typename PB>
	inline SO2(const SO2<PA>& a, const SO2<PB>& b) : my_matrix(a.get_matrix()*b.get_matrix()) {}

	Matrix<2,2,Precision> my_matrix;
};

/// Write an SO2 to a stream 
/// @relates SO2
template <typename Precision>
inline std::ostream& operator<< (std::ostream& os, const SO2<Precision> & rhs){
	return os << rhs.get_matrix();
}

/// Read from SO2 to a stream 
/// @relates SO2
template <typename Precision>
inline std::istream& operator>>(std::istream& is, SO2<Precision>& rhs){
	is >> rhs.my_matrix;
	rhs.coerce();
	return is;
}

/// Right-multiply by a Vector
/// @relates SO2
template<int D, typename P1, typename PV, typename Accessor>
inline Vector<2, typename Internal::MultiplyType<P1, PV>::type> operator*(const SO2<P1> & lhs, const Vector<D, PV, Accessor> & rhs){
	return lhs.get_matrix() * rhs;
}

/// Left-multiply by a Vector
/// @relates SO2
template<int D, typename P1, typename PV, typename Accessor>
inline Vector<2, typename Internal::MultiplyType<PV,P1>::type> operator*(const Vector<D, PV, Accessor>& lhs, const SO2<P1> & rhs){
	return lhs * rhs.get_matrix();
}

/// Right-multiply by a Matrix
/// @relates SO2
template <int R, int C, typename P1, typename P2, typename Accessor> 
inline Matrix<2,C,typename Internal::MultiplyType<P1,P2>::type> operator*(const SO2<P1> & lhs, const Matrix<R,C,P2,Accessor>& rhs){
	return lhs.get_matrix() * rhs;
}

/// Left-multiply by a Matrix
/// @relates SO2
template <int R, int C, typename P1, typename P2, typename Accessor>
inline Matrix<R,2,typename Internal::MultiplyType<P1,P2>::type> operator*(const Matrix<R,C,P1,Accessor>& lhs, const SO2<P2>& rhs){
	return lhs * rhs.get_matrix();
}

}

#endif
