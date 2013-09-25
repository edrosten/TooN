// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Gerhard Reitmayr (gr281@cam.ac.uk)

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

#ifndef TOON_INCLUDE_SL_H
#define TOON_INCLUDE_SL_H

#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <TooN/gaussian_elimination.h>
#include <TooN/determinant.h>
#include <cassert>

namespace TooN {

template <int N, typename P> class SL;
template <int N, typename P> std::istream & operator>>(std::istream &, SL<N, P> &);

/// represents an element from the group SL(n), the NxN matrices M with det(M) = 1.
/// This can be used to conveniently estimate homographies on n-1 dimentional spaces.
/// The implementation uses the matrix exponential function @ref exp for
/// exponentiation from an element in the Lie algebra and LU to compute an inverse.
/// 
/// The Lie algebra are the NxN matrices M with trace(M) = 0. The N*N-1 generators used
/// to represent this vector space are the following:
/// - diag(...,1,-1,...), n-1 along the diagonal
/// - symmetric generators for every pair of off-diagonal elements
/// - anti-symmetric generators for every pair of off-diagonal elements
/// This choice represents the fact that SL(n) can be interpreted as the product
/// of all symmetric matrices with det() = 1 times SO(n).
/// @ingroup gTransforms
template <int N, typename Precision = DefaultPrecision>
class SL {
	friend std::istream & operator>> <N,Precision>(std::istream &, SL &);
public:
	static const int size = N;          ///< size of the matrices represented by SL<N>
	static const int dim = N*N - 1;     ///< dimension of the vector space represented by SL<N>

	/// default constructor, creates identity element
	SL() : my_matrix(Identity) {}

	/// exp constructor, creates element through exponentiation of Lie algebra vector. see @ref SL::exp.
	template <int S, typename P, typename B>
	SL( const Vector<S,P,B> & v ) { *this = exp(v); }

	/// copy constructor from a matrix, coerces matrix to be of determinant = 1
	template <int R, int C, typename P, typename A>
	SL(const Matrix<R,C,P,A>& M) : my_matrix(M) { coerce(); }

	/// returns the represented matrix
	const Matrix<N,N,Precision> & get_matrix() const { return my_matrix; }
	/// returns the inverse using LU
	SL inverse() const { return SL(*this, Invert()); }

	/// multiplies to SLs together by multiplying the underlying matrices
	template <typename P>
	SL<N,typename Internal::MultiplyType<Precision, P>::type> operator*( const SL<N,P> & rhs) const { return SL<N,typename Internal::MultiplyType<Precision, P>::type>(*this, rhs); }

	/// right multiplies this SL with another one
	template <typename P>
	SL operator*=( const SL<N,P> & rhs) { *this = *this*rhs; return *this; }

	/// exponentiates a vector in the Lie algebra to compute the corresponding element
	/// @arg v a vector of dimension SL::dim
	template <int S, typename P, typename B>
	static inline SL exp( const Vector<S,P,B> &);

	inline Vector<N*N-1, Precision> ln() const ;

	/// returns one generator of the group. see SL for a detailed description of 
	/// the generators used.
	/// @arg i number of the generator between 0 and SL::dim -1 inclusive
	static inline Matrix<N,N,Precision> generator(int);

private:
	struct Invert {};
	SL( const SL & from, struct Invert ) {
		const Matrix<N> id = Identity;
		my_matrix = gaussian_elimination(from.my_matrix, id);
	}
	SL( const SL & a, const SL & b) : my_matrix(a.get_matrix() * b.get_matrix()) {}

	void coerce(){
		using std::abs;
		Precision det = determinant(my_matrix);
		assert(abs(det) > 0);
        using std::pow;
		my_matrix /= pow(det, 1.0/N);
	}

	/// these constants indicate which parts of the parameter vector 
	/// map to which generators
	///{
	static const int COUNT_DIAG = N - 1;
	static const int COUNT_SYMM = (dim - COUNT_DIAG)/2;
	static const int COUNT_ASYMM = COUNT_SYMM;
	static const int DIAG_LIMIT = COUNT_DIAG;
	static const int SYMM_LIMIT = COUNT_SYMM + DIAG_LIMIT;
	///}

	Matrix<N,N,Precision> my_matrix;
};

template <int N, typename Precision>
template <int S, typename P, typename B>
inline SL<N, Precision> SL<N, Precision>::exp( const Vector<S,P,B> & v){
	SizeMismatch<S,dim>::test(v.size(), dim);
	Matrix<N,N,Precision> t(Zeros);
	for(int i = 0; i < dim; ++i)
		t += generator(i) * v[i];
	SL<N, Precision> result;
	result.my_matrix = TooN::exp(t);
	return result;
}

template <int N, typename Precision>
inline Vector<N*N-1, Precision> SL<N, Precision>::ln() const {
	const Matrix<N> l = TooN::log(my_matrix);
	Vector<SL<N,Precision>::dim, Precision> v;
	Precision last = 0;
	for(int i = 0; i < DIAG_LIMIT; ++i){	// diagonal elements
		v[i] = l(i,i) + last;
		last = l(i,i);
	}
	for(int i = DIAG_LIMIT, row = 0, col = 1; i < SYMM_LIMIT; ++i) {	// calculate symmetric and antisymmetric in one go
		// do the right thing here to calculate the correct indices !
		v[i] = (l(row, col) + l(col, row))*0.5;
		v[i+COUNT_SYMM] = (-l(row, col) + l(col, row))*0.5;
		++col;
		if( col == N ){
			++row;
			col = row+1;
		}
	}
	return v;
}

template <int N, typename Precision>
inline Matrix<N,N,Precision> SL<N, Precision>::generator(int i){
	assert( i > -1 && i < dim );
	Matrix<N,N,Precision> result(Zeros);
	if(i < DIAG_LIMIT) { 				// first ones are the diagonal ones
		result(i,i) = 1;
		result(i+1,i+1) = -1;
	} else if(i < SYMM_LIMIT){			// then the symmetric ones
		int row = 0, col = i - DIAG_LIMIT + 1;
		while(col > (N - row - 1)){
			col -= (N - row - 1); 
			++row;
		}
		col += row;
		result(row, col) = result(col, row) = 1;
	} else {							// finally the antisymmetric ones
		int row = 0, col = i - SYMM_LIMIT + 1;
		while(col > N - row - 1){
			col -= N - row - 1; 
			++row;
		}
		col += row;
		result(row, col) = -1;
		result(col, row) = 1;
	}
	return result;
}

template <int S, typename PV, typename B, int N, typename P>
Vector<N, typename Internal::MultiplyType<P, PV>::type> operator*( const SL<N, P> & lhs, const Vector<S,PV,B> & rhs ){
	return lhs.get_matrix() * rhs;
}

template <int S, typename PV, typename B, int N, typename P>
Vector<N, typename Internal::MultiplyType<PV, P>::type> operator*( const Vector<S,PV,B> & lhs, const SL<N,P> & rhs ){
	return lhs * rhs.get_matrix();
}

template<int R, int C, typename PM, typename A, int N, typename P> inline
Matrix<N, C, typename Internal::MultiplyType<P, PM>::type> operator*(const SL<N,P>& lhs, const Matrix<R, C, PM, A>& rhs){
	return lhs.get_matrix() * rhs;
}

template<int R, int C, typename PM, typename A, int N, typename P> inline
Matrix<R, N, typename Internal::MultiplyType<PM, P>::type> operator*(const Matrix<R, C, PM, A>& lhs, const SL<N,P>& rhs){
	return lhs * rhs.get_matrix();
}

template <int N, typename P>
std::ostream & operator<<(std::ostream & out, const SL<N, P> & h){
	out << h.get_matrix();
	return out;
}

template <int N, typename P>
std::istream & operator>>(std::istream & in, SL<N, P> & h){
	in >> h.my_matrix;
	h.coerce();
	return in;
}

};

#endif
