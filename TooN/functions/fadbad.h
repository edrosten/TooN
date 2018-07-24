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

#ifndef TOON_FADBAD_INTEGATION_H
#define TOON_FADBAD_INTEGATION_H

#include <iostream>

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/se2.h>

#include <FADBAD++/fadiff.h>

namespace fadbad {

template <typename P, unsigned int N>
inline std::ostream & operator<<( std::ostream & out, const F<P, N> & val ){
	return out << val.val();
}

template <typename P, unsigned int N>
inline F<P, N> abs( const F<P, N> & val ){
    return (val.val() < 0) ? -val : val;
}

}

namespace TooN {

template<typename C, unsigned int N> struct IsField<fadbad::F<C, N> >
{
	static const int value = numeric_limits<C>::is_specialized; ///<Is C a field?
};

template <int N, typename T, typename A, unsigned D>
inline Vector<N, T> get_value( const Vector<N, fadbad::F<T, D>, A> & v ){
	Vector<N,T> result(v.size());
	for(int i = 0; i < result.size(); ++i)
		result[i] = v[i].val();
	return result;
}

template <typename P, int N, typename A>
inline Vector<N, fadbad::F<P> > make_fad_vector( const Vector<N, P, A> & val, const unsigned start = 0, const unsigned size = N ){
	Vector<N, fadbad::F<P> > result = val;
	for(unsigned i = 0, d = start; i < val.size() && d < size; ++i, ++d)
		result[i].diff(d,size);
	return result;
}

template <unsigned D, typename P, int N, typename A>
inline Vector<N, fadbad::F<P,D> > make_fad_vector( const Vector<N, P, A> & val, const unsigned start = 0 ){
	Vector<N, fadbad::F<P,D> > result = val;
	for(unsigned i = 0, d = start; i < unsigned(val.size()) && d < D; ++i, ++d)
		result[i].diff(d);
	return result;
}

template <unsigned D, typename P, int N, typename A>
inline Vector<N, P> get_derivative( const Vector<N, fadbad::F<P,D>, A> & val, const int dim ){
	Vector<N, P> r;
	for(int i = 0; i < N; ++i)
		r[i] = val[i].deriv(dim);
	return r;
}

template <unsigned D, typename P, int N, typename A>
inline Matrix<N, D, P> get_jacobian( const Vector<N, fadbad::F<P, D>,  A> & val ){
	Matrix<N, D, P> result(N, val[0].size());
	for(unsigned i = 0; i < val[0].size(); ++i)
		result.T()[i] = get_derivative( val, i );
	return result;
}

template <int R, int C, typename P, unsigned D, typename A>
inline Matrix<R, C, P> get_derivative( const Matrix<R,C, fadbad::F<P, D>, A> & val, const int dim ){
	Matrix<R, C, P> result;
	for(int r = 0; r < R; ++r)
		for(int c = 0; c < C; ++c)
			result[r][c] = val[r][c].deriv(dim);
	return result;
}

template <typename P>
inline SO3<fadbad::F<P> > make_fad_so3(int start = 0, int size = 3){
	// const Vector<3, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0), start, size);
	// return SO3<fadbad::F<double> >(p);
	SO3<fadbad::F<P> > r;
	// this is a hack
	Matrix<3,3,fadbad::F<P> > & m = const_cast<Matrix<3,3,fadbad::F<P> > &>(r.get_matrix());
	m(2,1).diff(start, size);
	m(1,2) = m(2,1) * -1;

	m(0,2).diff(start+1, size);
	m(2,0) = m(0,2) * -1;

	m(1,0).diff(start+2, size);
	m(0,1) = m(1,0) * -1;

	return r;
}

template <typename P, unsigned D>
inline SO3<fadbad::F<P, D> > make_fad_so3(int start = 0){
	// const Vector<3, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0), start, size);
	// return SO3<fadbad::F<double> >(p);
	SO3<fadbad::F<P, D> > r;
	// this is a hack
	Matrix<3,3,fadbad::F<P, D> > & m = const_cast<Matrix<3,3,fadbad::F<P, D> > &>(r.get_matrix());
	m(2,1).diff(start);
	m(1,2) = m(2,1) * -1;

	m(0,2).diff(start+1);
	m(2,0) = m(0,2) * -1;

	m(1,0).diff(start+2);
	m(0,1) = m(1,0) * -1;

	return r;
}

template <typename P>
inline SE3<fadbad::F<P> > make_fad_se3( int start = 0, int size = 6){
	return SE3<fadbad::F<P> >(make_fad_so3<P>( start+3, size ), make_fad_vector(makeVector<P>(0.0, 0, 0), start, size));
}

template <typename P, unsigned D>
inline SE3<fadbad::F<P, D> > make_fad_se3( int start = 0){
	return SE3<fadbad::F<P, D> >(make_fad_so3<P, D>( start+3 ), make_fad_vector<D>(makeVector<P>(0.0, 0, 0), start));
}

template <typename P>
inline SE2<fadbad::F<P> > make_fad_se2(int start = 0, int size = 3) {
	return SE2<fadbad::F<P> >(make_fad_vector(makeVector<P>(0.0, 0, 0), start, size));
}

template <typename P, unsigned D>
inline SE2<fadbad::F<P, D> > make_fad_se2(int start = 0) {
	return SE2<fadbad::F<P, D> >(make_fad_vector<D>(makeVector<P>(0.0, 0, 0), start));
}

template <typename P>
inline SO2<fadbad::F<P> > make_fad_so2(int start = 0, int size = 1) {
	fadbad::F<P> r = 0;
	r.diff(start,size) = 1;
	return SO2<fadbad::F<P> >(r);
}

template <typename P, unsigned D>
inline SO2<fadbad::F<P, D> > make_fad_so2(int start = 0) {
	fadbad::F<P, D> r = 0;
	r.diff(start) = 1;
	return SO2<fadbad::F<P, D> >(r);
}

template <typename P>
inline SO3<fadbad::F<P> > make_left_fad_so3( const SO3<P> & r, int start = 0, int size = 3 ){
	return make_fad_so3<P>(start, size) * r;
}

template <typename P, unsigned D>
inline SO3<fadbad::F<P, D> > make_left_fad_so3( const SO3<P> & r, int start = 0){
	return make_fad_so3<P, D>(start) * r;
}

template <typename P>
inline SE3<fadbad::F<P> > make_left_fad_se3( const SE3<P> & t,  int start = 0, int size = 6 ){
	return make_fad_se3<P>(start, size) * t;
}

template <typename P, unsigned D>
inline SE3<fadbad::F<P, D> > make_left_fad_se3( const SE3<P> & t,  int start = 0){
	return make_fad_se3<P, D>(start) * t;
}

template <typename P>
inline SO2<fadbad::F<P> > make_left_fad_so2( const SO2<P> & r, int start = 0, int size = 1 ){
	return make_fad_so2<P>(start, size) * r;
}

template <typename P, unsigned D>
inline SO2<fadbad::F<P, D> > make_left_fad_so2( const SO2<P> & r, int start = 0 ){
	return make_fad_so2<P, D>(start) * r;
}

template <typename P>
inline SE2<fadbad::F<P> > make_left_fad_se2( const SE2<P> & t, int start = 0, int size = 3 ){
	return make_fad_se2<P>(start, size) * t;
}

template <typename P, unsigned D>
inline SE2<fadbad::F<P, D> > make_left_fad_se2( const SE2<P> & t, int start = 0){
	return make_fad_se2<P, D>(start) * t;
}

}

#endif // TOON_FADBAD_INTEGATION_H

