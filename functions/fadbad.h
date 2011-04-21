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

#ifndef TOON_FADBAD_INTEGATION_H
#define TOON_FADBAD_INTEGATION_H

#include <iostream>

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/se2.h>

#include <FADBAD++/fadiff.h>

namespace fadbad {

template <typename P>
inline std::ostream & operator<<( std::ostream & out, const F<P> & val ){
	return out << val.val();
}

template <typename P>
inline F<P> abs( const F<P> & val ){
    return (val.val() < 0) ? -val : val;
}

}

namespace TooN {

template<class C> struct IsField<fadbad::F<C> >
{
	static const int value = numeric_limits<C>::is_specialized; ///<Is C a field?
};

template <int N, typename T, typename A>
inline Vector<N, T> get_value( const Vector<N, fadbad::F<T>, A> & v ){
    Vector<N,T> result(v.size());
    for(int i = 0; i < result.size(); ++i)
        result[i] = v[i].val();
    return result;
}

template <int N> 
inline Vector<N, fadbad::F<double> > make_fad_vector( const Vector<N, double> & val, const int start = 0, const int size = N ){
	using std::min;
	Vector<N, fadbad::F<double> > result = val;
	for(int i = 0, d = start; i < val.size() && d < size; ++i, ++d)
		result[i].diff(d,size);
	return result;
}

template <int N>
inline Vector<N, double> get_derivative( const Vector<N, fadbad::F<double> > & val, const int D ){
	Vector<N, double> r;
	for(int i = 0; i < N; ++i)
		r[i] = val[i].deriv(D);
	return r;
}

template <int N, int D>
inline Matrix<N, D, double> get_jacobian( const Vector<N, fadbad::F<double> > & val ){
	Matrix<N, D, double> result;
	for(int i = 0; i < D; ++i)
		result.T()[i] = get_derivative( val, i );
	return result;
}

template <int R, int C>
inline Matrix<R, C, double> get_derivative( const Matrix<R,C, fadbad::F<double> > & val, const int D ){
	Matrix<R, C, double> result;
	for(int r = 0; r < R; ++r)
		for(int c = 0; c < C; ++c)
			result[r][c] = val[r][c].deriv(D);
	return result;
}

inline SO3<fadbad::F<double> > make_fad_so3(int start = 0, int size = 3){
	const Vector<3, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0), start, size);
	return SO3<fadbad::F<double> >(p);
}

inline SE3<fadbad::F<double> > make_fad_se3( int start = 0, int size = 6){
	const Vector<6, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0, 0, 0, 0), start, size);
	return SE3<fadbad::F<double> >(p);
}

inline SE2<fadbad::F<double> > make_fad_se2(int start = 0, int size = 3) {
	return SE2<fadbad::F<double> >(make_fad_vector(makeVector(0.0, 0, 0), start, size));
}

inline SO2<fadbad::F<double> > make_fad_so2(int start = 0, int size = 1) {
	fadbad::F<double> r = 0.0;
	r.diff(start,size) = 1.0;
	return SO2<fadbad::F<double> >(r);
}

inline SO3<fadbad::F<double> > make_left_fad_so3( const SO3<double> & r, int start = 0, int size = 3 ){
	return make_fad_so3(start, size) * r;
}

inline SE3<fadbad::F<double> > make_left_fad_se3( const SE3<double> & t,  int start = 0, int size = 6 ){
	return make_fad_se3(start, size) * t;
}

inline SO2<fadbad::F<double> > make_left_fad_so2( const SO2<double> & r, int start = 0, int size = 1 ){
	return make_fad_so2(start, size) * r;
}

inline SE2<fadbad::F<double> > make_left_fad_se2( const SE2<double> & t, int start = 0, int size = 3 ){
	return make_fad_se2(start, size) * t;
}

}

#endif // TOON_FADBAD_INTEGATION_H

