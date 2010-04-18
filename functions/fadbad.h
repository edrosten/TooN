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

}

namespace TooN {

template<class C> struct IsField<fadbad::F<C> >
{
	static const int value = numeric_limits<C>::is_specialized; ///<Is C a field?
};

template <int N> 
inline Vector<N, fadbad::F<double> > make_fad_vector( const Vector<N, double> & val, const int start = 0, const int size = N ){
	using std::min;
	Vector<N, fadbad::F<double> > result = val;
	for(int i = start; i < min(size, start+N); ++i)
		result[i].diff(i,size);
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

inline SO3<fadbad::F<double> > make_fad_so3(){
	const Vector<3, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0));
	return SO3<fadbad::F<double> >(p);
}

inline SE3<fadbad::F<double> > make_fad_se3(){
	const Vector<6, fadbad::F<double> > p = make_fad_vector(makeVector(0.0, 0, 0, 0, 0, 0));
	return SE3<fadbad::F<double> >(p);
}

inline SE2<fadbad::F<double> > make_fad_se2() {
	return SE2<fadbad::F<double> >(make_fad_vector(makeVector(0.0, 0, 0)));
}

inline SO2<fadbad::F<double> > make_fad_so2() {
	fadbad::F<double> r = 0.0;
	r.diff(0,1) = 1.0;
	return SO2<fadbad::F<double> >(r);
}

inline SO3<fadbad::F<double> > make_left_fad_so3( const SO3<double> & r ){
	return make_fad_so3() * r;
}

inline SE3<fadbad::F<double> > make_left_fad_se3( const SE3<double> & t ){
	return make_fad_se3() * t;
}

inline SO2<fadbad::F<double> > make_left_fad_so2( const SO2<double> & r ){
	return make_fad_so2() * r;
}

inline SE2<fadbad::F<double> > make_left_fad_se2( const SE2<double> & t ){
	return make_fad_se2() * t;
}

}

#endif // TOON_FADBAD_INTEGATION_H

