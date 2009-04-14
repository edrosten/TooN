// -*- c++ -*-

// Copyright (C) 2005,2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk), Gerhard Reitmayr (gr281@cam.ac.uk)
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


#ifndef TOON_INCLUDE_HELPERS_H
#define TOON_INCLUDE_HELPERS_H

#include <TooN/TooN.h>
#include <cmath>

namespace TooN {


	template<int Size, class Precision, class Base> void Fill(Vector<Size, Precision, Base>& v, const Precision& p)
	{
		for(int i=0; i < v.size(); i++)
			v[i]= p;
	}

	template<int Rows, int Cols, class Precision, class Base> void Fill(Matrix<Rows, Cols, Precision, Base>& m, const Precision& p)
	{
		for(int i=0; i < m.num_rows(); i++)
			for(int j=0; j < m.num_cols(); j++)
				m[i][j] = p;
	}

	template<int Size, class Precision, class Base> inline Precision norm(const Vector<Size, Precision, Base>& v)
	{
		using std::sqrt;
		return sqrt(v*v);
	}

	template<int Size, class Precision, class Base> inline Precision norm_sq(const Vector<Size, Precision, Base>& v)
	{
		v*v;
	}

	template<int Size, class Precision, class Base> inline Vector<Size, Precision> unit(const Vector<Size, Precision, Base> & v)
	{
		using std::sqrt;
		return v * (1/sqrt(v*v));
	}

	template<int Size, typename Precision, typename Base> inline Vector<Size-1, Precision> project( const Vector<Size, Precision, Base> & v){
		return v.template slice<0,Size-1>() / v[Size-1];
	}
	
	template<typename Precision, typename Base> inline Vector<-1, Precision> project( const Vector<-1, Precision, Base> & v){
		return v.slice(0,v.size()-1) / v[v.size()-1];
	}
	
	template<int Size, typename Precision, typename Base> inline Vector<Size+1, Precision> unproject( const Vector<Size, Precision, Base> & v){
		Vector<Size+1, Precision> result;
		result.template slice<0,Size>() = v;
		result[Size] = 1;
		return result;
	}
	
	template<typename Precision, typename Base> inline Vector<-1, Precision> unproject( const Vector<-1, Precision, Base> & v){
		Vector<-1, Precision> result(v.size()+1);
		result.slice(0,v.size()) = v;
		result[v.size()] = 1;
		return result;
	}

	namespace Internal{
		struct Copy
		{
			template<int R, int C, class P, class B, class Data> static void eval(Matrix<R, C, P, B>& m, const Data * data)
			{
				for(int r=0; r < m.num_rows(); r++)
					for(int c=0; c < m.num_rows(); c++)
						m[r][c] = *data++;
			}
		};

		// dummy structs that are used in 0-ary operators
		struct Zero;
		struct SizedZero;
		struct RCZero;
		struct Identity;
		struct SizedIdentity;
	}

	template<> struct Operator<Internal::RCZero> {
		Operator(int r, int c) : my_rows(r), my_cols(c) {}

		const int my_rows;
		const int my_cols;

		int num_rows() const {return my_rows;}
		int num_cols() const {return my_cols;}

		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}
	};


	template<> struct Operator<Internal::SizedZero> {

		// no idea why this doesn't indent properly
		Operator(int s)	: my_size(s) {}
		
		const int my_size;
		

		int size() const {return my_size;}
		int num_rows() const {return my_size;}
		int num_cols() const {return my_size;}

		template<int Size, class Precision, class Base>
		void eval(Vector<Size, Precision, Base>& v) const {
			for(int i=0; i < v.size(); i++) {
				v[i]= 0;
			}
		}

		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}
		
	};


	template<> struct Operator<Internal::Zero> {
		template<int Size, class Precision, class Base>
		void eval(Vector<Size, Precision, Base>& v) const {
			for(int i=0; i < v.size(); i++) {
				v[i]= 0;
			}
		}

		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}

		Operator<Internal::SizedZero> operator()(int s){
			return Operator<Internal::SizedZero>(s);
		}

		Operator<Internal::RCZero> operator()(int r, int c){
			return Operator<Internal::RCZero>(r,c);
		}

	};

	template<> struct Operator<Internal::SizedIdentity> {
		Operator(int s)	: my_size(s) {}
		
		const int my_size;
		
		int num_rows() const {return my_size;}
		int num_cols() const {return my_size;}

		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());

			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
						
			for(int r=0; r < m.num_rows(); r++) {
				m(r,r) = 1;
			}
		}
	};
		


	template<> struct Operator<Internal::Identity> {

		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());

			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
						
			for(int r=0; r < m.num_rows(); r++) {
				m(r,r) = 1;
			}
		}

		Operator<Internal::SizedIdentity> operator()(int s){
			return Operator<Internal::SizedIdentity>(s);
		}
	};


	static Operator<Internal::Zero> Zero;
	static Operator<Internal::Identity> Identity;

	/// row sum norm of input matrix m
	/// computes the maximum of the sums of absolute values over rows
	template <int R, int C, typename P, typename B>
	P inline norm_inf( const Matrix<R,C,P,B> & m ){
		using std::abs;
		using std::max;
		P n = 0;
		for(int r = 0; r < m.num_rows(); ++r){
			P s = 0;
			for(int c = 0; c < m.num_cols(); ++c)
				s += abs(m(r,c));
			n = max(n,s);
		}
		return n;
	}
	
	/// col sum norm of input matrix m
	/// computes the maximum of the sums of absolute values over columns
	template <int R, int C, typename P, typename B>
	P inline norm_1( const Matrix<R,C,P,B> & m ){
		using std::abs;
		using std::max;
		P n = 0;
		for(int c = 0; c < m.num_cols(); ++c){
			P s = 0;
			for(int r = 0; r < m.num_rows(); ++r)
				s += abs(m(r,c));
			n = max(n,s);
		}
		return n;
	}

	namespace Internal {
		template <int R, int C, typename P, typename B>
		inline Matrix<R, C, P> exp_taylor( const Matrix<R,C,P,B> & m ){
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());
			Matrix<R,C,P> result = TooN::Zero(m.num_rows(), m.num_cols());
			Matrix<R,C,P> f = TooN::Identity(m.num_rows());
			P k = 1;
			while(norm_inf((result+f)-result) > 0){
				result += f;
				f = (m * f) / k;
				k += 1;
			}
			return result;
		}
	};
	
	/// computes the matrix exponential of a matrix m by 
	/// scaling m by 1/(powers of 2), using Taylor series and 
	/// squaring again.
	/// @param m input matrix, must be square
	/// @return result matrix of the same size/type as input
	template <int R, int C, typename P, typename B>
	inline Matrix<R, C, P> exp( const Matrix<R,C,P,B> & m ){
		using std::max;
		SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());
		const P l = log2(norm_inf(m));
		const int s = max(0,(int)ceil(l));
		Matrix<R,C,P> result = Internal::exp_taylor(m/(1<<s));
		for(int i = 0; i < s; ++i)
			result = result * result;
		return result;
	}
	
	/// Returns true if every element is finite
	template<int S, class P, class B> bool isfinite(const Vector<S, P, B>& v)
	{ 
		using std::isfinite;
		for(int i=0; i < v.size(); i++)
			if(!isfinite(v[i]))
				return 0;
		return 1;
	}

	/// Returns true if any element is NaN
	template<int S, class P, class B> bool isnan(const Vector<S, P, B>& v)
	{ 
		using std::isnan;
		for(int i=0; i < v.size(); i++)
			if(isnan(v[i]))
				return 1;
		return 0;
	}
}
#endif
