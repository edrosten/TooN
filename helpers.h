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
		return v*v;
	}

	template<int Size, class Precision, class Base> inline Vector<Size, Precision> unit(const Vector<Size, Precision, Base> & v)
	{
		using std::sqrt;
		return v * (1/sqrt(v*v));
	}
	
	//Note because of the overload later, this function will ONLY receive sliced vectors. Therefore
	//a copy can be made, which is still a slice, so operating on the copy operates on the original
	//data.
	///Normalize a vector in place
	///@param v Vector to normalize
	///@ingroup gLinAlg
	template<int Size, class Precision, class Base> inline void normalize(Vector<Size, Precision, Base> v)
	{
		using std::sqrt;
		v /= sqrt(v*v);
	}
	
	//This overload is required to operate on non-slice vectors
	///@overload
	template<int Size, class Precision> inline void normalize(Vector<Size, Precision> & v)
	{
		normalize(v.as_slice());
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
		template<class P> struct Identity;
		template<class P> struct SizedIdentity;
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

	template<class Precision> struct Operator<Internal::SizedIdentity<Precision> > {

		const Precision val;
		const int my_size;

		Operator(int s, const Precision& v=1)
		:my_size(s),val(v)
		{}

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
				m(r,r) = val;
			}
		}
	};

	template<class Precision> struct Operator<Internal::Identity<Precision> > {

		const Precision val;
		Operator(const Precision& v=1)
		:val(v)
		{}
		
		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());

			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
						
			for(int r=0; r < m.num_rows(); r++) {
				m(r,r) = val;
			}
		}

		Operator<Internal::SizedIdentity<Precision> > operator()(int s){
			return Operator<Internal::SizedIdentity<Precision> >(s);
		}
	};
	
	template<class P1, class P2> Operator<Internal::Identity<typename Internal::MultiplyType<P1, P2>::type> > operator*(const P1& p, const Operator<Internal::Identity<P2> >& i)
	{
		return Operator<Internal::Identity<typename Internal::MultiplyType<P1, P2>::type> >(p * i.val);
	}

	template<class P1, class P2> Operator<Internal::Identity<typename Internal::MultiplyType<P1, P2>::type> > operator*(const Operator<Internal::Identity<P2> >& i, const P1&p)
	{
		return Operator<Internal::Identity<typename Internal::MultiplyType<P1, P2>::type> >(p * i.val);
	}

	template<class P1, class P2> Operator<Internal::SizedIdentity<typename Internal::MultiplyType<P1, P2>::type> > operator*(const P1& p, const Operator<Internal::SizedIdentity<P2> >& i)
	{
		return Operator<Internal::SizedIdentity<typename Internal::MultiplyType<P1, P2>::type> >(i.my_size, p * i.val);
	}

	template<class P1, class P2> Operator<Internal::SizedIdentity<typename Internal::MultiplyType<P1, P2>::type> > operator*(const Operator<Internal::SizedIdentity<P2> >& i, const P1&p)
	{
		return Operator<Internal::SizedIdentity<typename Internal::MultiplyType<P1, P2>::type> >(i.my_size, p * i.val);
	}

	static Operator<Internal::Zero> Zero;
	static Operator<Internal::Identity<double> > Identity;

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

	/// Symmetrize a matrix
	template<int Rows, int Cols, typename Precision, typename Base>
	void Symmetrize(Matrix<Rows,Cols,Precision,Base>& m){
		SizeMismatch<Rows,Cols>::test(m.num_rows(), m.num_cols());
		for(int r=0; r<m.num_rows()-1; r++){
			for(int c=r+1; c<m.num_cols(); c++){
				const Precision temp=(m(r,c)+m(c,r))/2;
				m(r,c)=temp;
				m(c,r)=temp;
			}
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////
	//
	// Addition of scalars ro vectors and matrices
	//

	namespace Internal{
		template<int S, class P, class B, class Ps> class ScalarsVector;	
		template<int R, int C, class P, class B, class Ps> class ScalarsMatrix;	
		template<class P> class Scalars;	
	}
	
	//Operator to construct a new vector a a vector with a scalar added to every element
	template<int S, class P, class B, class Precision> struct Operator<Internal::ScalarsVector<S,P,B,Precision> >
	{
		const Precision s;
		const Vector<S,P,B>& v;
		Operator(Precision s_, const Vector<S,P,B>& v_)
		:s(s_),v(v_){}

		template<int S1, class P1, class B1>
		void eval(Vector<S1,P1,B1>& vv) const{
			for(int i=0; i < v.size(); i++)
				vv[i] = s + v[i];
		}

		int size() const
		{
			return v.size();
		}
	};

	//Operator to construct a new matrix a a matrix with a scalar added to every element
	template<int R, int C, class P, class B, class Precision> struct Operator<Internal::ScalarsMatrix<R,C,P,B,Precision> >
	{
		const Precision s;
		const Matrix<R,C,P,B>& m;
		Operator(Precision s_, const Matrix<R,C,P,B>& m_)
		:s(s_),m(m_){}
		template<int R1, int C1, class P1, class B1>
		void eval(Matrix<R1,C1,P1,B1>& mm) const{
			for(int r=0; r < m.num_rows(); r++)
				for(int c=0; c < m.num_cols(); c++)
					mm[r][c] = s + m[r][c];
		}

		int num_rows() const
		{
			return m.num_rows();
		}
		int num_cols() const
		{
			return m.num_cols();
		}
	};

	//Generic scalars object. Knows how to be added, knows how to deal with +=
	template<class P> struct Operator<Internal::Scalars<P> >
	{
		typedef P Precision;
		const Precision s;
		Operator(Precision s_)
		:s(s_){}

		template <int Size, typename P1, typename B1> 
		void plusequals(Vector<Size, P1, B1>& v) const
		{
			for(int i=0; i < v.size(); i++)
				v[i] += s;
		}

		template <int Size, typename P1, typename B1> 
		Operator<Internal::ScalarsVector<Size,P1,B1,Precision> > add(const Vector<Size, P1, B1>& v) const
		{
			return Operator<Internal::ScalarsVector<Size,P1,B1,Precision> >(s, v);
		}

		template <int Rows, int Cols, typename P1, typename B1> 
		void plusequals(Matrix<Rows,Cols, P1, B1>& m) const
		{
			for(int r=0; r < m.num_rows(); r++)
				for(int c=0; c < m.num_cols(); c++)
					m[r][c] += s;
		}

		template <int Rows, int Cols, typename P1, typename B1> 
		Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> > add(const Matrix<Rows,Cols, P1, B1>& v) const
		{
			return Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> >(s, v);
		}
	};
	
	/**This function us used to add a scalar to every element of a vector or
	matrix. For example:
	@code
		Vector<> v;
		...
		...
		v += Scalars(3); //Add 3 to every element of v;
	@endcode
	Both + and += are supported on vectors,matrices and slices.
	@param  s Scalar to add.
	*/
	template<class P> Operator<Internal::Scalars<P> > Scalars(const P& s)
	{
		return Operator<Internal::Scalars<P> > (s);
	}
}
#endif
