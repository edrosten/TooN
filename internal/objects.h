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

namespace TooN {

namespace Internal{
	// dummy structs that are used in 0-ary operators
	struct Zero;
	struct SizedZero;
	struct RCZero;
	template<class P> struct Identity;
	template<class P> struct SizedIdentity;

	template<int S, class P, class B, class Ps> class ScalarsVector;	
	template<int R, int C, class P, class B, class Ps> class ScalarsMatrix;	
	template<int R, int C, class P, class B, class Ps> class AddIdentity;	
	template<class P> class Scalars;	
	template<class P> class SizedScalars;	
	template<class P> class RCScalars;
}

////////////////////
// Zero
////////////////////



template<> struct Operator<Internal::SizedZero>;
template<> struct Operator<Internal::RCZero>;

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

	Operator<Internal::SizedZero> operator()(int s);
	Operator<Internal::RCZero> operator()(int r, int c);
	
};

template<> struct Operator<Internal::RCZero> : public Operator<Internal::Zero> {
	Operator(int r, int c) : my_rows(r), my_cols(c) {}

	const int my_rows;
	const int my_cols;

	int num_rows() const {return my_rows;}
	int num_cols() const {return my_cols;}
};

template<> struct Operator<Internal::SizedZero> : public Operator<Internal::Zero> {

	Operator(int s)	: my_size(s) {}
		
	const int my_size;
	
	int size() const {return my_size;}
	int num_rows() const {return my_size;}
	int num_cols() const {return my_size;}
};

inline Operator<Internal::SizedZero> Operator<Internal::Zero>::operator()(int s){
	return Operator<Internal::SizedZero>(s);
}

inline Operator<Internal::RCZero> Operator<Internal::Zero>::operator()(int r, int c){
	return Operator<Internal::RCZero>(r,c);
}


//////////////
// Identity
//////////////

//Operator to construct a new matrix with idendity added 
template<int R, int C, class P, class B, class Precision> struct Operator<Internal::AddIdentity<R,C,P,B,Precision> >
{
	const Precision s;
	const Matrix<R,C,P,B>& m;
	bool invert_m;

	Operator(Precision s_, const Matrix<R,C,P,B>& m_, bool b)
		:s(s_),m(m_),invert_m(b){}
	template<int R1, int C1, class P1, class B1>
	void eval(Matrix<R1,C1,P1,B1>& mm) const{
		for(int r=0; r < m.num_rows(); r++)
			for(int c=0; c < m.num_cols(); c++)
				if(invert_m)
					mm[r][c] = -m[r][c];
				else
					mm[r][c] = m[r][c];

		for(int i=0; i < m.num_rows(); i++)
				mm[i][i] += s;
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


template<class Pr> struct Operator<Internal::Identity<Pr> > {
	
	typedef Pr Precision;
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

	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> > add(const Matrix<Rows,Cols, P1, B1>& m) const
	{
		SizeMismatch<Rows, Cols>::test(m.num_rows(), m.num_cols());
		return Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> >(val, m, 0);
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> > rsubtract(const Matrix<Rows,Cols, P1, B1>& m) const
	{
		SizeMismatch<Rows, Cols>::test(m.num_rows(), m.num_cols());
		return Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> >(-val, m, 0);
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> > lsubtract(const Matrix<Rows,Cols, P1, B1>& m) const
	{
		SizeMismatch<Rows, Cols>::test(m.num_rows(), m.num_cols());
		return Operator<Internal::AddIdentity<Rows,Cols,P1,B1,Precision> >(val, m, 1);
	}


	Operator<Internal::SizedIdentity<Precision> > operator()(int s){
		return Operator<Internal::SizedIdentity<Precision> >(s);
	}

	template<class Pout, class Pmult> Operator<Internal::Identity<Pout> > scale_me(const Pmult& m) const
	{
		return Operator<Internal::Identity<Pout> >(val*m);
	}


};
	

template<class Precision> struct Operator<Internal::SizedIdentity<Precision> > 
	: public  Operator<Internal::Identity<Precision> > {

	using Operator<Internal::Identity<Precision> >::val;
	const int my_size;

	Operator(int s, const Precision& v=1)
		:Operator<Internal::Identity<Precision> > (v), my_size(s)
	{}

	int num_rows() const {return my_size;}
	int num_cols() const {return my_size;}

	template<class Pout, class Pmult> Operator<Internal::SizedIdentity<Pout> > scale_me(const Pmult& m) const
	{
		return Operator<Internal::SizedIdentity<Pout> >(my_size, val*m);
	}
};
////////////////////////////////////////////////////////////////////////////////
//
// Addition of scalars to vectors and matrices
//

	
//Operator to construct a new vector a a vector with a scalar added to every element
template<int S, class P, class B, class Precision> struct Operator<Internal::ScalarsVector<S,P,B,Precision> >
{
	const Precision s;
	const Vector<S,P,B>& v;
	const bool invert_v;
	Operator(Precision s_, const Vector<S,P,B>& v_, bool inv)
		:s(s_),v(v_),invert_v(inv){}

	template<int S1, class P1, class B1>
	void eval(Vector<S1,P1,B1>& vv) const{
		for(int i=0; i < v.size(); i++)
			if(invert_v)
				vv[i] = s - v[i];
			else
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
	const bool invert_m;
	Operator(Precision s_, const Matrix<R,C,P,B>& m_, bool inv)
		:s(s_),m(m_),invert_m(inv){}
	template<int R1, int C1, class P1, class B1>
	void eval(Matrix<R1,C1,P1,B1>& mm) const{
		for(int r=0; r < m.num_rows(); r++)
			for(int c=0; c < m.num_cols(); c++)
				if(invert_m)
					mm[r][c] = s - m[r][c];
				else
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

	////////////////////////////////////////
	//
	// All applications for vector
	//

	template <int Size, typename P1, typename B1> 
	void eval(Vector<Size, P1, B1>& v) const
	{
		for(int i=0; i < v.size(); i++)
			v[i] = s;
	}

	template <int Size, typename P1, typename B1> 
	void plusequals(Vector<Size, P1, B1>& v) const
	{
		for(int i=0; i < v.size(); i++)
			v[i] += s;
	}

	template <int Size, typename P1, typename B1> 
	Operator<Internal::ScalarsVector<Size,P1,B1,Precision> > add(const Vector<Size, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsVector<Size,P1,B1,Precision> >(s, v, 0);
	}

	template <int Size, typename P1, typename B1> 
	Operator<Internal::ScalarsVector<Size,P1,B1,Precision> > rsubtract(const Vector<Size, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsVector<Size,P1,B1,Precision> >(-s, v, 0);
	}

	template <int Size, typename P1, typename B1> 
	Operator<Internal::ScalarsVector<Size,P1,B1,Precision> > lsubtract(const Vector<Size, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsVector<Size,P1,B1,Precision> >(s, v, 1);
	}

	////////////////////////////////////////
	//
	// All applications for matrix
	//

	template <int Rows, int Cols, typename P1, typename B1> 
	void eval(Matrix<Rows,Cols, P1, B1>& m) const
	{
		for(int r=0; r < m.num_rows(); r++)
			for(int c=0; c < m.num_cols(); c++)
				m[r][c] = s;
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	void plusequals(Matrix<Rows,Cols, P1, B1>& m) const
	{
		for(int r=0; r < m.num_rows(); r++)
			for(int c=0; c < m.num_cols(); c++)
				m[r][c] += s;
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	void minusequals(Matrix<Rows,Cols, P1, B1>& m) const
	{
		for(int r=0; r < m.num_rows(); r++)
			for(int c=0; c < m.num_cols(); c++)
				m[r][c] -= s;
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> > add(const Matrix<Rows,Cols, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> >(s, v, 0);
	}


	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> > rsubtract(const Matrix<Rows,Cols, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> >(-s, v, 0);
	}

	template <int Rows, int Cols, typename P1, typename B1> 
	Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> > lsubtract(const Matrix<Rows,Cols, P1, B1>& v) const
	{
		return Operator<Internal::ScalarsMatrix<Rows,Cols,P1,B1,Precision> >(s, v, 1);
	}
	////////////////////////////////////////
	//
	// Create sized versions for initialization
	//

	Operator<Internal::SizedScalars<Precision> > operator()(int size) const
	{
		return Operator<Internal::SizedScalars<Precision> > (s,size);
	}

	Operator<Internal::RCScalars<Precision> > operator()(int r, int c) const
	{
		return Operator<Internal::RCScalars<Precision> > (s,r,c);
	}

	template<class Pout, class Pmult> Operator<Internal::Scalars<Pout> > scale_me(const Pmult& m) const
	{
		return Operator<Internal::Scalars<Pout> >(s*m);
	}
};

template<class P> struct Operator<Internal::SizedScalars<P> >: public Operator<Internal::Scalars<P> >
{
	using Operator<Internal::Scalars<P> >::s;
	const int my_size;
	int size() const {
		return my_size;
	}
		
	Operator(P s, int sz)
		:Operator<Internal::Scalars<P> >(s),my_size(sz){}
		
	template<class Pout, class Pmult> Operator<Internal::SizedScalars<Pout> > scale_me(const Pmult& m) const
	{
		return Operator<Internal::SizedScalars<Pout> >(s*m, my_size);
	}

private:
	void operator()(int);
	void operator()(int,int);
};

		

template<class P> struct Operator<Internal::RCScalars<P> >: public Operator<Internal::Scalars<P> >
{
	using Operator<Internal::Scalars<P> >::s;
	const int my_rows, my_cols;
	int num_rows() const {
		return my_rows;
	}
	int num_cols() const {
		return my_cols;
	}
		
	Operator(P s, int r, int c)
		:Operator<Internal::Scalars<P> >(s),my_rows(r),my_cols(c)
	{}
		
	template<class Pout, class Pmult> Operator<Internal::RCScalars<Pout> > scale_me(const Pmult& m) const
	{
		return Operator<Internal::RCScalars<Pout> >(s*m, my_rows, my_cols);
	}

private:
	void operator()(int);
	void operator()(int,int);
};


////////////////////////////////////////////////////////////////////////////////
//
// How to scale scalable operators
//
	
template<template<class> class Op, class Pl, class Pr> 
Operator<Op<typename Internal::MultiplyType<Pl, Pr>::type > >
operator*(const Pl& l, const Operator<Op<Pr> >& r)
{
	return r.template scale_me<typename Internal::MultiplyType<Pl, Pr>::type, Pl>(l); 
}

template<template<class> class Op, class Pl, class Pr> 
Operator<Op<typename Internal::MultiplyType<Pl, Pr>::type > >
operator*(const Operator<Op<Pl> >& l, const Pr&  r)
{
	return l.template scale_me<typename Internal::MultiplyType<Pl, Pr>::type, Pl>(r); 
}

template<template<class> class Op, class Pl, class Pr> 
Operator<Op<typename Internal::DivideType<Pl, Pr>::type > >
operator/(const Operator<Op<Pl> >& l, const Pr&  r)
{
	return l.template scale_me<typename Internal::MultiplyType<Pl, Pr>::type, Pl>(static_cast<typename Internal::DivideType<Pl,Pr>::type>(1)/r); 
}


template<class Op, class P>
Operator<Op> operator-(const Operator<Op>& o)
{
	return o.scalar(-1);
}


/**This function us used to add a scalar to every element of a vector or
   matrix. For example:
   @code
   Vector<3> v;
   ...
   ...
   v += Ones * 3; //Add 3 to every element of v;
   @endcode
   Both + and += are supported on vectors,matrices and slices.
*/

static Operator<Internal::Zero> Zeros;
static Operator<Internal::Identity<DefaultPrecision> > Identity;
static const Operator<Internal::Scalars<DefaultPrecision> > Ones(1);

}
