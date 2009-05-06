// -*- c++ -*-

// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)
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

template <int Rows=-1, int Cols=Rows, class Precision=DefaultPrecision, class Layout = RowMajor>
struct Matrix : public Layout::template MLayout<Rows, Cols, Precision>
{
public:

	using Layout::template MLayout<Rows, Cols, Precision>::my_data;
	using Layout::template MLayout<Rows, Cols, Precision>::num_rows;
	using Layout::template MLayout<Rows, Cols, Precision>::num_cols;

	//Use Tom's sneaky constructor hack...
	Matrix(){}

	Matrix(int rows, int cols) :
		Layout::template MLayout<Rows,Cols,Precision>(rows, cols)
	{}

	Matrix(Precision* p) :
		Layout::template MLayout<Rows, Cols, Precision>(p)
	{}

	Matrix(Precision* p, int r, int c) :
		Layout::template MLayout<Rows, Cols, Precision>(p, r, c)
	{}

	// Internal constructor used by GenericMBase::slice(...)
	Matrix(Precision* data, int rows, int cols, int rowstride, int colstride, Internal::Slicing)
	:Layout::template MLayout<Rows, Cols, Precision>(data, rows, cols, rowstride, colstride){}

	//See vector.hh and allocator.hh for details about why the
	//copy constructor should be default.
	template <class Op>
	inline Matrix(const Operator<Op>& op)
		:Layout::template MLayout<Rows,Cols,Precision>(op)
	{
		op.eval(*this);
	}

	// constructors to allow return value optimisations
	// construction from 1-ary operator
	template <class T, class Op>
	inline Matrix(const T& arg, int rows, int cols, const Operator<Op>&) 
	:Layout::template MLayout<Rows,Cols,Precision>(rows, cols) 
	{
	    Op::eval(*this,arg);
	}

	// constructor from 2-ary operator
	template <class LHS, class RHS, class Op>
	inline Matrix(const LHS& lhs, const RHS& rhs, int rows, int cols, const Operator<Op>&)
	:Layout::template MLayout<Rows,Cols,Precision>(rows, cols)
	{
	    Op::eval(*this,lhs,rhs);
	}

	// constructor from arbitrary matrix
	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	inline Matrix(const Matrix<Rows2, Cols2,Precision2,Base2>& from)
	:Layout::template MLayout<Rows,Cols,Precision>(from.num_rows(), from.num_cols())
	{
	    operator=(from);
	}

	// operator = from copy
	inline Matrix& operator= (const Matrix& from)
	{
		SizeMismatch<Rows, Rows>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	// operator = 0-ary operator
	template<class Op> inline Matrix& operator= (const Operator<Op>& op)
	{
		op.eval(*this);
		return *this;
	}

	// operator =
	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	Matrix& operator*=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] *= rhs;

		  return *this;
	}

	Matrix& operator/=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] /= rhs;

		  return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator+= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] += from[r][c];

	    return *this;
	}

	template<class Op>
	Matrix& operator+=(const Operator<Op>& op)
	{
		op.plusequals(*this);
		return *this;
	}

	template<class Op>
	Matrix& operator-=(const Operator<Op>& op)
	{
		op.minusequals(*this);
		return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator-= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] -= from[r][c];

	    return *this;
	}

  	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	bool operator== (const Matrix<Rows2, Cols2, Precision2, Base2>& rhs)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), rhs.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), rhs.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
		    if((*this)[r][c] != rhs[r][c])
		      return 0;
	    return 1;
	}

  	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	bool operator!= (const Matrix<Rows2, Cols2, Precision2, Base2>& rhs)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), rhs.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), rhs.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
		    if((*this)[r][c] != rhs[r][c])
		      return 1;
	    return 0;
	}


	Matrix& ref()
	{
		return *this;
	}
};

}
