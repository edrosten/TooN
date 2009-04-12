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

////////////////////////////////////////////////////////////////////////////////
//
// Helper classes for matrices constructed as references to foreign data
//

struct Reference
{

	template<int Size, typename Precision>
	struct VLayout
		: public Internal::GenericVBase<Size, Precision, 1, Internal::VectorSlice<Size, Precision> >
	{

		VLayout(Precision* p, int sz=0)
			: Internal::GenericVBase<Size, Precision, 1, Internal::VectorSlice<Size, Precision> >(p, sz, 0)
		{}
	};


	struct RowMajor
	{
		template<int Rows, int Cols, class Precision>
		struct MLayout
			: public Internal::GenericMBase<Rows, Cols, Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> >
		{

			MLayout(Precision* p)
				: Internal::GenericMBase<Rows,Cols,Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> > (p, 0, 0, 0, 0)
			{}
			MLayout(Precision* p, int r, int c)
				: Internal::GenericMBase<Rows,Cols,Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> > (p, r, c, 0, 0)
			{}
		};
	};

	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct MLayout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			MLayout(Precision* p)
				: Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p, 0, 0, 0, 0)
			{}
			MLayout(Precision* p, int r, int c)
				: Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p, r, c, 0, 0)
			{}
		};
	};
};

template<int R, int C, typename Precision=double, class Type=Reference::RowMajor> struct Wrap
{
	static Matrix<R, C, Precision, Type> wrap(Precision* p)
	{
		return Matrix<R, C, Precision, Type>(p);
	}
};


template<int R, typename Precision, class Type> struct Wrap<R, Dynamic, Precision, Type>
{
	static Matrix<R, Dynamic, Precision, Type> wrap(Precision* p, int cols)
	{
		return Matrix<R, Dynamic, Precision, Type>(p, 0, cols);
	}
};


template<int C, typename Precision, class Type> struct Wrap<Dynamic, C, Precision, Type>
{
	static Matrix<Dynamic, C, Precision, Type> wrap(Precision* p, int rows)
	{
		return Matrix<Dynamic, C, Precision, Type>(p, rows, 0);
	}
};


template<typename Precision, class Type> struct Wrap<Dynamic, Dynamic, Precision, Type>
{
	static Matrix<Dynamic, Dynamic, Precision, Type> wrap(Precision* p, int rows, int cols)
	{
		return Matrix<Dynamic, Dynamic, Precision, Type>(p, rows, cols);
	}
};

