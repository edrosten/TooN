//-*- c++ -*-
//
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

template<int Size=-1, typename Precision=DefaultPrecision, typename Base=Internal::VBase>
class DiagonalMatrix {
public:
	inline DiagonalMatrix() {}
	inline DiagonalMatrix(int size_in) : my_vector(size_in) {}
	inline DiagonalMatrix(Precision* data) : my_vector(data) {}
	inline DiagonalMatrix(Precision* data, int size_in) : my_vector(data,size_in) {}
	inline DiagonalMatrix(Precision* data_in, int size_in, int stride_in, Internal::Slicing)
		: my_vector(data_in, size_in, stride_in, Internal::Slicing() ) {}


	// constructors to allow return value optimisations
	// construction from 0-ary operator
	template <class Op>
	inline DiagonalMatrix(const Operator<Op>& op)
		: my_vector (op)
	{
		op.eval(*this);
	}

	// constructor from arbitrary vector
	template<int Size2, typename Precision2, typename Base2>
	inline DiagonalMatrix(const Vector<Size2,Precision2,Base2>& from)
		: my_vector(from.size())
	{
		my_vector=from;
	}


	Precision& operator[](int i){return my_vector[i];}
	const Precision& operator[](int i) const {return my_vector[i];}

	typename Vector<Size, Precision, Base>::as_slice_type diagonal_slice() {
		return my_vector.as_slice();
	}


	Vector<Size,Precision,Base> my_vector;
};
