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
class Vector : public Base::template VLayout<Size, Precision> {
public:
  // sneaky hack: only one of these constructors will work with any given base
  // class but they don't generate errors unless the user tries to use one of them
  // although the error message may be less than helpful - maybe this can be changed?
	inline Vector(){}
	inline Vector(int size_in) : Base::template VLayout<Size, Precision>(size_in) {}

	inline Vector(Precision* data) : Base::template VLayout<Size, Precision> (data) {}
	inline Vector(Precision* data, int size_in) : Base::template VLayout<Size, Precision> (data, size_in) {}

	// internal constructor
	inline Vector(Precision* data_in, int size_in, int stride_in, Internal::Slicing) : Base::template VLayout<Size, Precision>(data_in, size_in, stride_in) {}
	
	using Base::template VLayout<Size, Precision>::size;

	// constructors to allow return value optimisations
	// construction from 0-ary operator
	template <class Op>
	inline Vector(const Operator<Op>& op)
		: Base::template VLayout<Size, Precision> (op)
	{
		op.eval(*this);
	}

	// Copy construction is a very special case. Copy construction goes all the
	// way down to the bottom. GenericVBase has no idea how to copy itself.
	// However, the underlying allocator objects do.  In the case of static sized
	// objects, C++ automatically copies the data.  For slice objects, C++ copies
	// all parts (pointer and size), which is correct.  For dynamically sized
	// non-slice objects the copying has to be done by hand.
	
	// inline Vector(const Vector&from);

	// constructor from arbitrary vector
	template<int Size2, typename Precision2, typename Base2>
	inline Vector(const Vector<Size2,Precision2,Base2>& from):
		Base::template VLayout<Size, Precision>(from.size()) {
		operator=(from);
	}

	// assignment from a 0-ary operator
	template <class Op>
	inline Vector & operator=(const Operator<Op>& op){
		op.eval(*this);
		return *this;
	}

	// operator = from copy
	inline Vector& operator= (const Vector& from){
		SizeMismatch<Size,Size>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}

	// operator =
	template<int Size2, typename Precision2, typename Base2>
	Vector<Size,Precision,Base >& operator= (const Vector<Size2, Precision2, Base2>& from){
		SizeMismatch<Size,Size2>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}


	Vector& operator+=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]+=rhs;
		return *this;
	}
	

	Vector& operator-=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]-=rhs;
		return *this;
	}
	
	Vector& operator/=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]/=rhs;
		return *this;
	}
	

	Vector& operator*=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]*=rhs;
		return *this;
	}
	
	template<int Size2, class Precision2, class Base2>
	Vector& operator+=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]+=rhs[i];
		return *this;
	}

	template<int Size2, class Precision2, class Base2>
	Vector& operator-=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]-=rhs[i];
		return *this;
	}

	Vector& ref()
	{
		return *this;
	}

};
