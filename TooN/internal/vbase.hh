// -*- c++ -*-

// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

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

namespace TooN {

namespace Internal{
template<int Size, class Precision, int Stride, class Mem> struct GenericVBase;

////////////////////////////////////////////////////////////////////////////////
//
// Slice holding class
//
struct Default{};

template<int Stride, class Ptr=Default, class CPtr=Default, class Ref=Default, class CRef=Default>
struct SliceVBase {

	// this class is really just a typedef
	template<int Size, typename Precision>
	struct VLayout
		: public GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision, Ptr, CPtr, Ref, CRef> > {
		typedef typename VectorSlice<Size, Precision, Ptr, CPtr, Ref, CRef>::PointerType PointerType;
	
		VLayout(PointerType d, int length, int stride)
			:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision, Ptr, CPtr, Ref, CRef> >(d, length, stride){
		}

		template<class Op>
		VLayout(const Operator<Op>& op)
			:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(op) {}
	};

};

template<int Stride>
struct SliceVBase<Stride, Default, Default, Default, Default> {

	// this class is really just a typedef
	template<int Size, typename Precision>
	struct VLayout
		: public GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> > {

		typedef typename VectorSlice<Size, Precision>::PointerType PointerType;
	
		VLayout(PointerType d, int length, int stride)
			:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(d, length, stride){
		}

		template<class Op>
		VLayout(const Operator<Op>& op)
			:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(op) {}
	};

};

////////////////////////////////////////////////////////////////////////////////
//
// Classes for Vectors owning memory
//

struct VBase {

	// this class is really just a typedef
	template<int Size, class Precision>
	struct VLayout 
		: public GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> > {
	
		VLayout(){}

		VLayout(VLayout&&) = default;
		VLayout(const VLayout&) = default;

		VLayout(std::initializer_list<Precision> i)
			:GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >(i)
		{}
		
		template<typename Precision2, int Size2>
		VLayout(const Precision2(&i)[Size2])
			:GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >(i)
		{}

		VLayout(int s)
			:GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >(s)
		{}

		template<class Op>
		VLayout(const Operator<Op>& op)
			:GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >(op) {}
	};
};

////////////////////////////////////////////////////////////////////////////////
//
// Generic implementation
//

template<int Size, typename Precision, int Stride, typename Mem> struct GenericVBase: public Mem, public StrideHolder<Stride>
{	
	int stride() const{
		return StrideHolder<Stride>::stride();
	}

	//Optional constuctors
	GenericVBase(){}

	GenericVBase(GenericVBase&&) = default;
	GenericVBase(const GenericVBase&) = default;

	GenericVBase(int s)
	:Mem(s)
	{}

	template<typename Precision2, int Size2>
	GenericVBase(const Precision2(&i)[Size2])
	:Mem(i)
	{}

	GenericVBase(std::initializer_list<Precision> i)
	:Mem(i)
	{}

	typedef typename Mem::PointerType PointerType;
	typedef typename Mem::ConstPointerType ConstPointerType;
	typedef typename Mem::ReferenceType ReferenceType;
	typedef typename Mem::ConstReferenceType ConstReferenceType;

	GenericVBase(PointerType d, int length, int stride)
	:Mem(d, length),StrideHolder<Stride>(stride){
	}
	
	template<class Op>
	GenericVBase(const Operator<Op> & op) : Mem(op), StrideHolder<Stride>(op) {}

	using Mem::data;
	using Mem::size;

	ReferenceType operator[](int i) {
		Internal::check_index(size(), i);
		return data()[i * stride()];
	}

	ConstReferenceType operator[](int i) const {
		Internal::check_index(size(), i);
		return data()[i * stride()];
	}

	typedef SliceVBase<Stride, PointerType, ConstPointerType, ReferenceType, ConstReferenceType> SliceBase;
	typedef SliceVBase<Stride, ConstPointerType, ConstPointerType, ConstReferenceType, ConstReferenceType> ConstSliceBase;


	//Completely generic Vector slice operations below:
	template<int Start, int Length> 
	Vector<Length, Precision, SliceBase> slice(int start, int length){
		Internal::CheckSlice<Size, Start, Length>::check(size(), start, length);	
		return Vector<Length, Precision, SliceBase>(data() + stride() * (Start==Dynamic?start:Start), Length==Dynamic?length:Length, stride(), Slicing());
	}

	template<int Start, int Length> 
	const Vector<Length, const Precision, ConstSliceBase> slice(int start, int length) const{
		Internal::CheckSlice<Size, Start, Length>::check(size(), start, length);	
		return Vector<Length, const Precision, ConstSliceBase>(data() + stride() * (Start==Dynamic?start:Start), Length==Dynamic?length:Length, stride(), Slicing());
	}

	

	//Special case slice operations
	template<int Start, int Length> Vector<Length, Precision, SliceBase> slice(){
		Internal::CheckSlice<Size, Start, Length>::check();
		return slice<Start, Length>(Start, Length);
	}

	template<int Start, int Length> const Vector<Length, const Precision, ConstSliceBase> slice() const {
		Internal::CheckSlice<Size, Start, Length>::check();
		return slice<Start, Length>(Start, Length);
	}

	Vector<Dynamic, Precision, SliceBase> slice(int start, int length){
		return slice<Dynamic, Dynamic>(start, length);
	}

	const Vector<Dynamic, const Precision, ConstSliceBase> slice(int start, int length) const{
		return slice<Dynamic, Dynamic>(start, length);
	}
		
	//Other slices below
	const Matrix<1, Size, const Precision, Slice<1,Stride> > as_row() const{
		return Matrix<1, Size, const Precision, Slice<1,Stride> >(data(), 1, size(), 1, stride(), Slicing());
	}

	Matrix<1, Size, Precision, Slice<1,Stride> > as_row(){
		return Matrix<1, Size, Precision, Slice<1,Stride> >(data(), 1, size(), 1, stride(), Slicing());
	}

	const Matrix<Size, 1, const Precision, Slice<Stride,1> > as_col() const{
		return Matrix<Size, 1, const Precision, Slice<Stride,1> >(data(), size(), 1, stride(), 1, Slicing());
	}

	Matrix<Size, 1, Precision, Slice<Stride,1> > as_col(){
		return Matrix<Size, 1, Precision, Slice<Stride,1> >(data(), size(), 1, stride(), 1, Slicing());
	}

	typedef Vector<Size, Precision, SliceBase> as_slice_type;
	
	Vector<Size, Precision, SliceBase> as_slice(){                 
		return Vector<Size, Precision, SliceBase>(data(), size(), stride(), Slicing());         
	}

	const Vector<Size, const Precision, ConstSliceBase> as_slice() const {                 
		return Vector<Size, const Precision, ConstSliceBase>(data(), size(), stride(), Slicing());         
	}

	DiagonalMatrix<Size,Precision, SliceBase> as_diagonal() {
		return DiagonalMatrix<Size, Precision, SliceBase> (data(), size(), stride(), Slicing());
	}

	const DiagonalMatrix<Size,const Precision, ConstSliceBase> as_diagonal() const {
		return DiagonalMatrix<Size, const Precision, ConstSliceBase> (data(), size(), stride(), Slicing());
	}

};

}

}
