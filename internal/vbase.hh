// -*- c++ -*-

namespace Internal{
template<int Size, class Precision, int Stride, class Mem> struct GenericVBase;


////////////////////////////////////////////////////////////////////////////////
//
// Slice holding class
//

template<int Stride>
struct SliceVBase {

	// this class is really just a typedef
	template<int Size, typename Precision>
	struct VLayout
		: public GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> > {
	

		VLayout(Precision* d, int stride)
			:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(d, stride){
		}

		VLayout(Precision* d, int length, int stride)
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

	GenericVBase(int s)
	:Mem(s)
	{}

	GenericVBase(Precision* d, int length, int stride)
	:Mem(d, length),StrideHolder<Stride>(stride){
	}
	
	template<class Op>
	GenericVBase(const Operator<Op> & op) : Mem(op), StrideHolder<Stride>(op) {}

	using Mem::my_data;
	using Mem::size;

	Precision& operator[](int i) {
		Internal::check_index(size(), i);
		return my_data[i * stride()];
	}

	const Precision& operator[](int i) const {
		Internal::check_index(size(), i);
		return my_data[i * stride()];
	}


	template<int Start, int Length> 
	Vector<Length, Precision, SliceVBase<Stride> > slice(){
		Internal::CheckStaticSlice<Size, Start, Length>::check(size());
		return Vector<Length, Precision, SliceVBase<Stride> >(my_data + stride()*Start, Length, stride(), Slicing());
	}

	template<int Start, int Length> 
	const Vector<Length, Precision, SliceVBase<Stride> > slice() const {
		Internal::CheckStaticSlice<Size, Start, Length>::check(size());
		return Vector<Length, Precision, SliceVBase<Stride> >(const_cast<Precision*>(my_data + stride()*Start), Length, stride(), Slicing());
	}

	Vector<-1, Precision, SliceVBase<Stride> > slice(int start, int length){
		Internal::CheckDynamicSlice::check(size(), start, length);
		return Vector<-1, Precision, SliceVBase<Stride> >(my_data + stride()*start, length, stride(), Slicing());
	}

	const Vector<-1, Precision, SliceVBase<Stride> > slice(int start, int length) const {
		Internal::CheckDynamicSlice::check(size(), start, length);
		return Vector<-1, Precision, SliceVBase<Stride> >(const_cast<Precision*>(my_data + stride()*start), length, stride(), Slicing());
	}

	const Matrix<1, Size, Precision, Slice<1,Stride> > as_row() const{
		return Matrix<1, Size, Precision, Slice<1,Stride> >(const_cast<Precision*>(my_data), 1, Size, 1, stride(), Slicing());
	}

	Matrix<1, Size, Precision, Slice<1,Stride> > as_row(){
		return Matrix<1, Size, Precision, Slice<1,Stride> >(my_data, 1, Size, 1, stride(), Slicing());
	}

	const Matrix<Size, 1, Precision, Slice<Stride,1> > as_col() const{
		return Matrix<Size, 1, Precision, Slice<Stride,1> >(const_cast<Precision*>(my_data), Size, 1, stride(), 1, Slicing());
	}

	Matrix<Size, 1, Precision, Slice<Stride,1> > as_col(){
		return Matrix<Size, 1, Precision, Slice<Stride,1> >(my_data, Size, 1, stride(), 1, Slicing());
	}
};

}
