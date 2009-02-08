template<int,class,class> class Vector;
template<int Size, class Precision, int Stride, class Mem> struct GenericVBase;
template<typename T> class Operator;


////////////////////////////////////////////////////////////////////////////////
//
// A class similar to mem, but to hold the stride information. It is only needed
// for -1. For +int and -2, the stride is part fo teh type, or implicit.

template<int s> struct VStrideHolder
{
	//Constructos ignore superfluous arguments
	VStrideHolder(){}
	VStrideHolder(int){}

	int stride() const{
		return s;
	}
};

template<> struct VStrideHolder<-1>
{
	VStrideHolder(int s)
	:my_stride(s){}

	const int my_stride;
	int stride() const {
		return my_stride;
	}
};

////////////////////////////////////////////////////////////////////////////////
//
// Slice holding class
//

template<int Size, int Stride, class Precision> struct SliceVBase: public GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >{
	

	SliceVBase(Precision* d, int stride)
	:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(d, stride){
	}

	SliceVBase(Precision* d, int length, int stride)
	:GenericVBase<Size, Precision, Stride, VectorSlice<Size, Precision> >(d, length, stride){
	}

};

////////////////////////////////////////////////////////////////////////////////
//
// Classes for Matrices owning memory
//

template<int Size, class Precision> struct VBase : public GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >{
	//Optional 
	
	VBase(){}

	VBase(int s)
	:GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> >(s)
	{
	}
};

////////////////////////////////////////////////////////////////////////////////
//
// Generic implementation
//

template<int Size, typename Precision, int Stride, typename Mem> struct GenericVBase: public Mem
{	
	VStrideHolder<Stride> my_stride;
	int stride() const{
		return my_stride.stride();
	}

	//Optional constuctors
	GenericVBase(){}

	GenericVBase(int s)
	:Mem(s)
	{}

	GenericVBase(Precision* d, int stride)
	:Mem(d),my_stride(stride){
	}

	GenericVBase(Precision* d, int length, int stride)
	:Mem(d, length),my_stride(stride){
	}

	using Mem::my_data;

	Precision& operator[](int i) {
		return my_data[i * stride()];
	}

	const Precision& operator[](int i) const {
		return my_data[i * stride()];
	}

	template<int Start, int Length> 
	Vector<Length, Precision, SliceVBase<Length, Stride, Precision> > slice(){
		return Vector<Length, Precision, SliceVBase<Length, Stride, Precision> >(my_data + stride()*Start, stride(), Slicing());
	}

	Vector<-1, Precision, SliceVBase<-1, Stride, Precision> > slice(int start, int length){
		return Vector<-1, Precision, SliceVBase<-1, Stride, Precision> >(my_data + stride()*start, length, stride(), Slicing());
	}
};
