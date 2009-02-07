//Slicing and etc...
//SS = static size, static stride
//DS = dynamic size, static stride  (can use DD)
//Dl = dynamic size, linked stride: for RowMajor user instantiated dynamic matrices, stride==cols
//DD = dynamic size, dynamic stride
//SD = static size, dynamic stride.

//SS.slice<>() -> SS
//SS.slice()   -> DS
//DS.slice<>() -> SS
//DS.slice()   -> DS
//Dc.slice<>() -> SD
//Dl.silce()   -> DD
//DD.slice<>() -> SD
//DD.slice()   -> DD
//SD.slice<>() -> SD
//SD.slice()   -> DD


// As usual, a positive integer means static and -1 means dynamic.
// The new case is that for strides, -2 means that the stride is 
// the same as num_cols/num_rows, which must be dynamically sized.

template<int,int,class,template<int,int,class> class> class Matrix;
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericRowMajor;

//Closure used to acquire strides
//-1 means dynamic stride
//-2 means dynamic stride is tied to size
template<int Stride> struct Slice
{
	template<int Rows, int Cols, class Precision> struct RowMajor: public GenericRowMajor<Rows, Cols, Precision, Stride, MatrixSlice<Rows, Cols, Precision> >
	{
		//Optional constructors.
		
		RowMajor(Precision* p)
		:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p)
		{
		}

		RowMajor(Precision* p, int stride)
		:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, stride)
		{
		}

		RowMajor(Precision* p, int rows, int cols)
		:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols)
		{
		}

		RowMajor(Precision* p, int rows, int cols, int stride)
		:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols, stride)
		{
		}

	};
};


template<int Rows, int Cols, class Precision> struct RowMajor: public GenericRowMajor<Rows, Cols, Precision, (Cols==-1?-2:Cols), MatrixAlloc<Rows, Cols, Precision> >
{
	//Optional constructors.
	
	RowMajor(){}

	RowMajor(int rows, int cols)
	:GenericRowMajor<Rows, Cols, Precision, (Cols == -1 ? -2 : Cols), MatrixAlloc<Rows, Cols, Precision> >(rows, cols)
	{}
};


////////////////////////////////////////////////////////////////////////////////
//
// Given the rows, cols, and stride, what should the resulting vector type be?
//

template<int Rows, int Cols, typename Precision> struct RMVec
{
	typedef Vector<Cols, Precision, SVBase<Cols, 1, Precision> > Type;

	static Type vec(Precision* d, int /*cols*/)
	{
		return Type(d);
	}
};

template<typename Precision> struct RMVec<-1, -1, Precision>
{
	typedef Vector<-1, Precision, SDVBase<1, Precision> > Type;

	static Type vec(Precision* d, int cols)
	{
		return Type(d, cols);
	}
};

////////////////////////////////////////////////////////////////////////////////
//
// A class similar to mem, but to hold the stride information. It is only needed
// for -1. For +int and -2, the stride is part fo teh type, or implicit.

template<int s> struct StrideHolder
{
	//Constructos ignore superfluous arguments
	StrideHolder(){}
	StrideHolder(int){}

	int get(const void*) const{
		return s;
	}
};

template<> struct StrideHolder<-1>
{
	StrideHolder(int s)
	:stride(s){}

	const int stride;
	int get(const void*) const {
		return stride;
	}
};

template<> struct StrideHolder<-2>{
	template<class C> int get(const C* c) const {
		return c->tied_stride();
	}
};

////////////////////////////////////////////////////////////////////////////////
//
//Generic access to Row Major data.
//
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericRowMajor: public Mem
{
	//Slices can never have tied strides
	static const int SliceStride = Stride == -2?-1: Stride;

	//This little hack returns the stride value if it exists,
	//or one of the implied strides if they exist.
	int tied_stride() const{ 
		//Only valid if stride is -2
		return num_cols();
	}
	StrideHolder<Stride> my_stride;
	int stride() const {
		return my_stride.get(this);
	}


	//Optional constructors
	
	GenericRowMajor(){}

	GenericRowMajor(Precision* p)
	:Mem(p) {}

	GenericRowMajor(Precision* p, int s)
	:Mem(p),my_stride(s) {}

	GenericRowMajor(Precision* p, int r, int c)
	:Mem(p, r, c) {}

	GenericRowMajor(Precision* p, int r, int c, int stride)
	:Mem(p, r, c),my_stride(stride) {}

	GenericRowMajor(int r, int c)
	:Mem(r, c) {}

	using Mem::my_data;
	using Mem::num_cols;
	using Mem::num_rows;

	Precision& operator()(int r, int c){
		return my_data[r*stride() + c];
	}

	const Precision& operator()(int r, int c) const {
		return my_data[r*stride() + c];
	}
	
	typedef RMVec<Rows, Cols, Precision> Vec;
	typename Vec::Type operator[](int r)
	{
		return Vec::vec(my_data + stride()* r, num_cols());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, Slice<SliceStride>::template RowMajor> slice()
	{
		//Always pass the stride as a run-time parameter. It will be ignored
		//by SliceHolder (above) if it is statically determined.
		return Matrix<Rlength, Clength, Precision, Slice<SliceStride>::template RowMajor>(my_data+stride()*Rstart + Cstart, stride(), Slicing());
	}

	Matrix<-1, -1, Precision, Slice<SliceStride>::template RowMajor > slice(int rs, int cs, int rl, int cl){
		return Matrix<-1, -1, Precision, Slice<SliceStride>::template RowMajor >(my_data+stride()*rs +cs, rl, cl, stride(), Slicing());
	}
};
