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

template<int,int,class,class> class Matrix;
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericRowMajor;
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericColMajor;


////////////////////////////////////////////////////////////////////////////////
//Closure used to acquire strides
//-1 means dynamic stride
//-2 means dynamic stride is tied to size for a normal matrix
//-3 means dynamic stride is tied to size for a transposed matrix
template<int Stride> struct Slice
{
	struct RowMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public GenericRowMajor<Rows, Cols, Precision, Stride, MatrixSlice<Rows, Cols, Precision> >
		{
			//Optional constructors.
			
			Layout(Precision* p)
			:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p)
			{
			}

			Layout(Precision* p, int stride)
			:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, stride)
			{
			}

			Layout(Precision* p, int rows, int cols)
			:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols)
			{
			}

			Layout(Precision* p, int rows, int cols, int stride)
			:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols, stride)
			{
			}

		};
	};

	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public GenericColMajor<Rows, Cols, Precision, Stride, MatrixSlice<Rows, Cols, Precision> >
		{
			//Optional constructors.
			
			Layout(Precision* p)
			:GenericColMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p)
			{
			}

			Layout(Precision* p, int stride)
			:GenericColMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, stride)
			{
			}

			Layout(Precision* p, int rows, int cols)
			:GenericColMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols)
			{
			}

			Layout(Precision* p, int rows, int cols, int stride)
			:GenericColMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols, stride)
			{
			}

		};
	};
};


////////////////////////////////////////////////////////////////////////////////
//
// Classes for Matrices owning memory
//
//
struct RowMajor
{
	template<int Rows, int Cols, class Precision> struct Layout: public GenericRowMajor<Rows, Cols, Precision, (Cols==-1?-2:Cols), MatrixAlloc<Rows, Cols, Precision> >
	{
		//Optional constructors.
		
		Layout(){}

		Layout(int rows, int cols)
		:GenericRowMajor<Rows, Cols, Precision, (Cols == -1 ? -2 : Cols), MatrixAlloc<Rows, Cols, Precision> >(rows, cols)
		{}
	};
};

struct ColMajor
{
	template<int Rows, int Cols, class Precision> struct Layout: public GenericColMajor<Rows, Cols, Precision, (Rows==-1?-2:Rows), MatrixAlloc<Rows, Cols, Precision> >
	{
		//Optional constructors.
		
		Layout(){}

		Layout(int rows, int cols)
		:GenericColMajor<Rows, Cols, Precision, (Rows == -1 ? -2 : Rows), MatrixAlloc<Rows, Cols, Precision> >(rows, cols)
		{}
	};
};




////////////////////////////////////////////////////////////////////////////////
//
// Row major matrix implementation
//

template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericRowMajor: public Mem, StrideHolder<Stride>
{
	//Slices can never have tied strides
	static const int SliceStride = Stride == -2?-1: Stride;

	int stride() const {
		if(Stride == -2) //Normal tied stride
			return num_cols();
		else
			return StrideHolder<Stride>::stride();
	}

	//Optional constructors
	GenericRowMajor(){}

	GenericRowMajor(Precision* p)
	:Mem(p) {}

	GenericRowMajor(Precision* p, int s)
	:Mem(p),StrideHolder<Stride>(s) {}

	GenericRowMajor(Precision* p, int r, int c)
	:Mem(p, r, c) {}

	GenericRowMajor(Precision* p, int r, int c, int stride)
	:Mem(p, r, c),StrideHolder<Stride>(stride) {}

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


	typedef Vector<Cols, Precision, SliceVBase<1> > Vec;
	
	Vec operator[](int r) {
		Internal::check_index(num_rows(), r);
		return Vec(my_data + stride()* r, num_cols(), 1, Slicing());
	}

	const Vec operator[](int r) const {
		Internal::check_index(num_rows(), r);
		return Vec(const_cast<Precision*>(my_data + stride()* r), num_cols(), 1, Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::RowMajor> slice()
	{
		//Always pass the stride as a run-time parameter. It will be ignored
		//by SliceHolder (above) if it is statically determined.
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::RowMajor>(my_data+stride()*Rstart + Cstart, stride(), Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	const Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::RowMajor> slice() const
	{
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::RowMajor>(const_cast<Precision*>(my_data+stride()*Rstart + Cstart), stride(), Slicing());
	}

	Matrix<-1, -1, Precision, typename Slice<SliceStride>::RowMajor > slice(int rs, int cs, int rl, int cl){
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceStride>::RowMajor >(my_data+stride()*rs +cs, rl, cl, stride(), Slicing());
	}

	const Matrix<-1, -1, Precision, typename Slice<SliceStride>::RowMajor > slice(int rs, int cs, int rl, int cl) const {
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceStride>::RowMajor >(const_cast<Precision*>(my_data+stride()*rs +cs), rl, cl, stride(), Slicing());
	}


	Matrix<Cols, Rows, Precision, typename Slice<Stride>::ColMajor> T(){
		return Matrix<Cols, Rows, Precision, typename Slice<Stride>::ColMajor>(my_data, num_cols(), num_rows(), stride(), Slicing());
	}

	const Matrix<Cols, Rows, Precision, typename Slice<Stride>::ColMajor> T() const{
		return Matrix<Cols, Rows, Precision, typename Slice<Stride>::ColMajor>(const_cast<Precision*>(my_data), num_cols(), num_rows(), stride(), Slicing());
	}
};






////////////////////////////////////////////////////////////////////////////////
//
// Column major matrix implementation
//

template<int Rows, int Cols, class Precision, int Stride, class Mem> struct GenericColMajor: public Mem, public StrideHolder<Stride>
{
	//Slices can never have tied strides
	static const int SliceStride = Stride == -2?-1: Stride;

	int stride() const {
		if(Stride == -2) //Normal tied stride
			return num_rows();
		else
			return StrideHolder<Stride>::stride();
	}


	//Optional constructors
	
	GenericColMajor(){}

	GenericColMajor(Precision* p)
	:Mem(p) {}

	GenericColMajor(Precision* p, int s)
	:Mem(p),StrideHolder<Stride>(s) {}

	GenericColMajor(Precision* p, int r, int c)
	:Mem(p, r, c) {}

	GenericColMajor(Precision* p, int r, int c, int stride)
	:Mem(p, r, c),StrideHolder<Stride>(stride) {}

	GenericColMajor(int r, int c)
	:Mem(r, c) {}

	using Mem::my_data;
	using Mem::num_cols;
	using Mem::num_rows;

	Precision& operator()(int r, int c){
		return my_data[c*stride() + r];
	}

	const Precision& operator()(int r, int c) const {
		return my_data[c*stride() + r];
	}
	
	typedef Vector<Cols, Precision, SliceVBase<SliceStride> > Vec;
	Vec operator[](int r) {
		Internal::check_index(num_rows(), r);
		return Vec(my_data + r, num_cols(), stride(), Slicing());
	}

	const Vec operator[](int r) const {
		Internal::check_index(num_rows(), r);
		return Vec(const_cast<Precision*>(my_data + r), num_cols(), stride(), Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::ColMajor> slice()
	{
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		//Always pass the stride as a run-time parameter. It will be ignored
		//by SliceHolder (above) if it is statically determined.
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::ColMajor>(my_data+Rstart + stride()*Cstart, stride(), Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	const Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::ColMajor> slice() const
	{
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceStride>::ColMajor>(const_cast<Precision*>(my_data+Rstart + stride()*Cstart), stride(), Slicing());
	}
	
	Matrix<-1, -1, Precision, typename Slice<SliceStride>::ColMajor > slice(int rs, int cs, int rl, int cl){
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceStride>::ColMajor >(my_data+rs +stride()*cs, rl, cl, stride(), Slicing());
	}

	const Matrix<-1, -1, Precision, typename Slice<SliceStride>::ColMajor > slice(int rs, int cs, int rl, int cl)const{
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceStride>::ColMajor >(const_cast<Precision*>(my_data+rs +stride()*cs), rl, cl, stride(), Slicing());
	}

	Matrix<Cols, Rows, Precision, typename Slice<Stride>::RowMajor> T(){
		return Matrix<Cols, Rows, Precision, typename Slice<Stride>::RowMajor>(my_data, num_cols(), num_rows(), stride(), Slicing());
	}

	const Matrix<Cols, Rows, Precision, typename Slice<Stride>::RowMajor> T() const {
		return Matrix<Cols, Rows, Precision, typename Slice<Stride>::RowMajor>(const_cast<Precision*>(my_data), num_cols(), num_rows(), stride(), Slicing());
	}
};



