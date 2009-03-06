// As usual, a positive integer means static and -1 means dynamic.
// The new case is that for strides, -2 means that the stride is 
// the same as num_cols/num_rows, which must be dynamically sized.

template<int,int,class,class> class Matrix;

////////////////////////////////////////////////////////////////////////////////
//Closure used to acquire strides
//-1 means dynamic stride
//-2 means dynamic stride is tied to size for a normal matrix
template<int RowStride, int ColStride> struct Slice
{
	struct Layout{
		template<int Rows, int Cols, class Precision> struct Layout: public GenericMBase<Rows, Cols, Precision, RowStride, ColStride, MatrixSlice<Rows, Cols, Precision> >
		{
			//Optional constructors.
			
			Layout(Precision* p)
				:GenericMBase<Rows,Cols,Precision,RowStride,ColStride,MatrixSlice<Rows, Cols, Precision> >(p)
			{
			}
			
#if 0
			Layout(Precision* p, int stride)
				:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, stride)
			{
			}
			
			Layout(Precision* p, int rows, int cols)
				:GenericRowMajor<Rows,Cols,Precision,Stride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols)
			{
			}
#endif		
			
			
			Layout(Precision* p, int rows, int cols, int rowstride, int colstride)
				:GenericMBase<Rows,Cols,Precision,RowStride,ColStride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols, rowstride, colstride)
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

template<int Rows, int Cols, class Precision, int RowStride, int ColStride, class Mem> struct GenericMBase
	: public Mem, 
	RowStrideHolder<RowStride>,
	ColStrideHolder<ColStride>
{
	//Slices can never have tied strides
	static const int SliceRowStride = RowStride == -2?-1: RowStride;
	static const int SliceColStride = ColStride == -2?-1: ColStride;

	int rowstride() const {
		if(RowStride == -2) { //Normal tied stride
			return num_cols();
		} else {
			return RowStrideHolder<RowStride>::stride();
		}
	}

	int colstride() const {
		if(ColStride == -2) { //Normal tied stride
			return num_row();
		} else {
			return ColStrideHolder<ColStride>::stride();
		}
	}

	//Optional constructors
	GenericMBaser(){}

	GenericMBase(Precision* p)
	:Mem(p) {}

	// We think we don't need these any more
#if 0
	GenericRowMajor(Precision* p, int s)
	:Mem(p),StrideHolder<Stride>(s) {}

	GenericRowMajor(Precision* p, int r, int c)
	:Mem(p, r, c) {}
#endif


	GenericMBase(Precision* p, int r, int c, int rowstride, int colstride)
	:Mem(p, r, c),
	 RowStrideHolder<RowStride>(rowstride) 
	 ColStrideHolder<ColStride>(colstride) 
	{}

	GenericMBase(int r, int c)
	:Mem(r, c) {}

	using Mem::my_data;
	using Mem::num_cols;
	using Mem::num_rows;

	Precision& operator()(int r, int c){
		return my_data[r*rowstride() + c*colstride()];
	}

	const Precision& operator()(int r, int c) const {
		return my_data[r*rowstride() + c*colstride()];
	}

	// this is the type of vector obtained by [ ]
	typedef Vector<Cols, Precision, SliceVBase<SliceColStride> > Vec;
	
	Vec operator[](int r) {
		Internal::check_index(num_rows(), r);
		return Vec(my_data + rowstride()* r, num_cols(), colstride(), Slicing());
	}

	const Vec operator[](int r) const {
		Internal::check_index(num_rows(), r);
		return Vec(const_cast<Precision*>(my_data + rowstride()* r), num_cols(), colstride(), Slicing());
	}



	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, typename Slice<SliceRowStride,SliceColStride> > slice()
	{
		//Always pass the stride as a run-time parameter. It will be ignored
		//by SliceHolder (above) if it is statically determined.
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceRowStride,SliceColStride>::Layout>(my_data+rowstride()*Rstart + colstride()*Cstart, rowstride(), colstride(), Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	const Matrix<Rlength, Clength, Precision, typename Slice<SliceRowStride,SliceColStride>::Layout> slice() const
	{
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, typename Slice<SliceRowStride,SliceColStride>::Layout>(const_cast<Precision*>(my_data+rowstride()*Rstart + colstride()*Cstart), rowstride(), colstride(), Slicing());
	}

	Matrix<-1, -1, Precision, typename Slice<SliceRowStride,SliceColStride>::RowMajor > slice(int rs, int cs, int rl, int cl){
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceRowStride,SliceColStride>::RowMajor >(my_data+rowstride()*rs +colstride()*cs, rl, cl, rowstride(), colstride(), Slicing());
	}

	const Matrix<-1, -1, Precision, typename Slice<SliceRowStride,SliceColStride>::RowMajor > slice(int rs, int cs, int rl, int cl) const {
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, typename Slice<SliceRowStride,SliceColStride>::RowMajor >(const_cast<Precision*>(my_data+stride()*rs +colstride()*cs), rl, cl, rowstride(), colstride(), Slicing());
	}


	Matrix<Cols, Rows, Precision, typename Slice<SliceColStride,SliceRowStride>::Layout> T(){
		return Matrix<Cols, Rows, Precision, typename Slice<SliceColStride,SliceRowStride>::Layout>(my_data, num_cols(), num_rows(), colstride(), rowstride(), Slicing());
	}

	const Matrix<Cols, Rows, Precision, typename Slice<SliceColStride,SliceRowStride>::Layout> T() const{
		return Matrix<Cols, Rows, Precision, typename Slice<SliceColStride,SliceRowStride>::Layout>(const_cast<Precision*>(my_data), num_cols(), num_rows(), colstride(), rowstride(), Slicing());
	}
};






