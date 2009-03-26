namespace Internal
{
// As usual, a positive integer means static and -1 means dynamic.
// The new case is that for strides, -2 means that the stride is 
// the same as num_cols/num_rows, which must be dynamically sized.

template<int, int, class, int, int, class> class GenericMBase;

////////////////////////////////////////////////////////////////////////////////
//Closure used to acquire strides
//-1 means dynamic stride
//-2 means dynamic stride is tied to size for a normal matrix
template<int RowStride, int ColStride> struct Slice
{
  
	template<int Rows, int Cols, class Precision> struct Layout: public GenericMBase<Rows, Cols, Precision, RowStride, ColStride, MatrixSlice<Rows, Cols, Precision> >
	{
		//Optional constructors.
		Layout(Precision* p, int rowstride, int colstride)
			:GenericMBase<Rows,Cols,Precision,RowStride,ColStride,MatrixSlice<Rows, Cols, Precision> >(p, rowstride, colstride, Internal::SpecifyStride())
		{
		}
		
		Layout(Precision* p, int rows, int cols, int rowstride, int colstride)
			:GenericMBase<Rows,Cols,Precision,RowStride,ColStride,MatrixSlice<Rows, Cols, Precision> >(p, rows, cols, rowstride, colstride)
		{
		}
	};
};


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
			return num_rows();
		} else {
			return ColStrideHolder<ColStride>::stride();
		}
	}

	//Optional constructors
	GenericMBase(){}

	GenericMBase(Precision* p)
	:Mem(p)
	{}

	//Both size and stride require two integer parameters.
	//In the case where only one is required, the meaning is 
	//ambiguous. Therefore a tag is used to specify the meaning.
	GenericMBase(Precision* p, int r, int c, Internal::SpecifySize)
	:Mem(p, r, c)
	{}

	GenericMBase(Precision* p, int rs, int cs, Internal::SpecifyStride)
	:Mem(p),RowStrideHolder<RowStride>(rs),ColStrideHolder<ColStride>(cs) {}


	GenericMBase(Precision* p, int r, int c, int rowstride, int colstride)
	:Mem(p, r, c),
	 RowStrideHolder<RowStride>(rowstride),
	 ColStrideHolder<ColStride>(colstride) 
	{}

	GenericMBase(int r, int c)
	:Mem(r, c) {}

	//Unpack slice specifications in to specific constructors
	//These slice specifications are not allowed to ignore superfluous data.
	GenericMBase(Precision* p, const Spec_____& s):Mem(p                    )                                                                                                              {}
	GenericMBase(Precision* p, const Spec____C& s):Mem(p                    )                                                        ,ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec___R_& s):Mem(p                    ), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore())                                                      {}
	GenericMBase(Precision* p, const Spec___RC& s):Mem(p                    ), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore()),ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec__C__& s):Mem(p, s.c ,SpecifyCols())                                                                                                              {}
	GenericMBase(Precision* p, const Spec__C_C& s):Mem(p, s.c, SpecifyCols())                                                        ,ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec__CR_& s):Mem(p, s.c, SpecifyCols()), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore())                                                      {}
	GenericMBase(Precision* p, const Spec__CRC& s):Mem(p, s.c, SpecifyCols()), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore()),ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec_R___& s):Mem(p, s.r, SpecifyRows())                                                                                                              {}
	GenericMBase(Precision* p, const Spec_R__C& s):Mem(p, s.r, SpecifyRows())                                                        ,ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec_R_R_& s):Mem(p, s.r, SpecifyRows()), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore())                                                      {}
	GenericMBase(Precision* p, const Spec_R_RC& s):Mem(p, s.r, SpecifyRows()), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore()),ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec_RC__& s):Mem(p, s.r, s.c          )                                                                                                              {}
	GenericMBase(Precision* p, const Spec_RC_C& s):Mem(p, s.r, s.c          )                                                        ,ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}
	GenericMBase(Precision* p, const Spec_RCR_& s):Mem(p, s.r, s.c          ), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore())                                                      {}
	GenericMBase(Precision* p, const Spec_RCRC& s):Mem(p, s.r, s.c          ), RowStrideHolder<RowStride>(s.rs, Internal::NoIgnore()),ColStrideHolder<ColStride>(s.cs,Internal::NoIgnore()){}

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
	Matrix<Rlength, Clength, Precision, Slice<SliceRowStride,SliceColStride> > slice()
	{
		//Always pass the stride as a run-time parameter. It will be ignored
		//by SliceHolder (above) if it is statically determined.
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, Slice<SliceRowStride,SliceColStride> >(my_data+rowstride()*Rstart + colstride()*Cstart, rowstride(), colstride(), Slicing());
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	const Matrix<Rlength, Clength, Precision, Slice<SliceRowStride,SliceColStride> > slice() const
	{
		Internal::CheckStaticSlice<Rows, Rstart, Rlength>::check(num_rows());
		Internal::CheckStaticSlice<Cols, Cstart, Clength>::check(num_cols());
		return Matrix<Rlength, Clength, Precision, Slice<SliceRowStride,SliceColStride> >(const_cast<Precision*>(my_data+rowstride()*Rstart + colstride()*Cstart), rowstride(), colstride(), Slicing());
	}

	Matrix<-1, -1, Precision, Slice<SliceRowStride,SliceColStride> > slice(int rs, int cs, int rl, int cl){
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, Slice<SliceRowStride,SliceColStride> >(my_data+rowstride()*rs +colstride()*cs, rl, cl, rowstride(), colstride(), Slicing());
	}

	const Matrix<-1, -1, Precision, Slice<SliceRowStride,SliceColStride> > slice(int rs, int cs, int rl, int cl) const {
		Internal::CheckDynamicSlice::check(num_rows(), rs, rl);
		Internal::CheckDynamicSlice::check(num_cols(), cs, cl);
		return Matrix<-1, -1, Precision, Slice<SliceRowStride,SliceColStride> >(const_cast<Precision*>(my_data+rowstride()*rs +colstride()*cs), rl, cl, rowstride(), colstride(), Slicing());
	}


	Matrix<Cols, Rows, Precision, Slice<SliceColStride,SliceRowStride> > T(){
		return Matrix<Cols, Rows, Precision, Slice<SliceColStride,SliceRowStride> >(my_data, num_cols(), num_rows(), colstride(), rowstride(), Slicing());
	}

	const Matrix<Cols, Rows, Precision, Slice<SliceColStride,SliceRowStride> > T() const{
		return Matrix<Cols, Rows, Precision, Slice<SliceColStride,SliceRowStride> >(const_cast<Precision*>(my_data), num_cols(), num_rows(), colstride(), rowstride(), Slicing());
	}
};

}




////////////////////////////////////////////////////////////////////////////////
//
// Classes for Matrices owning memory
//
//
struct RowMajor
{
	template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixAlloc<Rows, Cols, Precision> >
	{
		//Optional constructors.
		
		Layout(){}

		Layout(int rows, int cols)
		:Internal::GenericMBase<Rows, Cols, Precision, (Cols == -1 ? -2 : Cols), 1, Internal::MatrixAlloc<Rows, Cols, Precision> >(rows, cols)
		{}
	};
};

struct ColMajor
{
	template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixAlloc<Rows, Cols, Precision> >
	{
		//Optional constructors.
		
		Layout(){}

		Layout(int rows, int cols)
		:Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows == -1 ? -2 : Rows), Internal::MatrixAlloc<Rows, Cols, Precision> >(rows, cols)
		{}
	};
};

////////////////////////////////////////////////////////////////////////////////
//
// Helper classes for matrices constructed as references to foreign data
//

namespace Reference
{

	struct RowMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, (Rows==-1?-2:Rows), 1, Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			template<class T> Layout(Precision* p, SliceSpec<T> spec)
			:Internal::GenericMBase<Rows, Cols, Precision, (Rows==-1?-2:Rows), 1, Internal::MatrixSlice<Rows, Cols, Precision> >(p, spec)
			{
			}
		};
	};
/*
	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			//Size fixed. Only data required.
			Layout(Precision* p)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p)
			{
			}

			//Size dynamic: size and data needed.
			Layout(Precision* p, int rows, int cols)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows == -1 ? -2 : Rows), Internal::MatrixAlloc<Rows, Cols, Precision> >(p, rows, cols, Internal::SpecifySize())
			{}
		};
	};


	template<int Stride=-1> struct ColMajorStride
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, 1, Stride, Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			
			//Size and stride fixed. Only data required.
			Layout(Precision* p)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, Stride, Internal::MatrixSlice<Rows, Cols, Precision> >(p)
			{
			}

			//Size fixed. Data and stride required.
			Layout(Precision* p, int stride)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, Stride, Internal::MatrixSlice<Rows, Cols, Precision> >(p, 1, stride, Internal::SpecifySize())
			{
			}

			//Stride fixed. Data and size required.
			Layout(Precision* p, int rows, int cols)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, Stride, Internal::MatrixAlloc<Rows, Cols, Precision> >(p, rows, cols, Internal::SpecifySize())
			{}

			//Dynamic: Stride,  data and size required.
			Layout(Precision* p, int rows, int cols, int stride)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, Stride, Internal::MatrixAlloc<Rows, Cols, Precision> >(p, rows, cols, 1, stride)
			{}
		};
	};*/
}





