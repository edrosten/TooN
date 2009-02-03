//Types:
//SS = static size, static stride
//DS = dynamic size, static stride  (can use DD)
//Du = dynamic size, unstrided (used for user instantiated Matrix), saves implicit stride parameter
//DD = dynamic size, dynamic stride
//SD = static size, dynamic stride.

//SS.slice<>() -> SS
//SS.slice()   -> DS
//DS.slice<>() -> SS
//DS.slice()   -> DS
//Du.slice<>() -> SD
//Du.silce()   -> DD
//DD.slice<>() -> SD
//DD.slice()   -> DD
//SD.slice<>() -> SD
//SD.slice()   -> DD

template<int,int,class,template<int,int,class> class> class Matrix;
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct SSRowMajor;
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct SSColMajor;

template<int Stride> struct SSSlice
{
	template<int Rows, int Cols, class Precision> struct RowMajor: public SSRowMajor<Rows, Cols, Precision, Stride, SliceHolder<Precision> >
	{
		RowMajor(Precision* p,  Slicing)
		:SSRowMajor<Rows,Cols,Precision,Stride,SliceHolder<Precision> >(p)
		{
		}
	};

	template<int Rows, int Cols, class Precision> struct ColMajor: public SSColMajor<Rows, Cols, Precision, Stride, SliceHolder<Precision> >
	{
		ColMajor(Precision* p, Slicing)
		:SSColMajor<Rows,Cols,Precision,Stride,SliceHolder<Precision> >(p)
		{
		}
	};
};

template<int Rows, int Cols, class Precision> struct RowMajor: public SSRowMajor<Rows, Cols, Precision, Cols, StaticSizedAllocator<Rows*Cols, Precision> >
{
};

template<int Rows, int Cols, class Precision> struct ColMajor: public SSColMajor<Rows, Cols, Precision, Rows, StaticSizedAllocator<Rows*Cols, Precision> >
{
};


//Generic access to Static sized, statically strided Row Major data.
template<int Rows, int Cols, class Precision, int Stride, class Mem> struct SSRowMajor: public Mem
{
	//Optional constructors
	
	SSRowMajor(){}
	SSRowMajor(Precision* p)
	:Mem(p)
	{}

	int num_rows()const {
		return Rows;
	}
	int num_cols()const{
		return Cols;
	}
	using Mem::my_data;

	Precision& operator()(int r, int c){
		Internal::check_index(Rows, r);	
		Internal::check_index(Cols, c);	
		return my_data[r*Stride + c];
	}

	const Precision& operator()(int r, int c) const {
		Internal::check_index(Rows, r);	
		Internal::check_index(Cols, c);	
		return my_data[r*Stride + c];
	}

	Vector<Cols, Precision, SVBase<Cols, 1, Precision> > operator[](int r)
	{
		Internal::check_index(Rows, r);	
		return Vector<Cols, Precision, SVBase<Cols, 1, Precision> >(my_data + Stride* r);
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template RowMajor> slice()
	{
		Internal::CheckSlice<Rows, Rstart, Rlength>::check();
		Internal::CheckSlice<Cols, Cstart, Clength>::check();

		return Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template RowMajor>(my_data+Stride*Rstart + Cstart, Slicing());
	}

	Matrix<Cols, Rows, Precision, SSSlice<Stride>::template ColMajor> T()
	{
		return Matrix<Cols, Rows, Precision, SSSlice<Stride>::template ColMajor>(my_data, Slicing());
	}
};

template<int Rows, int Cols, class Precision, int Stride, class Mem> struct SSColMajor: public Mem
{
	SSColMajor(){}
	SSColMajor(Precision* p)
	:Mem(p)
	{}

	int num_rows()const{
		return Rows;
	}
	int num_cols()const{
		return Cols;
	}

	using Mem::my_data;

	Precision& operator()(int r, int c){
		Internal::check_index(Rows, r);	
		Internal::check_index(Cols, c);	
		return my_data[c*Stride + r];
	}

	const Precision& operator()(int r, int c)const{
		Internal::check_index(Rows, r);	
		Internal::check_index(Cols, c);	
		return my_data[c*Stride + r];
	}

	Vector<Cols, Precision, SVBase<Cols, Stride, Precision> > operator[](int r)
	{
		Internal::check_index(Rows, r);	
		return Vector<Cols, Precision, SVBase<Cols, Stride, Precision> >(my_data + r);
	}

	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template ColMajor> slice()
	{
		Internal::CheckSlice<Rows, Rstart, Rlength>::check();
		Internal::CheckSlice<Cols, Cstart, Clength>::check();

		return Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template ColMajor>(my_data+Stride*Cstart + Rstart, Slicing());
	}

	Matrix<Cols, Rows, Precision, SSSlice<Stride>::template RowMajor> T()
	{
		return Matrix<Cols, Rows, Precision, SSSlice<Stride>::template RowMajor>(my_data, Slicing());
	}
};












