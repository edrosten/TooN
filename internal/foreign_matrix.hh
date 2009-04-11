// -*- c++ -*-

////////////////////////////////////////////////////////////////////////////////
//
// Helper classes for matrices constructed as references to foreign data
//

namespace Reference
{

	struct RowMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			Layout(Precision* p, int r=0, int c=0)
				: Internal::GenericMBase<Rows,Cols,Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> > (p, r, c, 0, 0)
			{}

// 			template<class T> Layout(Precision* p, SliceSpec<T> spec)
// 			:Internal::GenericMBase<Rows, Cols, Precision, (Rows==-1?-2:Rows), 1, Internal::MatrixSlice<Rows, Cols, Precision> >(p, spec)
// 			{
// 			}
		};
	};

	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			Layout(Precision* p, int r=0, int c=0)
				: Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p, r, c, 0, 0)
			{}
		};
	};
}

template<int R, int C, typename Precision=double, class Type=Reference::RowMajor> struct Wrap
{
	static Matrix<R, C, Precision, Type> wrap(Precision* p)
	{
		return Matrix<R, C, Precision, Type>(p);
	}
};


template<int R, typename Precision, class Type> struct Wrap<R, Dynamic, Precision, Type>
{
	static Matrix<R, Dynamic, Precision, Type> wrap(Precision* p, int cols)
	{
		return Matrix<R, Dynamic, Precision, Type>(p, 0, cols);
	}
};


template<int C, typename Precision, class Type> struct Wrap<Dynamic, C, Precision, Type>
{
	static Matrix<Dynamic, C, Precision, Type> wrap(Precision* p, int rows)
	{
		return Matrix<Dynamic, C, Precision, Type>(p, rows, 0);
	}
};


template<typename Precision, class Type> struct Wrap<Dynamic, Dynamic, Precision, Type>
{
	static Matrix<Dynamic, Dynamic, Precision, Type> wrap(Precision* p, int rows, int cols)
	{
		return Matrix<Dynamic, Dynamic, Precision, Type>(p, rows, cols);
	}
};

