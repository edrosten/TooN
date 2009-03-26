
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

	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct Layout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			template<class T> Layout(Precision* p, SliceSpec<T> spec)
			:Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p, spec)
			{
			}
		};
	};
}


template<int R, int C, typename Precision=double, class Type=Reference::RowMajor> struct Wrap
{
	static Matrix<R, C, Precision, Type> wrap(Precision* p)
	{
		SliceSpec<Internal::Spec_____> s;
		return Matrix<R, C, Precision, Type>(p, s);
	}
};


template<int R, typename Precision, class Type> struct Wrap<R, Dynamic, Precision, Type>
{
	static Matrix<R, Dynamic, Precision, Type> wrap(Precision* p, int cols)
	{
		SliceSpec<Internal::Spec__C__> s;
		s.c = cols;
		return Matrix<R, Dynamic, Precision, Type>(p, s);
	}
};


template<int C, typename Precision, class Type> struct Wrap<Dynamic, C, Precision, Type>
{
	static Matrix<Dynamic, C, Precision, Type> wrap(Precision* p, int rows)
	{
		SliceSpec<Internal::Spec_R___> s;
		s.r = rows;
		return Matrix<Dynamic, C, Precision, Type>(p, s);
	}
};


template<typename Precision, class Type> struct Wrap<Dynamic, Dynamic, Precision, Type>
{
	static Matrix<Dynamic, Dynamic, Precision, Type> wrap(Precision* p, int rows, int cols)
	{
		SliceSpec<Internal::Spec_RC__> s;
		s.r = rows;
		s.c = cols;
		return Matrix<Dynamic, Dynamic, Precision, Type>(p, s);
	}
};

