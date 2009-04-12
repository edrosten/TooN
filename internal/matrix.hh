template <int Rows=-1, int Cols=Rows, class Precision=double, class Layout = RowMajor>
class Matrix : public Layout::template Layout<Rows, Cols, Precision>
{
public:

	using Layout::template Layout<Rows, Cols, Precision>::my_data;

	using Layout::template Layout<Rows, Cols, Precision>::num_rows;
	using Layout::template Layout<Rows, Cols, Precision>::num_cols;

	//Use Tom's sneaky constructor hack...
	Matrix(){}

	Matrix(int rows, int cols) :
		Layout::template Layout<Rows,Cols,Precision>(rows, cols)
	{}

	Matrix(Precision* p) :
		Layout::template Layout<Rows, Cols, Precision>(p)
	{}

	Matrix(Precision* p, int r, int c) :
		Layout::template Layout<Rows, Cols, Precision>(p, r, c)
	{}

	// Internal constructor used by GenericMBase::slice(...)
	Matrix(Precision* data, int rows, int cols, int rowstride, int colstride, Internal::Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data, rows, cols, rowstride, colstride){}

	//See vector.hh and allocator.hh for details about why the
	//copy constructor should be default.
	template <class Op>
	inline Matrix(const Operator<Op>& op)
		:Layout::template Layout<Rows,Cols,Precision>(op)
	{
		op.eval(*this);
	}

	// constructors to allow return value optimisations
	// construction from 1-ary operator
	template <class T, class Op>
	inline Matrix(const T& arg, int rows, int cols, const Operator<Op>&) 
	:Layout::template Layout<Rows,Cols,Precision>(rows, cols) 
	{
	    Op::eval(*this,arg);
	}

	// constructor from 2-ary operator
	template <class LHS, class RHS, class Op>
	inline Matrix(const LHS& lhs, const RHS& rhs, int rows, int cols, const Operator<Op>&)
	:Layout::template Layout<Rows,Cols,Precision>(rows, cols)
	{
	    Op::eval(*this,lhs,rhs);
	}

	// constructor from arbitrary matrix
	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	inline Matrix(const Matrix<Rows2, Cols2,Precision2,Base2>& from)
	:Layout::template Layout<Rows,Cols,Precision>(from.num_rows(), from.num_cols())
	{
	    operator=(from);
	}

	// operator = from copy
	inline Matrix& operator= (const Matrix& from)
	{
		SizeMismatch<Rows, Rows>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	// operator = 0-ary operator
	template<class Op> inline Matrix& operator= (const Operator<Op>& op)
	{
		op.eval(*this);
		return *this;
	}

	// operator =
	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	Matrix& operator+=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] += rhs;

		  return *this;
	}

	Matrix& operator-=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] -= rhs;

		  return *this;
	}

	Matrix& operator*=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] *= rhs;

		  return *this;
	}

	Matrix& operator/=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_cols(); c++)
			  	(*this)[r][c] /= rhs;

		  return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator+= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] += from[r][c];

	    return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator-= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_cols(); c++)
	  	  	(*this)[r][c] -= from[r][c];

	    return *this;
	}

	Matrix& ref()
	{
		return *this;
	}
};
