template <int Rows=-1, int Cols=Rows, class Precision=double, class Layout = RowMajor>
class Matrix : public Layout::template Layout<Rows, Cols, Precision>
{
private:
	//using Layout::template Layout<Rows, Cols, Precision>::my_data;
public:

	using Layout::template Layout<Rows, Cols, Precision>::my_data;

	using Layout::template Layout<Rows, Cols, Precision>::num_rows;
	using Layout::template Layout<Rows, Cols, Precision>::num_cols;

	//Use Tom's sneaky constructor hack...
		
	Matrix(){}

	//The stride is always passed during a slice. If it is not
	//needed, it will be ignored later and not stored.
	Matrix(Precision* data, int rowstride, int colstride, Internal::Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data, rowstride, colstride){}

	Matrix(Precision* data, int rows, int cols, int rowstride, int colstride, Internal::Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data, rows, cols, rowstride, colstride){}

	Matrix(int rows, int cols)
	:Layout::template Layout<Rows,Cols,Precision>(rows, cols)
	{}

	
	//See vector.hh and allocator.hh for details about why the
	//copy constructor should be default.

	// constructors to allow return value optimisations
	// construction from 1-ary operator
	template <class T, class Op>
	inline Matrix(const T& arg, const Operator<Op>&, int rows, int cols) 
	:Layout::template Layout<Rows,Cols,Precision>(rows, cols) 
	{
	    Op::eval(*this,arg);
	}

	// constructor from 2-ary operator
	template <class LHS, class RHS, class Op>
	inline Matrix(const LHS& lhs, const RHS& rhs, const Operator<Op>&, int rows, int cols)
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
	  	  for(int c=0; c < num_rows(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	// operator =
	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_rows(); c++)
	  	  	(*this)[r][c] = from[r][c];

	    return *this;
	}

	Matrix& operator+=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_rows(); c++)
			  	(*this)[r][c] += rhs;

		  return *this;
	}

	Matrix& operator-=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_rows(); c++)
			  	(*this)[r][c] -= rhs;

		  return *this;
	}

	Matrix& operator*=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_rows(); c++)
			  	(*this)[r][c] *= rhs;

		  return *this;
	}

	Matrix& operator/=(const Precision& rhs)
	{
		  for(int r=0; r < num_rows(); r++)
			  for(int c=0; c < num_rows(); c++)
			  	(*this)[r][c] /= rhs;

		  return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator+= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_rows(); c++)
	  	  	(*this)[r][c] += from[r][c];

	    return *this;
	}

	template<int Rows2, int Cols2, typename Precision2, typename Base2>
	Matrix& operator-= (const Matrix<Rows2, Cols2, Precision2, Base2>& from)
	{
		SizeMismatch<Rows, Rows2>::test(num_rows(), from.num_rows());
		SizeMismatch<Cols, Cols2>::test(num_cols(), from.num_cols());

	    for(int r=0; r < num_rows(); r++)
	  	  for(int c=0; c < num_rows(); c++)
	  	  	(*this)[r][c] -= from[r][c];

	    return *this;
	}
};
