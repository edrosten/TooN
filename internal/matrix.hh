template <int Rows=-1, int Cols=Rows, class Precision=double, class Layout = RowMajor>
class Matrix : public Layout::template Layout<Rows, Cols, Precision>
{
  private:
	using Layout::template Layout<Rows, Cols, Precision>::my_data;
  public:
	//Use Tom's sneaky constructor hack...
		
	Matrix(){}

	Matrix(Precision* data, Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data){}
	
	//The stride is always passed during a slice. If it is not
	//needed, it will be ignored later and not stored.
	Matrix(Precision* data, int stride, Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data, stride){}

	Matrix(Precision* data, int rows, int cols, int stride, Slicing)
	:Layout::template Layout<Rows, Cols, Precision>(data, rows, cols, stride){}

	Matrix(int rows, int cols)
	:Layout::template Layout<Rows,Cols,Precision>(rows, cols)
	{}

	Precision* data() {
	  return my_data;
	}

	const Precision* data() const {
	  return my_data;
	}
};
