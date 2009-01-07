// -*- c++ -*-








template<int Rows, int Cols>
class MatrixAllocator{




};






template<int Rows, int Cols>
struct RowMajor{

  // has to handle operator[row] and operator(row,col)

  Vector<Cols,SDVBase<Cols,1> > operator[](int row){
    return Vector<Cols,SDVBase<Cols,1> >(my_data+Cols*row);
  }

  const Vector<Cols,SDVBase<Cols,1> > operator[](int row) const {
    return Vector<Cols,SDVBase<Cols,1> >(my_data+Cols*row);
  }

  double& operator()(int row, int col){
    return my_data[row*Cols+col];
  }
  const double& operator()(int row, int col) const {
    return my_data[row*Cols+col];
  }
};

template<int Rows, int Cols>
struct ColMajor {

  // has to handle operator[row] and operator(row,col)

};



template<int Stride>
struct StrideM {
  template <int Rows, int Cols>
  struct RowMajor{

  };
  template <int Rows, int Cols>
  struct ColMajor{


  };
};



template <int Stride>
struct MStride {
  template <int Rows, int Cols>
  struct RowMajor{
    // in here goes the maths to compute Matrix(r,c) etc.
  };
  struct ColMajor{
    // in here goes the maths to compute Matrix(r,c) etc.
  };
};


// can now make a matrix of form
// Matrix<Rows, Cols, MStride<Stride>::RowMajor>




template <int Rows, int Cols, template<int R, int C> Layout = RowMajor>
class Matrix : public Layout<Rows, Cols>;
