// -*- c++ -*-

//Forward declarations


template<int Rows, int Cols, class Precision, template<int, int, class> class Layout> class Matrix;

struct Slicing{};

//Static stride Static slice 
template<int Stride> struct SSSlice
{
  template<int Rows, int Cols, class Precision>
  struct RowMajor
  {
    private: 
      Precision* my_data;

    public:
      int num_rows() const { return Rows; }
      int num_cols() const { return Cols; }
      //Construction

      RowMajor(Precision* d)
      :my_data(d)
      {}

      //Indexing	

      Vector<Cols, Precision, SVBase<Cols,1,Precision> > operator[](int row){
        Internal::check_index(Rows, row);
        return Vector<Cols, Precision, SVBase<Cols,1,Precision> >(my_data+Stride*row);
      }

      Precision& operator()(int row, int col){
        Internal::check_index(Rows, row);
        Internal::check_index(Cols, col);
        return my_data[row*Stride+col];
      }

      const Precision& operator()(int row, int col) const {
        Internal::check_index(Rows, row);
        Internal::check_index(Cols, col);
        return my_data[row*Stride+col];
      }

      //Slices go here
      template<int Rstart, int Cstart, int Rlength, int Clength>
      Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template RowMajor> slice(){
        Internal::CheckSlice<Rows, Rstart, Rlength>::check();
        Internal::CheckSlice<Cols, Cstart, Clength>::check();

        return Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template RowMajor>(my_data + Stride*Rstart + Cstart, Slicing());
      }
      
      //Transpose goes here
  };
};

template<int Rows, int Cols, class Precision>
struct RowMajor: public StaticSizedAllocator<Rows*Cols, Precision>
{
  int num_rows() const { return Rows; }
  int num_cols() const { return Cols; }

  // has to handle operator[row] and operator(row,col)

  using StaticSizedAllocator<Rows*Cols, Precision>::my_data;

  Vector<Cols, Precision, SVBase<Cols,1,Precision> > operator[](int row){
    Internal::check_index(Rows, row);
    return Vector<Cols, Precision, SVBase<Cols,1,Precision> >(my_data+Cols*row);
  }


  Precision& operator()(int row, int col){
    Internal::check_index(Rows, row);
    Internal::check_index(Cols, col);
    return my_data[row*Cols+col];
  }
  const Precision& operator()(int row, int col) const {
    Internal::check_index(Rows, row);
    Internal::check_index(Cols, col);
    return my_data[row*Cols+col];
  }


  template<int Rstart, int Cstart, int Rlength, int Clength>
  Matrix<Rlength, Clength, Precision, SSSlice<Cols>::template RowMajor> slice(){
    Internal::CheckSlice<Rows, Rstart, Rlength>::check();
    Internal::CheckSlice<Cols, Cstart, Clength>::check();

    return Matrix<Rlength, Clength, Precision, SSSlice<Cols>::template RowMajor>(my_data + Cols*Rstart + Cstart, Slicing());
  }
  
  //Transpose goes here
};
