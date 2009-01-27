// -*- c++ -*-

//Forward declarations

//TODO?
//maybe shift to having Layout have a layout template, so Matrix is Matrix<int,int,class,class>
template<int Rows, int Cols, class Precision, template<int, int, class> class Layout> class Matrix;

struct Slicing{};

//Static stride Static slice. Use closure to get Stride parameter
//in to RowMajor and ColMajor
template<int Stride> struct SSSlice
{
  template<int Rows, int Cols, class Precision> struct ColMajor;

  template<int Rows, int Cols, class Precision>
  struct RowMajor
  {
    protected: 
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
	  Matrix<Cols, Rows, Precision, SSSlice<Stride>::template ColMajor> T() {
		return Matrix<Cols, Rows, Precision, SSSlice<Stride>::template ColMajor>(my_data, Slicing());
	  }
  };

  template<int Rows, int Cols, class Precision>
  struct ColMajor
  {
    protected: 
      Precision* my_data;

    public:
      int num_rows() const { return Rows; }
      int num_cols() const { return Cols; }
      //Construction

      ColMajor(Precision* d)
      :my_data(d)
      {}

      //Indexing	

      Vector<Cols, Precision, SVBase<Cols,Stride,Precision> > operator[](int row){
        Internal::check_index(Rows, row);
        return Vector<Cols, Precision, SVBase<Cols,Stride,Precision> >(my_data+row);
      }

      Precision& operator()(int row, int col){
        Internal::check_index(Rows, row);
        Internal::check_index(Cols, col);
        return my_data[row+Stride*col];
      }

      const Precision& operator()(int row, int col) const {
        Internal::check_index(Rows, row);
        Internal::check_index(Cols, col);
        return my_data[row+Stride*col];
      }

      //Slices go here
      template<int Rstart, int Cstart, int Rlength, int Clength>
      Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template ColMajor> slice(){
        Internal::CheckSlice<Rows, Rstart, Rlength>::check();
        Internal::CheckSlice<Cols, Cstart, Clength>::check();

        return Matrix<Rlength, Clength, Precision, SSSlice<Stride>::template ColMajor>(my_data + Stride*Cstart + Rstart, Slicing());
      }
      
      //Transpose goes here
	  Matrix<Cols, Rows, Precision, SSSlice<Stride>::template RowMajor> T() {
		return Matrix<Cols, Rows, Precision, SSSlice<Stride>::template RowMajor>(my_data, Slicing());
	  }
	}; 
};

////////////////////////////////////////////////////////////////////////////////
//
// Storage base classes 
//

template<int Rows, int Cols, class Precision>
struct RowMajor: public StaticSizedAllocator<Rows*Cols, Precision>
{

  // has to handle operator[row] and operator(row,col)
  protected:
	using StaticSizedAllocator<Rows*Cols, Precision>::my_data;

  public:
	int num_rows() const { return Rows; }
	int num_cols() const { return Cols; }

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
  Matrix<Cols, Rows, Precision, SSSlice<Cols>::template ColMajor> T() {
  	return Matrix<Cols, Rows, Precision, SSSlice<Cols>::template ColMajor>(my_data, Slicing());
  }
};

template<int Rows, int Cols, class Precision>
struct ColMajor: public StaticSizedAllocator<Rows*Cols, Precision>
{

  // has to handle operator[row] and operator(row,col)
  protected:
	using StaticSizedAllocator<Rows*Cols, Precision>::my_data;

  public:
	int num_rows() const { return Rows; }
	int num_cols() const { return Cols; }

	Vector<Cols, Precision, SVBase<Cols,Rows,Precision> > operator[](int row){
	  Internal::check_index(Rows, row);
	  return Vector<Cols, Precision, SVBase<Cols,Rows,Precision> >(my_data+row);
	}


	Precision& operator()(int row, int col){
	  Internal::check_index(Rows, row);
	  Internal::check_index(Cols, col);
	  return my_data[row+col*Rows];
	}
	const Precision& operator()(int row, int col) const {
	  Internal::check_index(Rows, row);
	  Internal::check_index(Cols, col);
	  return my_data[row+col*Rows];
	}


	template<int Rstart, int Cstart, int Rlength, int Clength>
	Matrix<Rlength, Clength, Precision, SSSlice<Rows>::template ColMajor> slice(){
	  Internal::CheckSlice<Rows, Rstart, Rlength>::check();
	  Internal::CheckSlice<Cols, Cstart, Clength>::check();

	  return Matrix<Rlength, Clength, Precision, SSSlice<Rows>::template ColMajor>(my_data + Rstart + Cstart*Rows, Slicing());
	}
	
  //Transpose goes here
  Matrix<Cols, Rows, Precision, SSSlice<Rows>::template RowMajor> T() {
	return Matrix<Cols, Rows, Precision, SSSlice<Rows>::template RowMajor>(my_data, Slicing());
  }
};
