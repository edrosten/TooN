
/*                       
	 Copyright (C) 2005 Tom Drummond

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __MBASE_HH
#define __MBASE_HH

template <class Accessor>
class MatrixBase : public Accessor {
 public:
  const double* get_data_ptr()const{return this->my_values;}
  double* get_data_ptr(){return this->my_values;}

  DynamicMatrix<Accessor>& down_cast() {return reinterpret_cast<DynamicMatrix<Accessor>&> (*this);}
  const DynamicMatrix<Accessor>& down_cast() const {return reinterpret_cast<const DynamicMatrix<Accessor>&> (*this);}
};

// operator ostream& <<
template <class Accessor>
std::ostream& operator << (std::ostream& os, const MatrixBase<Accessor>& m){
  for(int r=0; r<m.num_rows(); r++){
    os << m[r] << std::endl;
  }
  return os;
}

// operator istream& >>
template <class Accessor>
std::istream& operator >> (std::istream& is, MatrixBase<Accessor>& m){
  for(int r=0; r<m.num_rows(); r++){
    is >> m[r];
  }
  return is;
}

// Data copying

template <class Accessor1, class Accessor2>
struct MatrixCopy {
  inline static void eval(MatrixBase<Accessor1>& to, const MatrixBase<Accessor2>& from){
    for(int r=0; r<from.num_rows(); r++){
      for(int c=0; c<from.num_cols(); c++){
	to(r,c)=from(r,c);
      }
    }
  }
};

template<int Rows, int Cols, class Layout, class Zone1, class Zone2>
struct MatrixCopy<FixedMAccessor<Rows,Cols,Layout,Zone1>,FixedMAccessor<Rows,Cols,Layout,Zone2> > {
  inline static void eval(MatrixBase<FixedMAccessor<Rows,Cols,Layout,Zone1> >& to,
			  const MatrixBase<FixedMAccessor<Rows,Cols,Layout,Zone2> >& from) {
    memcpy(to.get_data_ptr(), from.get_data_ptr(), Rows*Cols*sizeof(double));
  }
};

template<class Layout>
struct MatrixCopy<DynamicMAccessor<Layout>,DynamicMAccessor<Layout> > {
  inline static void eval(MatrixBase<DynamicMAccessor<Layout> >& to,
			  const MatrixBase<DynamicMAccessor<Layout> >& from){
    memcpy(to.get_data_ptr(), from.get_data_ptr(), from.num_rows()*from.num_cols()*sizeof(double));
  }
};


////////////////////////////////////////////////////////
//                                                    //
// Fixed and Dynamic Matrix classes                   //
// all the arithmetic and assignment                  //
// operations are applied to these                    //
//                                                    //
////////////////////////////////////////////////////////

template<int Rows, int Cols, class Accessor>
class FixedMatrix : public MatrixBase<Accessor> {
public:
  // assignment from correct sized FixedMatrix
  template<class Accessor2>
  inline FixedMatrix& operator=(const FixedMatrix<Rows,Cols,Accessor2>& from){
    MatrixCopy<Accessor, Accessor2>::eval(*this,from);
    return *this;
  }

  // copy assignment
  inline FixedMatrix& operator=(const FixedMatrix& from){
    MatrixCopy<Accessor,Accessor>::eval(*this,from);
    return *this;
  }
  
  // assignment from any DynamicMatrix
  template<class Accessor2>
    inline FixedMatrix& operator=(const DynamicMatrix<Accessor2>& from){
    assert(from.num_rows() == Rows && from.num_cols() == Cols);
    MatrixCopy<Accessor,Accessor2>::eval(*this,from);
    return *this;
  }
    static void dummy() {}
};

template<class Accessor>
class DynamicMatrix : public MatrixBase<Accessor> {
public:
  // assignment from any MatrixBase
  template<class Accessor2>
  DynamicMatrix& operator=(const MatrixBase<Accessor2>& from){
    assert(from.num_rows() == this->num_rows() && from.num_cols() == this->num_cols());
    MatrixCopy<Accessor,Accessor2>::eval(*this,from);
    return *this;
  }


  // repeated for explicit copy assignment
  DynamicMatrix& operator=(const DynamicMatrix& from){
    assert(from.num_rows() == this->num_rows() && from.num_cols() == this->num_cols());
    MatrixCopy<Accessor,Accessor>::eval(*this,from);
    return *this;
  }

    operator DynamicMatrix& () { return *this; }
    operator const DynamicMatrix& () const { return *this; }

};


#endif
