
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
#ifndef __MCLASSES_HH
#define __MCLASSES_HH



template <int Rows, int Cols, class Layout>
class Matrix : public FixedMatrix<Rows, Cols, FixedMAccessor<Rows, Cols, Layout, typename SizeTraits<Rows*Cols>::get_zone> > {
  public:
  inline Matrix(){}

  // construction from 1 and 2 dimensional c-style arrays
  inline Matrix(const double d[Rows*Cols]){*this = (reinterpret_cast<const FixedMatrix<Rows,Cols,FixedMAccessor<Rows,Cols,RowMajor,Stack<Rows*Cols> > >&>(*d));}
  inline Matrix(const double d[Rows][Cols]){*this = (reinterpret_cast<const FixedMatrix<Rows,Cols,FixedMAccessor<Rows,Cols,RowMajor,Stack<Rows*Cols> > >&>(*d));}

  // construction from a fixed Matrix
  template <class Accessor>
  inline Matrix(const FixedMatrix<Rows,Cols,Accessor>& m){
    FixedMatrix<Rows,Cols, FixedMAccessor<Rows, Cols, Layout, typename SizeTraits<Rows*Cols>::get_zone> >::operator=(m);
  }

  // construction from a dynamic Matrix
  template <class Accessor>
  inline Matrix(const DynamicMatrix<Accessor>& m){
    assert(m.num_rows()==Rows && m.num_cols() == Cols);
    FixedMatrix<Rows,Cols, FixedMAccessor<Rows, Cols, Layout, typename SizeTraits<Rows*Cols>::get_zone> >::operator=(m);
  }

  // constructor from 1-ary operator
  template <class Arg, class Op>
  Matrix(const Arg& arg, const Operator<Op>&){Op::eval(*this,arg);}

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  Matrix(const LHS& lhs, const RHS& rhs, const Operator<Op>&){Op::eval(*this,lhs,rhs);}



};

template <class Layout>
class Matrix<General, General, Layout> : public DynamicMatrix<DynamicMAccessor<Layout> > {
  friend struct MSizer;
 public:
  inline Matrix(){
    this->my_num_rows = 0;
    this->my_num_cols = 0;
    this->my_values = NULL;
  }

  inline Matrix(int num_rows, int num_cols){
    this->my_num_rows = num_rows;
    this->my_num_cols = num_cols;
    this->my_values = new double[num_rows*num_cols];
  }

  inline ~Matrix(){
    delete[] this->my_values;}

  // construction from a double*
  inline Matrix(int num_rows, int num_cols, const double* data){
    this->my_num_rows=num_rows;
    this->my_num_cols = num_cols;
    this->my_values = new double[num_rows*num_cols];
    RefMatrix<RowMajor> temp (num_rows,num_cols,data);
    DynamicMatrix<DynamicMAccessor<Layout> >::operator=(temp);
  }

  // construction from a general Matrix
  template <class Accessor>
  inline Matrix(const MatrixBase<Accessor>& from){
    this->my_num_rows = from.num_rows();
    this->my_num_cols = from.num_cols();
    this->my_values = new double [this->my_num_rows*this->my_num_cols];
    DynamicMatrix<DynamicMAccessor<Layout> >::operator=(from);
  }

  inline void resize(int rows, int cols)
  {
  	if(rows != this->my_num_rows || cols != this->my_num_cols)
	{
		delete[] this->my_values;
		this->my_num_rows = rows;
		this->my_num_cols = cols;
		this->my_values = new double[this->my_num_rows*this->my_num_cols];
	}
  }
    
    template <class Accessor> inline void assign(const MatrixBase<Accessor>& m) {
	resize(m.num_rows(),m.num_cols());
	DynamicMatrix<DynamicMAccessor<Layout> >::operator=(m);	
    }
    
  // copy construction
  inline Matrix(const Matrix<General,General,Layout>& from){
    this->my_num_rows = from.my_num_rows;
    this->my_num_cols = from.my_num_cols;
    this->my_values = new double[this->my_num_rows*this->my_num_cols];
    DynamicMatrix<DynamicMAccessor<Layout> >::operator=(from);
  }

  // constructor from 1-ary operator
  template <class Arg, class Op>
  Matrix(const Arg& arg, const Operator<Op>&){
    Op::eval(*this,arg);
  }

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  Matrix(const LHS& lhs, const RHS& rhs, const Operator<Op>&){
    Op::eval(*this,lhs,rhs);
  }
};



struct MSizer {
  template<class Layout>
  static inline void set_size(Matrix<General,General,Layout>& ret, int rows, int cols){
    ret.set(rows, cols, new double[rows*cols]);
  }
};


#endif
