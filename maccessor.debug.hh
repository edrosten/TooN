
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
#ifndef __MACCESSOR_HH
#define __MACCESSOR_HH


///////////////////////////////
//                           //
//  Matrix Accessor classes  //
//                           //
///////////////////////////////

#define TOON_CHECK_ROW TOON_ASSERT(r < Rows && r >= 0, TooNError::BadRowIndex)
#define TOON_CHECK_COL TOON_ASSERT(c < Cols && c >= 0, TooNError::BadColIndex)

//TODO: error checking in slices


template <int Rows,int Cols, class AllocZone>
class FixedMAccessor<Rows,Cols,RowMajor,AllocZone> : public AllocZone {
 public:
  inline FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) TOON_THROW {
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(my_values[r*Cols]);
  }
  inline const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) const TOON_THROW {
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(my_values[r*Cols]);
  }

  inline double& operator()(int r, int c) TOON_THROW {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*Cols+c];
  }

  inline const double& operator()(int r, int c) const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*Cols+c];
  }

  static inline int num_rows() throw() {return Rows;}
  static inline int num_cols() throw() {return Cols;}
  static inline int num_skip() throw() {return Cols;}
  typedef RowMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >& T() {
    return reinterpret_cast<FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >&>(*my_values);
  }
  inline const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >& T() const {
    return reinterpret_cast<const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >&>(*my_values);
  }

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >& slice(){
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >&>
      (my_values[Rstart*Cols+Cstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >& slice() const {
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >&>
      (my_values[Rstart*Cols+Cstart]);
  }
};



template <int Rows,int Cols, class AllocZone>
class FixedMAccessor<Rows,Cols,ColMajor,AllocZone> : public AllocZone {
 public:
  FixedVector<Cols, SkipAccessor<Cols, Rows> >& operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols, SkipAccessor<Cols, Rows> >&>(my_values[r]);
  }

  const FixedVector<Cols, SkipAccessor<Cols, Rows> >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols, SkipAccessor<Cols, Rows> >&>(my_values[r]);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[c*Rows+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[c*Rows+r];
  }

  static inline int num_rows(){return Rows;}
  static inline int num_cols(){return Cols;}
  static inline int num_skip(){return Rows;}
  typedef ColMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >& T() {
    return reinterpret_cast<FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >&>(*my_values);
  }
  inline const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >& T() const {
    return reinterpret_cast<const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >&>(*my_values);
  }

  // slice()
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >& slice(){
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >&>
      (my_values[Cstart*Rows+Rstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >& slice()const{
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >&>
      (my_values[Cstart*Rows+Rstart]);
  }
};





template<int Rows, int Cols, int Skip>
class SkipMAccessor<Rows,Cols,Skip,RowMajor> {
public:
  inline FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(my_values[r*Skip]);
  }
  inline const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(my_values[r*Skip]);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*Skip+c];
  }

  inline const double& operator()(int r, int c) const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*Skip+c];
  }

  static inline int num_rows(){return Rows;}
  static inline int num_cols(){return Cols;}
  static inline int num_skip(){return Skip;}
  typedef RowMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,ColMajor> >& T() {
    return reinterpret_cast<FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,ColMajor> >&>(*this); 
  }

  inline const FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,ColMajor> >& T() const {
    return reinterpret_cast<const FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,ColMajor> >&>(*this);
  }

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >& slice(){
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >&>
      (my_values[Rstart*Skip+Cstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >& slice() const {
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >&>
      (my_values[Rstart*Skip+Cstart]);
  }

 protected:
  double my_values[Rows*Skip];
};

template<int Rows, int Cols, int Skip> 
class SkipMAccessor<Rows, Cols, Skip, ColMajor> {
public:

  inline FixedVector<Cols, SkipAccessor<Cols, Skip> >& operator[](int r)  TOON_THROW
  {
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols, SkipAccessor<Cols, Skip> >&>(my_values[r]);
  }

  inline const FixedVector<Cols, SkipAccessor<Cols, Skip> >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols, SkipAccessor<Cols, Skip> >&>(my_values[r]);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[c*Skip+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[c*Skip+r];
  }

  static inline int num_rows(){return Rows;}
  static inline int num_cols(){return Cols;}
  static inline int num_skip(){return Skip;}
  typedef ColMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,RowMajor> >& T() {
    return reinterpret_cast<FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,RowMajor> >&>(*this); 
  }

  inline const FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,RowMajor> >& T() const {
    return reinterpret_cast<const FixedMatrix<Cols,Rows,SkipMAccessor<Cols,Rows,Skip,RowMajor> >&>(*this);
  }

  // slice()
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >& slice(){
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >&>
      (my_values[Cstart*Skip+Rstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >& slice()const{
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >&>
      (my_values[Cstart*Skip+Rstart]);
  }

 protected:
  double my_values[Cols*Skip];
};


/////////////////////////
//                     //
//  Dynamic Accessors  //
//                     //
/////////////////////////
#undef TOON_CHECK_ROW
#undef TOON_CHECK_COL
#define TOON_CHECK_ROW TOON_ASSERT(r < my_rows && r >= 0, TooNError::BadRowIndex)
#define TOON_CHECK_COL TOON_ASSERT(c < my_cols && c >= 0, TooNError::BadColIndex)

template <>
class DynamicMAccessor<RowMajor> {
 public:
  const RefVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return RefVector(my_num_cols,my_values+r*my_num_cols);
  }

  RefVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return RefVector(my_num_cols,my_values+r*my_num_cols);
  }

  inline double& operator()(int r, int c)TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return my_values[r*my_num_cols+c];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return my_values[r*my_num_cols+c];
  }

  int num_rows()const{return my_num_rows;}
  int num_cols()const{return my_num_cols;}
  int num_skip()const{return my_num_cols;}
  typedef RowMajor layout;

  inline DynamicMatrix<DynamicMAccessor<ColMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<DynamicMAccessor<ColMajor> >&>(*this);
  }
  inline const DynamicMatrix<DynamicMAccessor<ColMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<DynamicMAccessor<ColMajor> >&>(*this);
  }

 protected:
  int my_num_rows;
  int my_num_cols;
  double* my_values;
};

template <>
class DynamicMAccessor<ColMajor> {
 public:  
  const RefSkipVector operator[](int r) const TOON_THROW{ 
  	TOON_CHECK_ROW;
    return RefSkipVector(my_num_cols,my_num_rows,my_values+r);
  }

  RefSkipVector operator[](int r) TOON_THROW
  {
  	TOON_CHECK_ROW;
    return RefSkipVector(my_num_cols,my_num_rows,my_values+r);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[c*my_num_rows+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return my_values[c*my_num_rows+r];
  }

  int num_rows()const{return my_num_rows;}
  int num_cols()const{return my_num_cols;}
  int num_skip()const{return my_num_rows;}
  typedef ColMajor layout;

  inline DynamicMatrix<DynamicMAccessor<RowMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<DynamicMAccessor<RowMajor> >&>(*this);
  }
  inline const DynamicMatrix<DynamicMAccessor<RowMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<DynamicMAccessor<RowMajor> >&>(*this);
  }

 protected:
  int my_num_cols;  // representation of rows and cols transposed from above
  int my_num_rows;  // so that transpose operation is a simple cast
  double* my_values;
};

///////// RefSkipMAccessor ////////////

template <>
class RefSkipMAccessor<RowMajor> {
  friend class RefSkipVector;
 public:
  RefSkipMAccessor(){};
  
  const RefVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return RefVector(my_num_cols,my_values+r*my_skip);
  }

  RefVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return RefVector(my_num_cols,my_values+r*my_skip);
  }

  inline double& operator()(int r, int c) TOON_THROW{
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*my_skip+c];	
  }

  inline const double& operator()(int r, int c) const TOON_THROW {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return my_values[r*my_skip+c];	
  }

  int num_rows()const{return my_num_rows;}
  int num_cols()const{return my_num_cols;}
  int num_skip()const{return my_skip;}
  typedef RowMajor layout;

  inline DynamicMatrix<RefSkipMAccessor<ColMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<RefSkipMAccessor<ColMajor> >&>(*this);
  }

  inline const DynamicMatrix<RefSkipMAccessor<ColMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<RefSkipMAccessor<ColMajor> >&>(*this);
  }


 protected:
  int my_num_rows;
  int my_num_cols;
  int my_skip;
  double* my_values;
};

template <>
class RefSkipMAccessor<ColMajor> {
  friend class RefSkipVector;
 public:
  RefSkipMAccessor(){};
  
  const RefSkipVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return RefSkipVector(my_num_cols,my_skip,my_values+r);
  }

  RefSkipVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return RefSkipVector(my_num_cols,my_skip,my_values+r);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {	
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return my_values[c*my_skip+r];
  }

  inline const double& operator()(int r, int c) const TOON_THROW
  {	
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return my_values[c*my_skip+r];

  int num_rows()const{return my_num_rows;}
  int num_cols()const{return my_num_cols;}
  int num_skip()const{return my_skip;}
  typedef ColMajor layout;

  inline DynamicMatrix<RefSkipMAccessor<RowMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<RefSkipMAccessor<RowMajor> >&>(*this);
  }

  inline const DynamicMatrix<RefSkipMAccessor<RowMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<RefSkipMAccessor<RowMajor> >&>(*this);
  }


 protected:
  int my_num_cols;  // representation of rows and cols transposed from above
  int my_num_rows;  // so that transpose operation is a simple cast
  int my_skip;
  double* my_values;
};




#undef TOON_CHECK_ROW
#undef TOON_CHECK_COL


#endif
