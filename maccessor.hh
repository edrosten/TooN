
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

/////////////////////////
//                     //
//  Dynamic Accessors  //
//                     //
/////////////////////////
#define TOON_CHECK_ROW TOON_ASSERT(r < my_num_rows && r >= 0, TooNError::BadRowIndex)
#define TOON_CHECK_COL TOON_ASSERT(c < my_num_cols && c >= 0, TooNError::BadColIndex)


///////// RefSkipMAccessor ////////////

template <>
class RefSkipMAccessor<RowMajor> {
 public:
  typedef DynamicMatrix<RefSkipMAccessor<RowMajor> > RefSkipMatrixRM;
  RefSkipMAccessor(){};
  
  const RefVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefVector(my_num_cols,this->my_values+r*my_skip);
  }

  RefVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefVector(my_num_cols,this->my_values+r*my_skip);
  }

  inline double& operator()(int r, int c) TOON_THROW{
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*my_skip+c];	
  }

  inline const double& operator()(int r, int c) const TOON_THROW{
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*my_skip+c];	
  }

  int num_rows()const throw() {return my_num_rows;}
  int num_cols()const throw() {return my_num_cols;}
  int num_skip()const throw() {return my_skip;}
  typedef RowMajor layout;

  inline DynamicMatrix<RefSkipMAccessor<ColMajor> >& T();

  inline const DynamicMatrix<RefSkipMAccessor<ColMajor> >& T() const;

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline RefSkipMatrixRM slice(){
    RefSkipMatrixRM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart*this->my_skip + Cstart);
    return ret;
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const RefSkipMatrixRM slice() const {
    RefSkipMatrixRM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart*this->my_skip + Cstart);
    return ret;
  }

  inline RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    RefSkipMatrixRM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart*this->my_skip + Cstart);
    return ret;
  }

  inline const RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    RefSkipMatrixRM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart*this->my_skip + Cstart);
    return ret;
  }
  inline void set(int nr, int nc, int sk, double* v) { my_num_rows = nr; my_num_cols = nc; my_skip = sk; my_values = v; }
 protected:
  int my_num_rows;
  int my_num_cols;
  int my_skip;
  double* my_values;
};

template <>
class RefSkipMAccessor<ColMajor> {
 public:
  typedef DynamicMatrix<RefSkipMAccessor<ColMajor> > RefSkipMatrixCM;

  RefSkipMAccessor(){};
  
  const RefSkipVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefSkipVector(my_num_cols,my_skip,this->my_values+r);
  }

  RefSkipVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefSkipVector(my_num_cols,my_skip,this->my_values+r);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {	
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[c*my_skip+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[c*my_skip+r];
  }


  int num_rows()const throw() {return my_num_rows;}
  int num_cols()const throw() {return my_num_cols;}
  int num_skip()const throw() {return my_skip;}
  typedef ColMajor layout;

  inline DynamicMatrix<RefSkipMAccessor<RowMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<RefSkipMAccessor<RowMajor> >&>(*this);
  }

  inline const DynamicMatrix<RefSkipMAccessor<RowMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<RefSkipMAccessor<RowMajor> >&>(*this);
  }

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline RefSkipMatrixCM slice(){
    RefSkipMatrixCM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart + Cstart*this->my_skip);
    return ret;
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const RefSkipMatrixCM slice() const {
    RefSkipMatrixCM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart + Cstart*this->my_skip);
    return ret;
  }

  inline RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    RefSkipMatrixCM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart + Cstart*this->my_skip);
    return ret;
  }

  inline const RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    RefSkipMatrixCM ret;
    ret.set(Rsize, Csize, this->my_skip, this->my_values + Rstart + Cstart*this->my_skip);
    return ret;
  }
  inline void set(int nr, int nc, int sk, double* v) { my_num_rows = nr; my_num_cols = nc; my_skip = sk; my_values = v; }

 protected:
  int my_num_cols;  // representation of rows and cols transposed from above
  int my_num_rows;  // so that transpose operation is a simple cast
  int my_skip;
  double* my_values;
};

inline DynamicMatrix<RefSkipMAccessor<ColMajor> >& RefSkipMAccessor<RowMajor>::T() {
    return reinterpret_cast<DynamicMatrix<RefSkipMAccessor<ColMajor> >&>(*this);
}

inline const DynamicMatrix<RefSkipMAccessor<ColMajor> >& RefSkipMAccessor<RowMajor>::T() const {
    return reinterpret_cast<const DynamicMatrix<RefSkipMAccessor<ColMajor> >&>(*this);
}

typedef RefSkipMAccessor<RowMajor>::RefSkipMatrixRM RefSkipMatrixRM;
typedef RefSkipMAccessor<ColMajor>::RefSkipMatrixCM RefSkipMatrixCM;
inline RefSkipMatrixRM makeRefSkipMatrixRM(int nr, int nc, int sk, double* v) { RefSkipMatrixRM ret; ret.set(nr,nc,sk,v); return ret; }
inline RefSkipMatrixCM makeRefSkipMatrixCM(int nr, int nc, int sk, double* v) { RefSkipMatrixCM ret; ret.set(nr,nc,sk,v); return ret; }

template <>
class DynamicMAccessor<RowMajor> {
 public:
  const RefVector operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefVector(my_num_cols,this->my_values+r*my_num_cols);
  }

  RefVector operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return makeRefVector(my_num_cols,this->my_values+r*my_num_cols);
  }

  inline double& operator()(int r, int c)TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return this->my_values[r*my_num_cols+c];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return this->my_values[r*my_num_cols+c];
  }

  int num_rows()const throw() {return my_num_rows;}
  int num_cols()const throw() {return my_num_cols;}
  int num_skip()const throw() {return my_num_cols;}
  typedef RowMajor layout;

  inline DynamicMatrix<DynamicMAccessor<ColMajor> >& T();
  inline const DynamicMatrix<DynamicMAccessor<ColMajor> >& T() const;

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline RefSkipMatrixRM slice(){
    return makeRefSkipMatrixRM(Rsize, Csize, this->my_num_cols, this->my_values + Rstart*this->my_num_cols + Cstart);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const RefSkipMatrixRM slice() const {
    return makeRefSkipMatrixRM(Rsize, Csize, this->my_num_cols, this->my_values + Rstart*this->my_num_cols + Cstart);
  }

  inline RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixRM(Rsize, Csize, this->my_num_cols, this->my_values + Rstart*this->my_num_cols + Cstart);
  }

  inline const RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixRM(Rsize, Csize, this->my_num_cols, this->my_values + Rstart*this->my_num_cols + Cstart);
  }
  inline void set(int nr, int nc, double* v) { my_num_rows = nr; my_num_cols = nc; my_values = v; }
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
    return makeRefSkipVector(my_num_cols,my_num_rows,this->my_values+r);
  }

  RefSkipVector operator[](int r) TOON_THROW
  {
  	TOON_CHECK_ROW;
    return makeRefSkipVector(my_num_cols,my_num_rows,this->my_values+r);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[c*my_num_rows+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
  	TOON_CHECK_COL;
  	return this->my_values[c*my_num_rows+r];
  }

  int num_rows()const throw() {return my_num_rows;}
  int num_cols()const throw() {return my_num_cols;}
  int num_skip()const throw() {return my_num_rows;}
  typedef ColMajor layout;

  inline DynamicMatrix<DynamicMAccessor<RowMajor> >& T() {
    return reinterpret_cast<DynamicMatrix<DynamicMAccessor<RowMajor> >&>(*this);
  }
  inline const DynamicMatrix<DynamicMAccessor<RowMajor> >& T() const {
    return reinterpret_cast<const DynamicMatrix<DynamicMAccessor<RowMajor> >&>(*this);
  }

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline RefSkipMatrixCM slice(){
    return makeRefSkipMatrixCM(Rsize, Csize, this->my_num_rows, this->my_values + Rstart + Cstart*this->my_num_rows);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const RefSkipMatrixCM slice() const {
    return makeRefSkipMatrixCM(Rsize, Csize, this->my_num_rows, this->my_values + Rstart + Cstart*this->my_num_rows);
  }

  inline RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixCM(Rsize, Csize, this->my_num_rows, this->my_values + Rstart + Cstart*this->my_num_rows);
  }

  inline const RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixCM(Rsize, Csize, this->my_num_rows, this->my_values + Rstart + Cstart*this->my_num_rows);
  }

  inline void set(int nr, int nc, double* v) { my_num_rows = nr; my_num_cols = nc; my_values = v; }
 protected:
  int my_num_cols;  // representation of rows and cols transposed from above
  int my_num_rows;  // so that transpose operation is a simple cast
  double* my_values;
};

inline DynamicMatrix<DynamicMAccessor<ColMajor> >& DynamicMAccessor<RowMajor>::T() {
    return reinterpret_cast<DynamicMatrix<DynamicMAccessor<ColMajor> >&>(*this);
}
inline const DynamicMatrix<DynamicMAccessor<ColMajor> >& DynamicMAccessor<RowMajor>::T() const {
    return reinterpret_cast<const DynamicMatrix<DynamicMAccessor<ColMajor> >&>(*this);
}


typedef DynamicMatrix<DynamicMAccessor<RowMajor> > RefMatrixRM;
typedef DynamicMatrix<DynamicMAccessor<ColMajor> > RefMatrixCM;
inline RefMatrixRM makeRefMatrixRM(int nr, int nc, double* v) { RefMatrixRM ret; ret.set(nr,nc,v); return ret; }
inline RefMatrixCM makeRefMatrixCM(int nr, int nc, double* v) { RefMatrixCM ret; ret.set(nr,nc,v); return ret; }

#undef TOON_CHECK_ROW
#undef TOON_CHECK_COL
#define TOON_CHECK_ROW TOON_ASSERT(r < Rows && r >= 0, TooNError::BadRowIndex)
#define TOON_CHECK_COL TOON_ASSERT(c < Cols && c >= 0, TooNError::BadColIndex)




template<int Rows, int Cols, int Skip>
class SkipMAccessor<Rows,Cols,Skip,RowMajor> {
public:
  inline FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(this->my_values[r*Skip]);
  }
  inline const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(this->my_values[r*Skip]);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*Skip+c];
  }

  inline const double& operator()(int r, int c) const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*Skip+c];
  }

  static inline int num_rows() throw() {return Rows;}
  static inline int num_cols() throw() {return Cols;}
  static inline int num_skip() throw() {return Skip;}
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
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >& slice() {
      FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >::dummy();
     util::Assert<(Rstart+Rsize <= Rows)>();
     util::Assert<(Cstart+Csize <= Cols)>();
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >&>
      (this->my_values[Rstart*Skip+Cstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >& slice() const {
     FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >::dummy();
     util::Assert<(Rstart+Rsize <= Rows)>();
     util::Assert<(Cstart+Csize <= Cols)>();
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,RowMajor> >&>
      (this->my_values[Rstart*Skip+Cstart]);
  }

  inline RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixRM(Rsize, Csize, Skip, this->my_values + Rstart*Skip + Cstart);
  }

  inline const RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixRM(Rsize, Csize, Skip, this->my_values + Rstart*Skip + Cstart);
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
    return reinterpret_cast<FixedVector<Cols, SkipAccessor<Cols, Skip> >&>(this->my_values[r]);
  }

  inline const FixedVector<Cols, SkipAccessor<Cols, Skip> >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols, SkipAccessor<Cols, Skip> >&>(this->my_values[r]);
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
  	return this->my_values[c*Skip+r];
  }

  static inline int num_rows() throw() {return Rows;}
  static inline int num_cols() throw() {return Cols;}
  static inline int num_skip() throw() {return Skip;}
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
      (this->my_values[Cstart*Skip+Rstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >& slice()const{
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Skip,ColMajor> >&>
      (this->my_values[Cstart*Skip+Rstart]);
  }

  inline RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixCM(Rsize, Csize, Skip, this->my_values + Rstart + Cstart*Skip);
  }

  inline const RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixCM(Rsize, Csize, Skip, this->my_values + Rstart + Cstart*Skip);
  }

 protected:
  double my_values[Cols*Skip];
};


template <int Rows,int Cols, class AllocZone>
class FixedMAccessor<Rows,Cols,RowMajor,AllocZone> : public AllocZone {
 public:
  inline FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) TOON_THROW {
  	TOON_CHECK_ROW;
    return reinterpret_cast<FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(this->my_values[r*Cols]);
  }
  inline const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >& operator[](int r) const TOON_THROW {
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols,FixedVAccessor<Cols,Stack<Cols> > >&>(this->my_values[r*Cols]);
  }

  inline double& operator()(int r, int c) TOON_THROW {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*Cols+c];
  }

  inline const double& operator()(int r, int c) const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[r*Cols+c];
  }

  static inline int num_rows() throw() {return Rows;}
  static inline int num_cols() throw() {return Cols;}
  static inline int num_skip() throw() {return Cols;}
  typedef RowMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >& T() {
      FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >::dummy();
    return reinterpret_cast<FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >&>(*this->my_values);
  }
  inline const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >& T() const {
      FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >::dummy();
      return reinterpret_cast<const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,ColMajor,Stack<Rows*Cols> > >&>(*this->my_values);
  }

  // slice
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >& slice(){
      typedef FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> > ST;
      util::Assert<(Rstart+Rsize <= Rows)>();
      util::Assert<(Cstart+Csize <= Cols)>();
      ST::dummy();
      return reinterpret_cast<ST&>(this->my_values[Rstart*Cols+Cstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >& slice() const {
      util::Assert<(Rstart+Rsize <= Rows)>();
      util::Assert<(Cstart+Csize <= Cols)>();
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Cols,RowMajor> >&>
      (this->my_values[Rstart*Cols+Cstart]);
  }

  inline RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixRM(Rsize, Csize, Cols, this->my_values + Rstart*Cols + Cstart);
  }

  inline const RefSkipMatrixRM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixRM(Rsize, Csize, Cols, this->my_values + Rstart*Cols + Cstart);
  }

};



template <int Rows,int Cols, class AllocZone>
class FixedMAccessor<Rows,Cols,ColMajor,AllocZone> : public AllocZone {
 public:
  FixedVector<Cols, SkipAccessor<Cols, Rows> >& operator[](int r) TOON_THROW{
  	TOON_CHECK_ROW;
	FixedVector<Cols, SkipAccessor<Cols, Rows> >::dummy();
	return reinterpret_cast<FixedVector<Cols, SkipAccessor<Cols, Rows> >&>(this->my_values[r]);
  }

  const FixedVector<Cols, SkipAccessor<Cols, Rows> >& operator[](int r) const TOON_THROW{
  	TOON_CHECK_ROW;
    return reinterpret_cast<const FixedVector<Cols, SkipAccessor<Cols, Rows> >&>(this->my_values[r]);
  }

  inline double& operator()(int r, int c) TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[c*Rows+r];
  }

  inline const double& operator()(int r, int c)const TOON_THROW
  {
  	TOON_CHECK_ROW;
	TOON_CHECK_COL;
  	return this->my_values[c*Rows+r];
  }

  static inline int num_rows() throw() {return Rows;}
  static inline int num_cols() throw() {return Cols;}
  static inline int num_skip() throw() {return Rows;}
  typedef ColMajor layout;

  // Transpose operations
  inline FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >& T() {
      FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >::dummy();
      return reinterpret_cast<FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >&>(*this->my_values);
  }
  inline const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >& T() const {
      FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >::dummy();
      return reinterpret_cast<const FixedMatrix<Cols,Rows,FixedMAccessor<Cols,Rows,RowMajor,Stack<Rows*Cols> > >&>(*this->my_values);
  }

  // slice()
  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >& slice(){
      FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >::dummy();
    return reinterpret_cast<FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >&>
      (this->my_values[Cstart*Rows+Rstart]);
  }

  template<int Rstart, int Cstart, int Rsize, int Csize>
  inline const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >& slice()const{
      FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >::dummy();
    return reinterpret_cast<const FixedMatrix<Rsize,Csize,SkipMAccessor<Rsize,Csize,Rows,ColMajor> >&>
      (this->my_values[Cstart*Rows+Rstart]);
  }

  inline RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) {
    return makeRefSkipMatrixCM(Rsize, Csize, Rows, this->my_values + Rstart + Cstart*Rows);
  }

  inline const RefSkipMatrixCM slice(int Rstart, int Cstart, int Rsize, int Csize) const {
    return makeRefSkipMatrixCM(Rsize, Csize, Rows, this->my_values + Rstart + Cstart*Rows);
  }

};



#undef TOON_CHECK_ROW
#undef TOON_CHECK_COL


#endif
