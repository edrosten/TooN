
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
#ifndef __LINOPERATORS_HH
#define __LINOPERATORS_HH

inline RefMatrixRM DynamicVAccessor::as_row(){
  return makeRefMatrixRM(1,my_size,my_values);
}

inline RefMatrixCM DynamicVAccessor::as_col(){
  return makeRefMatrixCM(my_size,1,my_values);
}

inline RefSkipMatrixCM DynamicSkipAccessor::as_row(){
  return makeRefSkipMatrixCM(1,my_size,my_skip,my_values);
}

inline RefSkipMatrixRM DynamicSkipAccessor::as_col(){
  return makeRefSkipMatrixRM(my_size,1,my_skip,my_values);
}

//////////////////////
// unary operator - //
//    -Vector       //
//////////////////////

template<int Size, class Accessor>
struct FixedVectorNeg{
  inline static void eval(Vector<Size>& ret, const FixedVector<Size,Accessor>& arg){
    for(int i=0; i<Size; i++){
      ret[i]=-arg[i];
    }
  }
};

template<int Size, class Accessor> inline
Vector<Size> operator-(const FixedVector<Size,Accessor>& arg){
  return Vector<Size>(arg,Operator<FixedVectorNeg<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVectorNeg : public VSizer {
  inline static void eval(Vector<>& ret, const DynamicVector<Accessor>& arg){
    const int Size = arg.size();
    set_size(ret,Size);
    for(int i=0; i<Size; i++){
      ret[i]=-arg[i];
    }
  }
};


template <class Accessor>
Vector<> operator-(const DynamicVector<Accessor>& arg){
  return Vector<>(arg,Operator<DynamicVectorNeg<Accessor> >());
}

/////////////////////
//                 //
//   operator ==   //
//   operator !=   //
// Vector == Vector//
/////////////////////


template<int Size, class Accessor1, class Accessor2>
inline bool operator==(const FixedVector<Size, Accessor1>& lhs, const FixedVector<Size, Accessor2>& rhs){
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 0;
  return 1;
}

template<int Size, class Accessor1, class Accessor2>
inline bool operator==(const FixedVector<Size, Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(rhs.size() == Size);
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 0;
  return 1;
}

template<class Accessor1, class Accessor2>
inline bool operator==(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(rhs.size() == lhs.size());
  for(int i=0; i < rhs.size(); i++)
  	if(lhs[i] != rhs[i])
		return 0;
  return 1;
}

template<int Size, class Accessor1, class Accessor2>
inline bool operator==(const DynamicVector<Accessor1>& lhs, const FixedVector<Size, Accessor2>& rhs){
  assert(Size == lhs.size());
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 0;
  return 1;
}



template<int Size, class Accessor1, class Accessor2>
inline bool operator!=(const FixedVector<Size, Accessor1>& lhs, const FixedVector<Size, Accessor2>& rhs){
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 1;
  return 0;
}

template<int Size, class Accessor1, class Accessor2>
inline bool operator!=(const FixedVector<Size, Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(rhs.size() == Size);
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 1;
  return 0;
}

template<class Accessor1, class Accessor2>
inline bool operator!=(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(rhs.size() == lhs.size());
  for(int i=0; i < rhs.size(); i++)
  	if(lhs[i] != rhs[i])
		return 1;
  return 0;
}

template<int Size, class Accessor1, class Accessor2>
inline bool operator!=(const DynamicVector<Accessor1>& lhs, const FixedVector<Size, Accessor2>& rhs){
  assert(Size == lhs.size());
  for(int i=0; i < Size; i++)
  	if(lhs[i] != rhs[i])
		return 1;
  return 0;
}




/////////////////////
//                 //
//   operator +    //
// Vector + Vector //
/////////////////////

template <int Size, class LHS, class RHS>
struct FixedVAdd {
  inline static void eval(Vector<Size>& ret,
		   const LHS& lhs,
		   const RHS& rhs){
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i]+rhs[i];
    }
  }
};


template <int Size, class Accessor1, class Accessor2>
inline Vector<Size> operator+(const FixedVector<Size,Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs) {
  return Vector<Size>(lhs, rhs,
		      Operator<FixedVAdd<Size,
		      FixedVector<Size,Accessor1>,
		      FixedVector<Size,Accessor2> > >());
}


template <class LHS, class RHS>
struct DynamicVAdd : public VSizer{
  static void eval(Vector<>& ret,
		   const LHS& lhs,
		   const RHS& rhs){
    assert(lhs.size() == rhs.size());
    const int Size = lhs.size();
    set_size(ret,Size);
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i]+rhs[i];
    }
  }
};

// have to specify operator + three more times
// for any combination of Fixed and Dynamic except Fixed-Fixed

template <class Accessor1,class Accessor2>
inline Vector<> operator+(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  return Vector<>(lhs,rhs,
		  Operator<DynamicVAdd<
		  DynamicVector<Accessor1>,
		  DynamicVector<Accessor2> > >());
}

template <int Size, class Accessor1,class Accessor2>
inline Vector<Size> operator+(const FixedVector<Size,Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
    assert(Size == rhs.size());
    return Vector<Size>(lhs,rhs,
		  Operator<FixedVAdd<Size,
		  FixedVector<Size,Accessor1>,
		  DynamicVector<Accessor2> > >());
}

template <class Accessor1, int Size, class Accessor2>
inline Vector<Size> operator+(const DynamicVector<Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs){
  assert(Size == lhs.size());
  return Vector<Size>(lhs,rhs,
		  Operator<FixedVAdd<Size,
		  DynamicVector<Accessor1>,
		  FixedVector<Size,Accessor2> > >());
}

//////////////////////
//                  //
// operator +=      //
// Vector += Vector //
//                  //
//////////////////////

template <int Size, class LHAccessor, class RHAccessor>
inline FixedVector<Size,LHAccessor>& operator+= (FixedVector<Size,LHAccessor>& lhs,
						 const FixedVector<Size,RHAccessor>& rhs){
    util::AddV<0,Size-1>::eval(lhs,rhs);
    return lhs;
}

template <int Size, class LHAccessor, class RHAccessor>
inline FixedVector<Size,LHAccessor>& operator+= (FixedVector<Size,LHAccessor>& lhs,
						 const DynamicVector<RHAccessor>& rhs){
    assert(rhs.size() == Size);
    util::AddV<0,Size-1>::eval(lhs,rhs);
    return lhs;
}

template <int Size, class LHAccessor, class RHAccessor>
inline DynamicVector<LHAccessor>& operator+= (DynamicVector<LHAccessor>& lhs,
					      const FixedVector<Size,RHAccessor>& rhs){
  assert(lhs.size()==Size);
  util::AddV<0,Size-1>::eval(lhs,rhs);
  return lhs;
}

template <class LHAccessor, class RHAccessor>
inline DynamicVector<LHAccessor>& operator+= (DynamicVector<LHAccessor>& lhs,
					      const DynamicVector<RHAccessor>& rhs){
  assert(lhs.size()==rhs.size());
  for(int i=0; i<lhs.size(); i++){
    lhs[i]+=rhs[i];
  }
  return lhs;
}


////////////////
//            //
// operator - //
//            //
////////////////

template <int Size, class LHS, class RHS>
struct FixedVSub {
  static void eval(Vector<Size>& ret,
		   const LHS& lhs,
		   const RHS& rhs){
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i]-rhs[i];
    }
  }
};


template <int Size, class Accessor1, class Accessor2>
inline Vector<Size> operator-(const FixedVector<Size,Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs) {
  return Vector<Size>(lhs, rhs,
		      Operator<FixedVSub<Size,
		      FixedVector<Size,Accessor1>,
		      FixedVector<Size,Accessor2> > >());
}


template <class LHS, class RHS>
struct DynamicVSub : public VSizer{
  static void eval(Vector<>& ret,
		   const LHS& lhs,
		   const RHS& rhs){
    assert(lhs.size() == rhs.size());
    const int Size = lhs.size();
    set_size(ret,Size);
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i]-rhs[i];
    }
  }
};

// have to specify operator - three more times
// for any combination of Fixed and Dynamic except Fixed-Fixed

template <class Accessor1,class Accessor2>
inline Vector<> operator-(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  return Vector<>(lhs,rhs,
		  Operator<DynamicVSub<
		  DynamicVector<Accessor1>,
		  DynamicVector<Accessor2> > >());
}

template <int Size, class Accessor1,class Accessor2>
inline Vector<Size> operator-(const FixedVector<Size,Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(Size == rhs.size());
  return Vector<Size>(lhs,rhs,
		  Operator<FixedVSub<Size,
		  FixedVector<Size,Accessor1>,
		  DynamicVector<Accessor2> > >());
}

template <class Accessor1, int Size, class Accessor2>
inline Vector<Size> operator-(const DynamicVector<Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs){
  assert(Size == lhs.size());
  return Vector<Size>(lhs,rhs,
		  Operator<FixedVSub<Size,
		  DynamicVector<Accessor1>,
		  FixedVector<Size,Accessor2> > >());
}


/////////////////
//             //
// operator -= //
//             //
/////////////////

template <int Size, class LHAccessor, class RHAccessor>
inline FixedVector<Size,LHAccessor>& operator-= (FixedVector<Size,LHAccessor>& lhs,
						 const FixedVector<Size,RHAccessor>& rhs){
  for(int i=0; i<Size; i++){
    lhs[i]-=rhs[i];
  }
  return lhs;
}

template <int Size, class LHAccessor, class RHAccessor>
inline FixedVector<Size,LHAccessor>& operator-= (FixedVector<Size,LHAccessor>& lhs,
						 const DynamicVector<RHAccessor>& rhs){
  assert(rhs.size() == Size);
  for(int i=0; i<Size; i++){
    lhs[i]-=rhs[i];
  }
  return lhs;
}

template <int Size, class LHAccessor, class RHAccessor>
inline DynamicVector<LHAccessor>& operator-= (DynamicVector<LHAccessor>& lhs,
					      const FixedVector<Size,RHAccessor>& rhs){
  assert(lhs.size()==Size);
  for(int i=0; i<Size; i++){
    lhs[i]-=rhs[i];
  }
  return lhs;
}

template <class LHAccessor, class RHAccessor>
inline DynamicVector<LHAccessor>& operator-= (DynamicVector<LHAccessor>& lhs,
					      const DynamicVector<RHAccessor>& rhs){
  assert(lhs.size()==rhs.size());
  for(int i=0; i<lhs.size(); i++){
    lhs[i]-=rhs[i];
  }
  return lhs;
}


/////////////////////
// operator *      //
// Vector * double //
// (or v.v.)       //
/////////////////////

template <int Size, class Accessor>
struct FixedVdMult {
  static inline void eval(Vector<Size>& ret, const FixedVector<Size,Accessor>& lhs, double rhs){
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i] * rhs;
    }
  }
};

template <int Size, class Accessor>
inline Vector<Size> operator*(const FixedVector<Size,Accessor>& lhs, double rhs){
  return Vector<Size>(lhs,rhs,
		      Operator<FixedVdMult<Size,Accessor> >());
}
template <int Size, class Accessor>
inline Vector<Size> operator*(double lhs, const FixedVector<Size,Accessor>& rhs){
  return Vector<Size>(rhs,lhs,
		      Operator<FixedVdMult<Size,Accessor> >());
}

template <class Accessor>
struct DynamicVdMult : public VSizer{
  static inline void eval(Vector<>& ret, const DynamicVector<Accessor>& lhs, double rhs){
    set_size(ret,lhs.size());
    for(int i=0; i<lhs.size(); i++){
      ret[i] = lhs[i] * rhs;
    }
  }
};

template <class Accessor>
inline Vector<> operator*(const DynamicVector<Accessor>& lhs, double rhs){
  return Vector<>(lhs,rhs,
		  Operator<DynamicVdMult<Accessor> >());
}
template <class Accessor>
inline Vector<> operator*(double lhs, const DynamicVector<Accessor>& rhs){
  return Vector<>(rhs,lhs,
		  Operator<DynamicVdMult<Accessor> >());
}

//////////////////////
// operator *=      //
// Vector *= double //
//////////////////////

template <int Size, class Accessor>
inline FixedVector<Size,Accessor>& operator*=(FixedVector<Size,Accessor>& lhs, const double& rhs){
  for(int i=0; i<Size; i++){
    lhs[i]*=rhs;
  }
  return lhs;
}

template <class Accessor>
inline DynamicVector<Accessor>& operator*=(DynamicVector<Accessor>& lhs, const double& rhs){
  const int Size = lhs.size();
  for(int i=0; i<Size; i++){
    lhs[i]*=rhs;
  }
  return lhs;
}

/////////////////////
// operator /      //
// Vector / double //
/////////////////////

template <int Size, class Accessor>
inline Vector<Size> operator/(const FixedVector<Size,Accessor>& lhs, double rhs){
  return Vector<Size>(lhs,1/rhs,
		      Operator<FixedVdMult<Size,Accessor> >());
}

template <class Accessor>
inline Vector<> operator/(const DynamicVector<Accessor>& lhs, double rhs){
  return Vector<>(lhs,1/rhs,
		  Operator<DynamicVdMult<Accessor> >());
}


//////////////////////
// operator /=      //
// Vector /= double //
//////////////////////

template <int Size, class Accessor>
inline FixedVector<Size,Accessor>& operator/=(FixedVector<Size,Accessor>& lhs, const double& rhs){
  return lhs*=(1/rhs);
}

template <class Accessor>
inline DynamicVector<Accessor>& operator/=(DynamicVector<Accessor>& lhs, const double& rhs){
  return lhs*=(1/rhs);
}

///////////////////
//  operator *   //
// (dot product) //
///////////////////

template <int Size, class Accessor1, class Accessor2>
inline double operator*(const FixedVector<Size,Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs) {
    return util::Dot<0,Size-1>::eval(lhs,rhs);
}

template <class Accessor1,class Accessor2>
inline double operator*(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(lhs.size() == rhs.size());
  double dot=0;
  for(int i=0; i<lhs.size(); i++){
      dot += lhs[i]*rhs[i];
  }
  return dot;
}

template <int Size, class Accessor1,class Accessor2>
inline double operator*(const FixedVector<Size,Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(lhs.size() == rhs.size());
  return util::Dot<0,Size-1>::eval(lhs,rhs);
}

template <class Accessor1,int Size,class Accessor2>
inline double operator*(const DynamicVector<Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs){
  assert(lhs.size() == rhs.size());
  return util::Dot<0,Size-1>::eval(lhs,rhs);
}

///////////////////////////
// operator ^            //
// Vector<3> ^ Vector<3> //
// cross product         //
///////////////////////////

template<class Acc1, class Acc2>
struct VXProduct {
  inline static void eval(Vector<3>& ret, const FixedVector<3,Acc1>& lhs, const FixedVector<3,Acc2>& rhs){
    ret[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
    ret[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
    ret[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
  }
};

template<class Acc1, class Acc2>
Vector<3> operator^(const FixedVector<3,Acc1>& lhs, const FixedVector<3,Acc2>& rhs){
  return Vector<3>(lhs,rhs,Operator<VXProduct<Acc1,Acc2> >());
}



/////////////////////
// operator *      //
// Matrix * Vector //
/////////////////////

template <int Rows, int Cols, class LHS, class RHS>
struct FixedMVMult {
  inline static void eval(Vector<Rows>& ret, const LHS& lhs, const RHS& rhs){
    for(int r=0; r<Rows; r++){
	ret[r]=lhs[r]*rhs;
    }
  }
};

template <int Rows, int Cols, class MAccessor, class VAccessor>
inline Vector<Rows> operator*(const FixedMatrix<Rows, Cols, MAccessor>& lhs,
			      const FixedVector<Cols, VAccessor>& rhs){
  return Vector<Rows>(lhs,rhs,
		      Operator<FixedMVMult<Rows,Cols,
		      FixedMatrix<Rows, Cols, MAccessor>,
		      FixedVector<Cols, VAccessor> > >());
}


template <class LHS, class RHS>
struct DynamicMVMult : public VSizer {
  inline static void eval(Vector<>& ret, const LHS& lhs, const RHS& rhs){
    set_size(ret,lhs.num_rows());
    for(int r=0; r<lhs.num_rows(); r++){
      ret[r]=0;
      for(int c=0; c<lhs.num_cols(); c++){
	ret[r]+=lhs(r,c)*rhs[c];
      }
    }
  }
};

template <int Rows, int Cols, class MAccessor, class VAccessor>
inline Vector<Rows> operator*(const FixedMatrix<Rows,Cols,MAccessor>& lhs,
			  const DynamicVector<VAccessor>& rhs){
  assert(rhs.size() == Cols);
  return Vector<Rows>(lhs,rhs,
		  Operator<FixedMVMult<Rows, Cols,
		  FixedMatrix<Rows,Cols,MAccessor>,
		  DynamicVector<VAccessor> > >());

}


template <int Size, class MAccessor, class VAccessor>
inline Vector<> operator*(const DynamicMatrix<MAccessor>& lhs,
			  const FixedVector<Size,VAccessor>& rhs){
  assert(lhs.num_cols() == Size);
  return Vector<>(lhs,rhs,
		  Operator<DynamicMVMult<
		  DynamicMatrix<MAccessor>,
		  FixedVector<Size,VAccessor> > >());
}

template <class MAccessor, class VAccessor>
inline Vector<> operator*(const DynamicMatrix<MAccessor>& lhs,
			  const DynamicVector<VAccessor>& rhs){
  assert(lhs.num_cols() == rhs.size());
  return Vector<>(lhs,rhs,
		  Operator<DynamicMVMult<
		  DynamicMatrix<MAccessor>,
		  DynamicVector<VAccessor> > >());
}


/////////////////////
// operator *      //
// Vector * Matrix //
/////////////////////

template <int Rows, int Cols, class MAccessor, class VAccessor>
inline Vector<Cols> operator*(const FixedVector<Rows, VAccessor>& lhs,
			      const FixedMatrix<Rows, Cols, MAccessor>& rhs){
  return (rhs.T() * lhs);
}



template <int Rows, int Cols, class MAccessor, class VAccessor>
inline Vector<Cols> operator*(const DynamicVector<VAccessor>& lhs,
			  const FixedMatrix<Rows,Cols,MAccessor>& rhs){
  return (rhs.T()*lhs);
}


template <int Size, class MAccessor, class VAccessor>
inline Vector<> operator*(const FixedVector<Size,VAccessor>& lhs,
			  const DynamicMatrix<MAccessor>& rhs){
  return (rhs.T()*lhs);
}

template <class MAccessor, class VAccessor>
inline Vector<> operator*(const DynamicVector<VAccessor>& lhs,
			  const DynamicMatrix<MAccessor>& rhs){
  return (rhs.T()*lhs);
}




/////////////////////
// operator *      //
// Matrix * double //
// (and v.v.)      //
/////////////////////

template <int Rows, int Cols, class Accessor>
struct FixedMdMult {
  inline static void eval(Matrix<Rows,Cols>& ret, const FixedMatrix<Rows,Cols,Accessor>& lhs, double rhs){
    for(int r=0; r<Rows; r++){
      for(int c=0; c<Cols; c++){
	ret[r][c]=lhs[r][c]*rhs;
      }
    }
  }
};

template <int Rows, int Cols, class Accessor>
inline Matrix<Rows,Cols> operator*(const FixedMatrix<Rows,Cols,Accessor>& lhs, double rhs){
  return Matrix<Rows,Cols>(lhs,rhs, Operator<FixedMdMult<Rows,Cols,Accessor> >());
}

template <int Rows, int Cols, class Accessor>
inline Matrix<Rows,Cols> operator*(double lhs, const FixedMatrix<Rows,Cols,Accessor>& rhs){
  return Matrix<Rows,Cols>(rhs,lhs, Operator<FixedMdMult<Rows,Cols,Accessor> >());
}


template <class Accessor>
struct DynamicMdMult : public MSizer {
  inline static void eval(Matrix<>& ret, const DynamicMatrix<Accessor>& lhs, double rhs){
    set_size(ret,lhs.num_rows(),lhs.num_cols());
    for(int r=0; r<lhs.num_rows(); r++){
      for(int c=0; c<lhs.num_cols(); c++){
	ret[r][c]=lhs[r][c]*rhs;
      }
    }
  }
};


template <class Accessor>
inline Matrix<> operator*(const DynamicMatrix<Accessor>& lhs, double rhs){
  return Matrix<>(lhs,rhs, Operator<DynamicMdMult<Accessor> >());
}

template <class Accessor>
inline Matrix<> operator*(double lhs, const DynamicMatrix<Accessor>& rhs){
  return Matrix<>(rhs,lhs, Operator<DynamicMdMult<Accessor> >());
}


//////////////////////
// operator *=      //
// Matrix *= double //
//////////////////////


template <class Accessor>
inline void operator*=(MatrixBase<Accessor>& lhs, double rhs){
  for(int r=0; r<lhs.num_rows(); r++){
    for(int c=0; c<lhs.num_cols(); c++){
      lhs[r][c]*=rhs;
    }
  }
}

/////////////////////
// operator /      //
// Matrix / double //
/////////////////////

template <int Rows, int Cols, class Accessor>
inline Matrix<Rows,Cols> operator/(const FixedMatrix<Rows,Cols,Accessor>& lhs, double rhs){
  return lhs * (1.0/rhs);
}

template <class Accessor>
inline Matrix<> operator/(const DynamicMatrix<Accessor>& lhs, double rhs){
  return lhs * (1.0/rhs);
}


//////////////////////
// operator /=      //
// Matrix /= double //
//////////////////////

template <class Accessor>
inline void operator/=(MatrixBase<Accessor>& lhs, double rhs){
  lhs*=(1.0/rhs);
}

/////////////////////
// operator +      //
// Matrix + Matrix //
/////////////////////

template <int Rows, int Cols, class LHAccessor, class RHAccessor>
struct FixedMMSum {
  inline static void eval(Matrix<Rows,Cols>& ret,
			  const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
			  const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
    for(int r=0; r<Rows; r++){
      for(int c=0; c<Cols; c++){
	ret[r][c]= lhs[r][c]+rhs[r][c];
      }
    }
  }
};

template <int Rows, int Cols, class LHAccessor, class RHAccessor>
Matrix<Rows,Cols> operator+(const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
			    const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
  return Matrix<Rows,Cols>(lhs,rhs, Operator<FixedMMSum<Rows,Cols,LHAccessor,RHAccessor> >());
}

template <class LHS, class RHS>
struct DynamicMMSum : public MSizer {
  inline static void eval(Matrix<>& ret, const LHS& lhs, const RHS& rhs){
    set_size(ret,lhs.num_rows(),lhs.num_cols());
    for(int r=0; r<lhs.num_rows(); r++){
      for(int c=0; c<lhs.num_cols(); c++){
	ret[r][c]= lhs[r][c]+rhs[r][c];
      }
    }
  }
};


template <int Rows, int Cols, class LHAccessor, class RHAccessor>
Matrix<> operator+(const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
		   const DynamicMatrix<RHAccessor>& rhs){
  assert(rhs.num_rows()==Rows && rhs.num_cols()==Cols);
  return Matrix<> (lhs,rhs,
		   Operator<DynamicMMSum<
		   FixedMatrix<Rows,Cols,LHAccessor>,
		   DynamicMatrix<RHAccessor> > >());
}


template <int Rows, int Cols, class LHAccessor, class RHAccessor>
Matrix<> operator+(const DynamicMatrix<LHAccessor>& lhs,
		   const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
  assert(lhs.num_rows()==Rows && lhs.num_cols()==Cols);
  return Matrix<> (lhs,rhs,
		   Operator<DynamicMMSum<
		   DynamicMatrix<LHAccessor>,
		   FixedMatrix<Rows,Cols,RHAccessor> > >());
}

template <class LHAccessor, class RHAccessor>
Matrix<> operator+(const DynamicMatrix<LHAccessor>& lhs,
		   const DynamicMatrix<RHAccessor>& rhs){
  assert(lhs.num_rows()==rhs.num_rows() && lhs.num_cols()==rhs.num_cols());
  return Matrix<> (lhs,rhs,
		   Operator<DynamicMMSum<
		   DynamicMatrix<LHAccessor>,
		   DynamicMatrix<RHAccessor> > >());
}

//////////////////////
// operator +=, -=  //
// Matrix += Matrix //
//////////////////////

//fixed fixed
//fixed dynamic
//dynamic fixed
//dynamic dynamic
//RefCM fixed
//RefCM dynamic
//RefRM fixed
//RefRM dynamic

#define TOON_MAKE_ELEMENT_OPS(OP) \
 \
template <int Rows, int Cols, class MAccessor1, class MAccessor2> inline FixedMatrix<Rows,Cols,MAccessor1>& operator OP ( FixedMatrix<Rows,Cols,MAccessor1>& lhs, const FixedMatrix<Rows,Cols,MAccessor2>& rhs){ \
    for(int r=0; r<Rows; r++) \
		lhs[r] OP rhs[r]; \
    return lhs; \
} \
\
template <int Rows, int Cols, class MAccessor1, class MAccessor2> inline FixedMatrix<Rows,Cols,MAccessor1>& operator OP ( FixedMatrix<Rows,Cols,MAccessor1>& lhs, const DynamicMatrix<MAccessor2>& rhs){ \
  assert(rhs.num_rows() == Rows && rhs.num_cols() == Cols); \
  for(int r=0; r<Rows; r++) \
      lhs[r] OP rhs[r]; \
  return lhs; \
} \
 \
template <int Rows, int Cols, class MAccessor1, class MAccessor2> inline DynamicMatrix<MAccessor1>& operator OP ( DynamicMatrix<MAccessor1>& lhs, const FixedMatrix<Rows,Cols,MAccessor2>& rhs){ \
  assert(lhs.num_rows() == Rows && lhs.num_cols() == Cols); \
  for(int r=0; r<Rows; r++){ \
      lhs[r] OP rhs[r]; \
  } \
  return lhs; \
} \
 \
template <class MAccessor1, class MAccessor2> inline DynamicMatrix<MAccessor1>& operator OP ( DynamicMatrix<MAccessor1>& lhs, const DynamicMatrix<MAccessor2>& rhs){ \
  assert(lhs.num_rows() == rhs.num_rows() && lhs.num_cols() == rhs.num_cols()); \
  for(int r=0; r<lhs.num_rows(); r++) \
    for(int c=0; c <lhs.num_cols(); c++) \
      lhs[r][c] OP rhs[r][c]; \
  return lhs; \
} \
 \
template <int Rows, int Cols, class MAccessor> inline RefSkipMatrixRM operator OP (RefSkipMatrixRM lhs, const FixedMatrix<Rows,Cols,MAccessor>& rhs) \
{ \
  assert(lhs.num_rows() == Rows && lhs.num_cols() == Cols); \
  for(int r=0; r<Rows; r++) \
	  for(int c=0; c<Cols; c++) \
		  lhs[r][c] OP rhs[r][c]; \
   \
  return lhs; \
} \
 \
template <class MAccessor> inline RefSkipMatrixRM operator OP (RefSkipMatrixRM lhs, const DynamicMatrix<MAccessor>& rhs) \
{ \
  assert(lhs.num_rows() == rhs.num_rows() && lhs.num_cols() == rhs.num_cols()); \
  for(int r=0; r < lhs.num_rows(); r++) \
	  for(int c=0; c< lhs.num_cols(); c++) \
      lhs[r][c] OP rhs[r][c]; \
   \
  return lhs; \
} \
 \
 \
template <int Rows, int Cols, class MAccessor> inline RefSkipMatrixCM operator OP (RefSkipMatrixCM lhs, const FixedMatrix<Rows,Cols,MAccessor>& rhs) { \
  assert(lhs.num_rows() == Rows && lhs.num_cols() == Cols); \
  for(int r=0; r<Rows; r++) \
	  for(int c=0; c<Cols; c++) \
      lhs[r][c] OP rhs[r][c]; \
   \
  return lhs; \
} \
 \
template <class MAccessor> inline RefSkipMatrixCM operator OP (RefSkipMatrixCM lhs, const DynamicMatrix<MAccessor>& rhs) \
{ \
  assert(lhs.num_rows() == rhs.num_rows() && lhs.num_cols() == rhs.num_cols()); \
  for(int r=0; r < lhs.num_rows(); r++) \
	  for(int c=0; c< lhs.num_cols(); c++) \
      lhs[r][c] OP rhs[r][c]; \
   \
  return lhs; \
} 

TOON_MAKE_ELEMENT_OPS(+=)
TOON_MAKE_ELEMENT_OPS(-=)





/////////////////////
// operator -      //
// Matrix - Matrix //
/////////////////////

template <int Rows, int Cols, class LHAccessor, class RHAccessor>
struct FixedMMSub {
  inline static void eval(Matrix<Rows,Cols>& ret,
			  const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
			  const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
    for(int r=0; r<Rows; r++){
      for(int c=0; c<Cols; c++){
	ret[r][c]= lhs[r][c]-rhs[r][c];
      }
    }
  }
};

template <int Rows, int Cols, class LHAccessor, class RHAccessor> inline
Matrix<Rows,Cols> operator-(const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
			    const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
  return Matrix<Rows,Cols>(lhs,rhs, Operator<FixedMMSub<Rows,Cols,LHAccessor,RHAccessor> >());
}

template <class LHS, class RHS>
struct DynamicMMSub : public MSizer {
  inline static void eval(Matrix<>& ret, const LHS& lhs, const RHS& rhs){
    set_size(ret,lhs.num_rows(),lhs.num_cols());
    for(int r=0; r<lhs.num_rows(); r++){
      for(int c=0; c<lhs.num_cols(); c++){
	ret[r][c]= lhs[r][c]-rhs[r][c];
      }
    }
  }
};


template <int Rows, int Cols, class LHAccessor, class RHAccessor>
Matrix<Rows, Cols> operator-(const FixedMatrix<Rows,Cols,LHAccessor>& lhs,
		   const DynamicMatrix<RHAccessor>& rhs){
  assert(rhs.num_rows()==Rows && rhs.num_cols()==Cols);
  return Matrix<Rows, Cols> (lhs,rhs,
		   Operator<FixedMMSub<Rows, Cols,
		   FixedMatrix<Rows,Cols,LHAccessor>,
		   DynamicMatrix<RHAccessor> > >());
}


template <int Rows, int Cols, class LHAccessor, class RHAccessor>
Matrix<Rows, Cols> operator-(const DynamicMatrix<LHAccessor>& lhs,
		   const FixedMatrix<Rows,Cols,RHAccessor>& rhs){
  assert(lhs.num_rows()==Rows && lhs.num_cols()==Cols);
  return Matrix<Rows, Cols> (lhs,rhs,
		   Operator<FixedMMSub<Rows, Cols,
		   DynamicMatrix<LHAccessor>,
		   FixedMatrix<Rows,Cols,RHAccessor> > >());
}

template <class LHAccessor, class RHAccessor>
Matrix<> operator-(const DynamicMatrix<LHAccessor>& lhs,
		   const DynamicMatrix<RHAccessor>& rhs){
  assert(lhs.num_rows()==rhs.num_rows() && lhs.num_cols()==rhs.num_cols());
  return Matrix<> (lhs,rhs,
		   Operator<DynamicMMSub<
		   DynamicMatrix<LHAccessor>,
		   DynamicMatrix<RHAccessor> > >());
}

/////////////////////
// operator *      //
// Matrix * Matrix //
/////////////////////

template <class RET, class LHS, class RHS>
inline void cppmmmult (RET& ret, const LHS& lhs, const RHS& rhs){
  const int Rows = lhs.num_rows();
  const int Cols = rhs.num_cols();
  const int Inter = lhs.num_cols();
  for(int r=0; r<Rows; r++){
    for(int c=0; c<Cols; c++){
      double temp=0;
      for(int i=0; i<Inter; i++){
	temp+=lhs(r,i)*rhs(i,c);
      }
      ret(r,c)=temp;
    }
  }
}


template <int Rows, int Inter, int Cols, class LHS, class RHS, int MultPolicy>
struct FixedMMMult;

template <int Rows, int Inter, int Cols, class LHS, class RHS>
struct FixedMMMult<Rows,Inter,Cols,LHS,RHS,NUMERICS::BlasMult> {
  inline static void eval(Matrix<Rows,Cols>& ret, const LHS& lhs, const RHS& rhs){

    //blasmmmult(ret,lhs,rhs);
      //cppmmmult(ret, lhs, rhs);
      util::matrix_multiply<Rows,Inter,Cols>(lhs, rhs, ret);
  }
};

template <int Rows, int Inter, int Cols, class LHS, class RHS>
struct FixedMMMult<Rows,Inter,Cols,LHS,RHS,NUMERICS::CPPMult> {
  inline static void eval(Matrix<Rows,Cols>& ret, const LHS& lhs, const RHS& rhs){
      //cppmmmult(ret,lhs,rhs);
      util::matrix_multiply<Rows,Inter,Cols>(lhs, rhs, ret);
  }
};


template <int Rows, int Inter, int Cols, class LMAccessor, class RMAccessor>
inline Matrix<Rows,Cols> operator*(const FixedMatrix<Rows,Inter,LMAccessor>& lhs,
				   const FixedMatrix<Inter,Cols,RMAccessor>& rhs){
  return Matrix<Rows,Cols>(lhs,rhs,
			   Operator<FixedMMMult<Rows,Inter,Cols,
			   FixedMatrix<Rows,Inter,LMAccessor>,
			   FixedMatrix<Inter,Cols,RMAccessor>,
			   (Rows*Inter*Cols>NUMERICS::MaxCPPMultCount ? NUMERICS::BlasMult : NUMERICS::CPPMult) > >());
}


template<class LHS, class RHS>
struct DynamicMMMult : public MSizer{
  inline static void eval(Matrix<>& ret, const LHS& lhs, const RHS& rhs){
    set_size(ret,lhs.num_rows(),rhs.num_cols());
    //if(lhs.num_rows()*lhs.num_cols()*rhs.num_cols() > NUMERICS::MaxCPPMultCount){
     // blasmmmult(ret,lhs,rhs);
    //} else {
      cppmmmult(ret,lhs,rhs);
    //}
  }
};

template<int Rows, int Cols, class LMAccessor, class RMAccessor>
inline Matrix<> operator*(const FixedMatrix<Rows,Cols,LMAccessor>& lhs,
			  const DynamicMatrix<RMAccessor>& rhs){
  assert(rhs.num_rows()==Cols);
  return Matrix<>(lhs,rhs,
		  Operator<DynamicMMMult<
		  FixedMatrix<Rows,Cols,LMAccessor>,
		  DynamicMatrix<RMAccessor> > > ());
}


template<int Rows, int Cols, class LMAccessor, class RMAccessor>
inline Matrix<> operator*(const DynamicMatrix<LMAccessor>& lhs,
			  const FixedMatrix<Rows,Cols,RMAccessor>& rhs){
  assert(lhs.num_cols()==Rows);
  return Matrix<>(lhs,rhs,
		  Operator<DynamicMMMult<
		  DynamicMatrix<LMAccessor>,
		  FixedMatrix<Rows,Cols,RMAccessor> > > ());
}

template <class LMAccessor, class RMAccessor>
inline Matrix<> operator*(const DynamicMatrix<LMAccessor>& lhs,
			  const DynamicMatrix<RMAccessor>& rhs){
  assert(lhs.num_cols() == rhs.num_rows());
  return Matrix<>(lhs,rhs,
		  Operator<DynamicMMMult<
		  DynamicMatrix<LMAccessor>,
		  DynamicMatrix<RMAccessor> > >());
}



///////////////////////////////////////////////////////////////////////////////////////
// Multiplication treating a vector as a diagonal matrix

/////////////////////////
//                     //
// diagmult()          //
// Vector * Vector     //
//                     //
/////////////////////////

template <class LHS, class RHS, int Size>
struct DiagVVMult {
  inline static void eval(Vector<Size>& ret, const LHS& lhs, const RHS& rhs){
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i]*rhs[i];
    }
  }
};

template <class Accessor1, class Accessor2, int Size>
Vector<Size> diagmult(const FixedVector<Size,Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs){
  return Vector<Size>(lhs,rhs,Operator<DiagVVMult<FixedVector<Size,Accessor1>,FixedVector<Size,Accessor2>,Size> >());
}


template<class LHS, class RHS>
struct DynamicDiagVVMult : public VSizer {
  inline static void eval(Vector<>& ret, const LHS& lhs, const RHS& rhs){
    assert(lhs.size() == rhs.size());
    const int Size = lhs.size();
    set_size(ret,Size);
    for(int i=0; i<Size; i++){
      ret[i] = lhs[i] * rhs[i];
    }
  }
};

template<class Accessor1, class Accessor2>
Vector<> diagmult(const DynamicVector<Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  return Vector<>(lhs,rhs,Operator<DynamicDiagVVMult<DynamicVector<Accessor1>,DynamicVector<Accessor2> > >());
}

template<class Accessor1, class Accessor2, int Size>
Vector<Size> diagmult(const FixedVector<Size,Accessor1>& lhs, const DynamicVector<Accessor2>& rhs){
  assert(lhs.size() == rhs.size());
  return Vector<Size>(lhs,rhs,Operator<DiagVVMult<FixedVector<Size,Accessor1>,DynamicVector<Accessor2>, Size > >());
}

template<class Accessor1, class Accessor2, int Size>
Vector<Size> diagmult(const DynamicVector<Accessor1>& lhs, const FixedVector<Size,Accessor2>& rhs){
  assert(lhs.size() == rhs.size());
  return Vector<Size>(lhs,rhs,Operator<DiagVVMult<DynamicVector<Accessor1>,FixedVector<Size,Accessor2>, Size > >());
}

/////////////////////////
//                     //
// diagmult()          //
// Vector * Matrix     //
//                     //
/////////////////////////


template <class LHS, class RHS, int Rows, int Cols>
struct DiagVMMult {
  inline static void eval(Matrix<Rows,Cols>& ret, const LHS& lhs, const RHS& rhs){
    for(int c=0; c<Cols; c++){
      ret.T()[c] = diagmult(lhs,rhs.T()[c]);
    }
  }
};

template <class Accessor1, class Accessor2, int Rows, int Cols>
inline Matrix<Rows,Cols> diagmult(const FixedVector<Rows, Accessor1>& lhs, const FixedMatrix<Rows,Cols,Accessor2>& rhs){
  return Matrix<Rows,Cols> (lhs,rhs,Operator<DiagVMMult<FixedVector<Rows, Accessor1>,FixedMatrix<Rows,Cols,Accessor2>,Rows,Cols> >());
}



template<class LHS, class RHS>
struct DynamicDiagVMMult : public MSizer {
  inline static void eval(Matrix<>& ret, const LHS& lhs, const RHS& rhs){
    assert(lhs.size() == rhs.num_rows());
    const int Rows = rhs.num_rows();
    const int Cols = rhs.num_cols();
    set_size(ret,Rows,Cols);
    for(int r=0; r<Rows; r++){
      for(int c=0; c<Cols; c++){
	ret(r,c) = lhs[r] * rhs(r,c);
      }
    }
  }
};

template<class Accessor1, class Accessor2>
inline Matrix<> diagmult(const DynamicVector<Accessor1>& lhs, const DynamicMatrix<Accessor2>& rhs){
  return Matrix<>(lhs,rhs,Operator<DynamicDiagVMMult<DynamicVector<Accessor1>,DynamicMatrix<Accessor2> > >());
}

template<class Accessor1, class Accessor2, int Size>
inline Matrix<> diagmult(const FixedVector<Size,Accessor1>& lhs, const DynamicMatrix<Accessor2>& rhs){
  return Matrix<>(lhs,rhs,Operator<DynamicDiagVMMult<FixedVector<Size,Accessor1>,DynamicMatrix<Accessor2> > >());
}

template<class Accessor1, class Accessor2, int Rows, int Cols>
inline Matrix<Rows, Cols> diagmult(const DynamicVector<Accessor1>& lhs, const FixedMatrix<Rows,Cols,Accessor2>& rhs){
  assert(lhs.size() == Rows);
  return Matrix<Rows, Cols>(lhs,rhs,Operator<DiagVMMult<DynamicVector<Accessor1>,FixedMatrix<Rows,Cols,Accessor2>,Rows,Cols > >());
}

/////////////////////////
//                     //
// diagmult()          //
// Matrix * Vector     //
//                     //
/////////////////////////


template <class LHS, class RHS, int Rows, int Cols>
struct DiagMVMult {
  inline static void eval(Matrix<Rows,Cols>& ret, const LHS& lhs, const RHS& rhs){
    for(int r=0; r<Rows; r++){
      ret[r] = diagmult(lhs[r],rhs);
    }
  }
};

template <class Accessor1, class Accessor2, int Rows, int Cols>
inline Matrix<Rows,Cols> diagmult(const FixedMatrix<Rows,Cols,Accessor1>& lhs, const FixedVector<Cols,Accessor2>& rhs){
  return Matrix<Rows,Cols> (lhs,rhs,Operator<DiagMVMult<FixedMatrix<Rows,Cols,Accessor1>,FixedVector<Cols,Accessor2>,Rows,Cols > >());
}



template<class LHS, class RHS>
struct DynamicDiagMVMult : public MSizer {
  inline static void eval(Matrix<>& ret, const LHS& lhs, const RHS& rhs){
    assert(rhs.size() == lhs.num_cols());
    const int Rows = lhs.num_rows();
    const int Cols = lhs.num_cols();
    set_size(ret,Rows,Cols);
    for(int r=0; r<Rows; r++){
      for(int c=0; c<Cols; c++){
	ret(r,c) = lhs(r,c) * rhs[c];
      }
    }
  }
};

template<class Accessor1, class Accessor2>
inline Matrix<> diagmult(const DynamicMatrix<Accessor1>& lhs,  const DynamicVector<Accessor2>& rhs){
  return Matrix<>(lhs,rhs,Operator<DynamicDiagMVMult<DynamicMatrix<Accessor1>, DynamicVector<Accessor2> > >());
}

template<class Accessor1, class Accessor2, int Size>
inline Matrix<> diagmult(const DynamicMatrix<Accessor1>& lhs,  const FixedVector<Size,Accessor2>& rhs){
  return Matrix<>(lhs,rhs,Operator<DynamicDiagMVMult<DynamicMatrix<Accessor1>, FixedVector<Size,Accessor2> > >());
}

template<class Accessor1, class Accessor2, int Rows, int Cols>
inline Matrix<Rows, Cols> diagmult(const FixedMatrix<Rows,Cols,Accessor1>& lhs,  const DynamicVector<Accessor2>& rhs){
  assert(rhs.size() == Cols);
  return Matrix<Rows, Cols>(lhs,rhs,Operator<DiagMVMult<FixedMatrix<Rows,Cols,Accessor1>, DynamicVector<Accessor2>,Rows,Cols > >());
}

#endif
