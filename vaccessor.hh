
/*                       
	 Copyright (C) 2005 The Authors

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
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __VACCESSOR_HH
#define __VACCESSOR_HH


///////////////////////////////
//                           //
//  Vector Accessor classes  //
//                           //
///////////////////////////////

/////////////  FIXED SIZED ACCESSORS ////////////////

template <int Size, class AllocZone>
class FixedVAccessor : public AllocZone {
 public:
  typedef AllocZone parent;
  inline const double& operator[] (int i)const 
  {
    return parent::my_values[i];
  }

  inline double& operator[] (int i) 
  {
    return parent::my_values[i];
  }


  inline static int size() {return Size;}

  template<int Start, int Length>
  inline FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >& slice() 
  {
    return reinterpret_cast<FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >&> (parent::my_values[Start]);
  }

  
  template<int Start, int Length>
  inline const FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >& slice() const 
  {
    return reinterpret_cast<const FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >&> (parent::my_values[Start]);
  }

  // convert to Matrices
  inline FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >& as_col() 
  {
    return reinterpret_cast<FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >&>(*parent::my_values);
  }

  inline const FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >& as_col() const 
  {
    return reinterpret_cast<const FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >&> (*parent::my_values);
  }

  inline FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >& as_row() 
  {
    return reinterpret_cast<FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >&> (*parent::my_values);
  }
  
  inline const FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >& as_row() const 
  {
    return reinterpret_cast<const FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >&> (*parent::my_values);
  }
  
};


template <int Size, int Skip>
class SkipAccessor : public Stack<Size*Skip>{
 public:
  typedef Stack<Size*Skip> parent;
  inline const double& operator[] (int i) const  
  {
    return parent::my_values[i*Skip];
  }

  inline double& operator[] (int i) 
  {
	return parent::my_values[i*Skip];
  }

  inline static int size() 
  {
    return Size;
  }

  template<int Start, int Length>
  inline FixedVector<Size, SkipAccessor<Size, Skip> >& slice() 
  {
	return reinterpret_cast<FixedVector<Size, SkipAccessor<Size, Skip> >&>(parent::my_values[Start*Skip]);
  }

  template<int Start, int Length>
  inline const FixedVector<Size, SkipAccessor<Size, Skip> >& slice() const 
  {
    return reinterpret_cast<const FixedVector<Size, SkipAccessor<Size, Skip> >&>(parent::my_values[Start*Skip]);
  }

  // convert to Matrices
  inline FixedMatrix<Size,1,SkipMAccessor<Size,1,Skip,RowMajor> >& as_col() 
  {
    return reinterpret_cast<FixedMatrix<Size,1,SkipMAccessor<Size,1,Skip,RowMajor> >&>(*parent::my_values);
  }

  inline const FixedMatrix<Size,1,SkipMAccessor<Size,1,Skip,RowMajor> >& as_col() const 
  {
    return reinterpret_cast<FixedMatrix<Size,1,SkipMAccessor<Size,1,Skip,RowMajor> >&>(*parent::my_values);
  }

  inline FixedMatrix<1,Size,SkipMAccessor<1,Size,Skip,ColMajor> >& as_row() 
  {
    return reinterpret_cast<FixedMatrix<1,Size,SkipMAccessor<1,Size,Skip,ColMajor> >&>(*parent::my_values);
  }

  inline const FixedMatrix<1,Size,SkipMAccessor<1,Size,Skip,ColMajor> >& as_row() const 
  {
    return reinterpret_cast<const FixedMatrix<1,Size,SkipMAccessor<1,Size,Skip,ColMajor> >&>(*parent::my_values);
  }

};

/////////////  DYNAMIC SIZED ACCESSORS ////////////////


class DynamicVAccessor {
  friend class VSizer;
 public:
  template<int Start, int Length>
  inline FixedVector<Length, FixedVAccessor<Length, Stack<Length> > >& slice() 
  {
    return reinterpret_cast<FixedVector<Length, FixedVAccessor<Length, Stack<Length> > >&> (my_values[Start]);
  }

  template<int Start, int Length>
  inline const FixedVector<Length, FixedVAccessor<Length, Stack<Length> > >& slice() const 
  {
    return reinterpret_cast<const FixedVector<Length, FixedVAccessor<Length, Stack<Length> > >&> (my_values[Start]);
  }

  inline const double& operator[] (int i)const 
  {
    return my_values[i];
  }

  inline double& operator[] (int i) 
  {
    return my_values[i];
  }

  inline int size() const 
  {
    return my_size;
  }

  inline RefMatrix<RowMajor> as_row(); // implemented in linoperators.hh
  inline RefMatrix<ColMajor> as_col(); //

 protected:
  int my_size;
  double* my_values;
};


#endif
