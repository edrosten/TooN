
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
#ifndef __VACCESSOR_HH
#define __VACCESSOR_HH


///////////////////////////////
//                           //
//  Vector Accessor classes  //
//                           //
///////////////////////////////

template <int Size, class AllocZone> class FixedVAccessor;

/////////////  DYNAMIC SIZED ACCESSORS ////////////////

class DynamicVAccessor {
  friend struct VSizer;
public:
  typedef DynamicVector<DynamicVAccessor> RefVector;
    
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

  inline RefVector slice(int start, int size) 
  {
    assert(0 <= start && start < this->my_size && size >=0 && start+size <= this->my_size);
    RefVector ret;
    ret.set(size,this->my_values + start);
    return ret;
  }

  inline const RefVector slice(int start, int size) const 
  {
    assert(0 <= start && start < this->my_size && size >=0 && start+size <= this->my_size);
    RefVector ret;
    ret.set(size,this->my_values + start);
    return ret;
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

  inline DynamicMatrix<DynamicMAccessor<RowMajor> > as_row(); // implemented in linoperators.hh
  inline DynamicMatrix<DynamicMAccessor<ColMajor> > as_col(); //
  inline void set(int size, double* values) { my_size = size;  my_values = values; }

  typedef double* iterator;
  typedef const double* const_iterator;

  iterator begin()
  {
  	return my_values;
  }

  iterator end()
  {
	return my_values + my_size;
  }

  const_iterator begin() const
  {
	return my_values;
  }	
  const_iterator end() const
  {
	return my_values + my_size;
  }
  
 protected:  
  int my_size;
  double* my_values;
};

typedef DynamicVAccessor::RefVector RefVector;
inline RefVector makeRefVector(int size, double* values) { RefVector ret; ret.set(size,values); return ret; }

class DynamicSkipAccessor{
public:
  typedef DynamicVector<DynamicSkipAccessor> RefSkipVector;
  
  //CHECK THIS
  template<int Start, int Length>
  inline DynamicVector<DynamicSkipAccessor> slice() 
  {
    //DynamicSkipAccessors do not own memory, so destruction of one will not free it 
	
	DynamicVector<DynamicSkipAccessor> r;

	r.my_size = Length;
	r.my_skip = my_skip;
	r.my_values = my_values + my_skip * Start;
	
	return r;
  }

  template<int Start, int Length>
  inline const DynamicVector<DynamicSkipAccessor> slice() const 
  {
    //DynamicSkipAccessors do not own memory, so destruction of one will not free it 
	
	DynamicVector<DynamicSkipAccessor> r;
	r.my_size = Length;
	r.my_skip = my_skip;
	r.my_values = my_values + my_skip * Start;
	
	return r;
  }

  inline RefSkipVector slice(int start, int size) 
  {
    assert(0 <= start && start < this->my_size && size >=0 && start+size <= this->my_size);
    RefSkipVector ret;
    ret.set(size,my_skip,this->my_values + start*my_skip);
    return ret;
  }

  inline const RefSkipVector slice(int start, int size) const 
  {
    assert(0 <= start && start < this->my_size && size >=0 && start+size <= this->my_size);
    RefSkipVector ret;
    ret.set(size,my_skip,this->my_values + start*my_skip);
    return ret;
  }

  inline const double& operator[] (int i)const 
  {
    return my_values[i*my_skip];
  }

  inline double& operator[] (int i) 
  {
    return my_values[i*my_skip];	
  }
  
  inline int size() const 
  {
    return my_size;
  }

  inline DynamicMatrix<RefSkipMAccessor<ColMajor> > as_row(); // implemented in linoperators.hh
  inline DynamicMatrix<RefSkipMAccessor<RowMajor> > as_col(); //

  inline void set(int size, int skip, double* values) { my_size = size; my_skip = skip; my_values = values; }

  private:
	  template<class T> class It: public std::iterator<std::random_access_iterator_tag, double>
	  {
		private:
		  T* d;
		  int skip;
		public:
		  bool operator!=(const It&i){return d != i.d;}
		  bool operator==(const It&i){return d == i.d;}
		  bool operator<=(const It&i){return d <= i.d;}
		  bool operator>=(const It&i){return d >= i.d;}
		  bool operator<(const It&i){return d < i.d;}
		  bool operator>(const It&i){return d > i.d;}
		  It& operator+=(ptrdiff_t i) {d += i*skip;return *this;}
		  It& operator-=(ptrdiff_t i) {d -= i*skip;return *this;}
		  It& operator++() {d += skip;return *this;}
		  It& operator--() {d -= skip;return *this;}
		  It operator++(int) {T* t=d;d += skip;return iterator(t,skip);}
		  It operator--(int) {T* t=d;d -= skip;return iterator(t,skip);}
		  It operator+(ptrdiff_t i) {return It(d + i*skip,skip);}
		  It operator-(ptrdiff_t i) {return It(d - i*skip,skip);}
		  ptrdiff_t operator-(const It& i) {return (d - i.d)/skip;};
		  T& operator*(){return *d;}
		  It(T*a, int s):d(a),skip(s){}
	  };

  public:

  typedef It<double> iterator;
  typedef It<const double> const_iterator;

  iterator begin()
  {
  	return iterator(my_values, my_skip);
  }

  iterator end()
  {
  	return iterator(my_values, my_skip) + my_size;
  }

  const_iterator begin() const
  {
  	return const_iterator(my_values, my_skip);
  }	
  const_iterator end() const
  {
  	return const_iterator(my_values, my_skip) + my_size;
  }

 protected:
  int my_size;
  int my_skip;
  double* my_values;
};

typedef DynamicSkipAccessor::RefSkipVector RefSkipVector;

inline RefSkipVector makeRefSkipVector(int size, int skip, double* values) { RefSkipVector ret; ret.set(size,skip,values); return ret; }


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
      util::Assert<(Start+Length <= Size)>();
      FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >::dummy();
      return reinterpret_cast<FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >&> (parent::my_values[Start]);
  }

  typedef FixedVAccessor<Size,AllocZone> this_type;
  inline RefVector slice(int start, int size) 
  {
    assert(0 <= start && start < Size && size >=0 && start+size <= Size);
    return makeRefVector(size, parent::my_values + start);
  }

  inline const RefVector slice(int start, int size) const
  {
    assert(0 <= start && start < Size && size >=0 && start+size <= Size);
    return makeRefVector(size, const_cast<double*>( parent::my_values + start));
  }
  
  template<int Start, int Length>
  inline const FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >& slice() const 
  {
      util::Assert<(Start+Length <= Size)>();
    return reinterpret_cast<const FixedVector<Length,FixedVAccessor<Length,Stack<Length> > >&> (parent::my_values[Start]);
  }

  // convert to Matrices
  inline FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >& as_col() 
  {
      FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >::dummy();
    return reinterpret_cast<FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >&>(*parent::my_values);
  }

  inline const FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >& as_col() const 
  {
      FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >::dummy();
    return reinterpret_cast<const FixedMatrix<Size,1,FixedMAccessor<Size,1,ColMajor,Stack<Size> > >&> (*parent::my_values);
  }

  inline FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >& as_row() 
  {
      FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >::dummy();
    return reinterpret_cast<FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >&> (*parent::my_values);
  }
  
  inline const FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >& as_row() const 
  {
      FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >::dummy();
    return reinterpret_cast<const FixedMatrix<1,Size,FixedMAccessor<1,Size,RowMajor,Stack<Size> > >&> (*parent::my_values);
  }
  
  typedef double* iterator;
  typedef const double* const_iterator;

  iterator begin()
  {
  	return parent::my_values;
  }

  iterator end()
  {
	return parent::my_values + Size;
  }

  const_iterator begin() const
  {
	return parent::my_values;
  }	
  const_iterator end() const
  {
	return parent::my_values + Size;
  }

};


template <int Size, int Skip>
class SkipAccessor : public Stack<Size*Skip>{
 public:
  typedef Stack<Size*Skip> parent;
  typedef SkipAccessor<Size,Skip> this_type;
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
  inline FixedVector<Length, SkipAccessor<Size, Skip> >& slice() 
  {
      util::Assert<(Start+Length <= Size)>();
      return reinterpret_cast<FixedVector<Length, SkipAccessor<Size, Skip> >&>(parent::my_values[Start*Skip]);
  }

  template<int Start, int Length>
  inline const FixedVector<Length, SkipAccessor<Size, Skip> >& slice() const 
  {
      util::Assert<(Start+Length <= Size)>();
    return reinterpret_cast<const FixedVector<Length, SkipAccessor<Size, Skip> >&>(parent::my_values[Start*Skip]);
  }

  RefSkipVector slice(int start, int size) 
  {
    assert(0 <= start && start < Size && size >=0 && start+size <= Size);
    return makeRefSkipVector(size, Skip, parent::my_values + start*Skip);
  }
  const RefSkipVector slice(int start, int size) const 
  {
    assert(0 <= start && start < Size && size >=0 && start+size <= Size);
    return makeRefSkipVector(size, Skip, parent::my_values + start*Skip);
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

  private:
	  template<class T> class It: public std::iterator<std::random_access_iterator_tag, double>
	  {
		private:
		  T* d;
		public:
		  bool operator!=(const It&i){return d != i.d;}
		  bool operator==(const It&i){return d == i.d;}
		  bool operator<=(const It&i){return d <= i.d;}
		  bool operator>=(const It&i){return d >= i.d;}
		  bool operator<(const It&i){return d < i.d;}
		  bool operator>(const It&i){return d > i.d;}
		  It& operator+=(ptrdiff_t i) {d += i*Skip;return *this;}
		  It& operator-=(ptrdiff_t i) {d -= i*Skip;return *this;}
		  It& operator++() {d += Skip;return *this;}
		  It& operator--() {d -= Skip;return *this;}
		  It operator++(int) {T* t=d;d += Skip;return iterator(t);}
		  It operator--(int) {T* t=d;d -= Skip;return iterator(t);}
		  It operator+(ptrdiff_t i) {return It(d + i*Skip);}
		  It operator-(ptrdiff_t i) {return It(d - i*Skip);}
		  ptrdiff_t operator-(const It& i) {return (d - i.d)/Skip;};
		  T& operator*(){return *d;}
		  It(T*a):d(a){}
	  };

  public:

  typedef It<double> iterator;
  typedef It<const double> const_iterator;

  iterator begin()
  {
  	return iterator(parent::my_values);
  }

  iterator end()
  {
	return iterator(parent::my_values) + Size;
  }

  const_iterator begin() const
  {
	return iterator(parent::my_values);
  }	
  const_iterator end() const
  {
	return iterator(parent::my_values) + Size;
  }

};


#endif
