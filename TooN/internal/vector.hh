//-*- c++ -*-
//
// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND OTHER CONTRIBUTORS ``AS IS''
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR OTHER CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.



namespace TooN {


/**
@class TooN::Vector vector.hh TooN/toon.h
A vector. 
Support is provided for all the usual vector operations: 
- elements can be accessed using the familiar [] notation with the index starting at 0
- they can be added or subtracted
- they can be printed or loaded from streams
- they can be multiplied (on either side) or divided by a scalar on the right:
- the vector dot product can be computed
- subvectors can be extracted using the templated slice() member function
- the vector cross product can be computed for statically sized 3-vectors

See individual member function documentation for examples of usage.


\par Statically- and dynamically-sized vectors

The library provides classes for both statically- and
dynamically-sized vectors. If you know what dimension of vector you're
going to use (e.g. 3 to represent a point in 3D space), it's more
efficient to statically sized vectors. The size of static vectors is
determined at compile time; that of dynamically-sized vectors at
run-time.

To create a 3-dimensional vector, use:
@code
Vector<3> v;
@endcode

and to create a vector of some other dimensionality just replace 3
with the positive integer of your choice, or some expression which the
compiler can evaluate to an integer at compile time.

The preferred way of initialising a vector is to use makeVector. The
%makeVector function constructs a static vector initialised to the
size and the contents of the comma-separated list of argments.  The
%makeVector vectors are real Vectors and so can be used anywhere where
a vector is needed, not just in initialisations. For example

@code
// Create a vector initialised to [1 2 3];
Vector<3> v = makeVector(1, 2, 3);
// Calculate the dot product with the vector [4 0 6]
double dot = v * makeVector(4, 0, 6);
@endcode

Because the %make_Vector syntax creates actual vectors, compile-time
checking is done to ensure that all vectors defined in this way have
the correct number of elements.

\par Dynamically-sized vectors

To create a dynamically sized vector, use:
@code
Vector<> v(size);
@endcode
where size is an integer which will be evaluated at run time.

Vector<> is actually a synonym for Vector<Dynamic> which is Vector<-1>
being a template specialisation of Vector<N> with a special
implementation that allows the size to be determined at runtime.


\par Row vectors and column vectors

This library makes no distinction between row vectors and column
vectors. Vectors that appear on the left of a multiplication are
treated as row vectors while those that appear on the right are
treated as column vectors (thus <code>v1*v2</code> means the dot
product). This means that sometimes you have to be careful to include
prarentheses since it is possible to write obscure stuff like

@code
Vector<4> v4 = v1 * v2 * v3;
@endcode

which in the absence of any extra parentheses means 'compute the dot
product between <code>v1</code> and <code>v2</code> and then multiply
<code>v3</code> by this scalar and assign to <code>v4</code>'.

If the row-column distinction is important, then vectors can be turned
into matrices with one row or column by using as_row() or as_col():

@code
double d[3] = {1,2,3};
Vector<3> v(d);
Matrix<3,3> M = v.as_col() * v.as_row(); // creates a symmetric rank 1 matrix from v 
@endcode

@ingroup gLinAlg
**/
template<int Size=Dynamic, typename Precision=DefaultPrecision, typename Base=Internal::VBase>
struct Vector : public Base::template VLayout<Size, Precision> {
protected:
public:
	
	///Value of template Size parameter. Used by internal classes to reduce
	///parameter list sizes.
	///@internal
	static const int SizeParameter = Size;

	typedef typename Base::template VLayout<Size, Precision>::PointerType PointerType;
  // sneaky hack: only one of these constructors will work with any given base
  // class but they don't generate errors unless the user tries to use one of them
  // although the error message may be less than helpful - maybe this can be changed?

	/// @name Constructors
	//@{

	/// Default constructor for vectors.  For fixed-sized vectors,
	/// this does nothing, i.e. does not guarantee to initialise the
	/// vector to any particular values.  For dynamically sized
	/// vectors, this sets the vector to have a length of 0 which
	/// renders the vector useless because vectors can't be resized
	inline Vector(){}

	Vector(Vector&&) = default;
	Vector(const Vector&) = default;

	/// Constructor for dynamically-size vectors.  This can also be
	/// used for statically sized vectors in which case the argument
	/// is ignored.  The values of the vector are uninitialised
	explicit inline Vector(int size_in) : Base::template VLayout<Size, Precision>(size_in) {}

	/// Constructor used when constructing a vector which references
	/// other data, e.g.
	/// @code 
	/// double[] d = {1,2,3};
	/// Vector<3,double,Reference> v(d);
	/// @endcode
	explicit inline Vector(PointerType data) : Base::template VLayout<Size, Precision> (data) {}

	/// Construct from an initializer list
	/// @code 
	/// Vector<Resizable> v = {1, 2, 3};
	/// @endcode
	template<typename X=Precision>
	inline  Vector(std::initializer_list<Precision> init, typename std::enable_if<sizeof(X)&&(Size <0)>::type* = 0)
	:Base::template VLayout<Size, Precision> (init) 
	{}


	/// Construct from an array by reference
	/// @code 
	/// Vector<3> v = {1, 2, 3};
	/// @endcode
	template<typename Precision2, int Size2>
	inline  Vector(const Precision2(&init)[Size2])
	:Base::template VLayout<Size, Precision> (init) 
	{}

	/// Constructor used when constructing a dynamic vector which references
	/// other data, e.g.
	/// @code 
	/// double[] d = {1,2,3};
	/// Vector<Dynamic,double,Reference> v(d,3);
	/// @endcode
	inline Vector(PointerType data, int size_in) : Base::template VLayout<Size, Precision> (data, size_in) {}

	/// internal constructor
	inline Vector(PointerType data_in, int size_in, int stride_in, Internal::Slicing)
  : Base::template VLayout<Size, Precision>(data_in, size_in, stride_in) {}
	
	using Base::template VLayout<Size, Precision>::size;
	using Base::template VLayout<Size, Precision>::try_destructive_resize;

	/// construction from Operator object
	///
	/// This is used to implement return value optimisation for
	/// vectors created from the product of a matrix and a vector, or
	/// another object like Ones
	template <class Op>
	inline Vector(const Operator<Op>& op)
		: Base::template VLayout<Size, Precision> (op)
	{
		op.eval(*this);
	}

	// Copy construction is a very special case. Copy construction goes all the
	// way down to the bottom. GenericVBase has no idea how to copy itself.
	// However, the underlying allocator objects do.  In the case of static sized
	// objects, C++ automatically copies the data.  For slice objects, C++ copies
	// all parts (pointer and size), which is correct.  For dynamically sized
	// non-slice objects the copying has to be done by hand.
	
	// inline Vector(const Vector&from);

	/// constructor from arbitrary vector
	template<int Size2, typename Precision2, typename Base2>
	inline Vector(const Vector<Size2,Precision2,Base2>& from):
		Base::template VLayout<Size, Precision>(from.size()) {
		operator=(from);
	}

	/// @}

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

	/// @name Accessing elements
	/// @{

	/// access an element of the vector
	///
	/// can be used as an l-value ie
	/// @code
	/// Vector<3> v;
	/// v[0] = 10;
	/// @endcode
	///
	/// @internal
	/// This method is not defined by Vector: it is inherited.
	Precision& operator[] (int i);

	/**
		@overload
	*/
	const Precision& operator[] (int i) const;

	/// @}

#endif

	
	/// @name Assignment
	/// @{

	/// operator = from copy
	/// A size mismatch is a fatal error, unless the destination
	/// is resizable.
	inline Vector& operator= (const Vector& from){
		try_destructive_resize(from.size());
		SizeMismatch<Size,Size>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}

	/// operator = another Vector
	/// A size mismatch is a fatal error, unless the destination
	/// is resizable.
	template<int Size2, typename Precision2, typename Base2>
	Vector<Size,Precision,Base >& operator= (const Vector<Size2, Precision2, Base2>& from){
		try_destructive_resize(from.size());
		SizeMismatch<Size,Size2>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}

	/// assignment from an Operator object
	/// Assignment from sized operators causes a resize
	/// of Resizable Vectors. Assignment from unsized
	/// operators dows not.
	template <class Op>
	inline Vector & operator=(const Operator<Op>& op){
		try_destructive_resize(op);
		op.eval(*this);
		return *this;
	}
	/// @}

	/// @name Operators on the vector
	/// @{

	/// divide this vector by a constant
	Vector& operator/=(const Precision rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]/=rhs;
		return *this;
	}
	
	/// multiply this vector by a constant
	Vector& operator*=(const Precision rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]*=rhs;
		return *this;
	}
	
	/// add another vector onto this one
	template<int Size2, class Precision2, class Base2>
	Vector& operator+=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]+=rhs[i];
		return *this;
	}
	
	/// add an Operator object onto this vector
	///
	/// this is used to handle cases such as:
	/// @code
	/// Vector<3> v;
	/// v+=Ones
	/// @endcode
	template<class Op>
	Vector& operator+=(const Operator<Op>& op)
	{
		op.plusequals(*this);
		return *this;
	}		

	template<class Op>
	Vector& operator-=(const Operator<Op>& op)
	{
		op.minusequals(*this);
		return *this;
	}		

	/// subtract another vector from this one
	template<int Size2, class Precision2, class Base2>
	Vector& operator-=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]-=rhs[i];
		return *this;
	}

	/// @}

	/// @name Comparison
	/// @{

	/// Test for equality with another vector
	template<int Size2, class Precision2, class Base2>
	bool operator==(const Vector<Size2, Precision2, Base2>& rhs) const {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
		  if((*this)[i]!=rhs[i])
		    return 0;
		return 1;
	}

	/// Test for inequality with another vector
	template<int Size2, class Precision2, class Base2>
	bool operator!=(const Vector<Size2, Precision2, Base2>& rhs) const {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
		  if((*this)[i]!=rhs[i])
		    return 1;
		return 0;
	}
	
	
	template<class Op>
	bool operator!=(const Operator<Op>& op)
	{
		return op.notequal(*this);
	}		
	
	/// @}

	/// @name Misc
	/// @{

	/// return me as a non const reference - useful for temporaries
	Vector& ref()
	{
		return *this;
	}

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

	/// What is the size of this vector?
	int size() const;
	
	/// Resize the vector. This is only provided if the vector is
	/// declared as Resizable. Existing elements are retained, new
	/// elements are uninitialized. Resizing has the same efficiency
	/// guarantees as <code>std::vector</code>.
	/// @param s The new size.
	///
	/// @internal
	/// This method is not defined by Vector: it is inherited.
	void resize(int s);

	/// Return a pointer to the first element of the vector. This method
	/// is only provided for non-slice vectors, i.e. a subset of the cases
	/// where the memory is guaranteed to be contiguous.
	///
	/// @internal
	/// This method is not defined by Vector: it is inherited.
	Precision* get_data_ptr();



	/// @}

	/// @name Reshaping, sub-vectors and matrices
	//@{
	/**
	   Convert this vector into a 1-by-Size matrix, i.e. a matrix which has this
	   vector as its only row.
	   @code
	   Vector<3> a = makeVector(1,2,3);
	   Matrix<1,3> m = a.as_row();  // now m = [1 2 3]
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	Matrix<1, Size, Precision> as_row();
  
	/**
	   Convert this vector into a Size-by-1 matrix, i.e. a matrix which has this
	   vector as its only column.
	   @code
	   Vector<3> a = makeVector(1,2,3);
	   Matrix<3,1> m = a.as_col();   // now m = [1 2 3]'
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	Matrix<Size, 1, Precision> as_col();
  
	/**
	   Convert this vector into a Diagonal Size-by-Size matrix, i.e. a matrix which is
	   zero everywhere except on the diagonal and the diagonal contains the values from this vector
	   @code
	   Vector<3> v = makeVector(1,2,3);
	   Vector<3> v2 = makeVector(2,3,4);
	   Vector<3> v3 = v.as_diagonal() * v2; // now v3 = (2,6,12)
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	DiagonalMatrix<Size,Precision> as_diagonal();

	/**
	   Extract a sub-vector. The vector extracted will be begin at element Start
	   and will contain the next Length elements.
	   @code
	   Vector<5> a = makeVector(1,2,3,4,5);
	   Extract the three elements starting from element 2
	   Vector<3> b = a.slice<2,3>();  /// b = [3 4 5]
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	template<Start, Length>
	const Vector<Length,Precision>& slice() const;
  
	/**
	   Extract a sub-vector. The vector extracted will be begin at element Start
	   and will contain the next Length elements. This version can be used as an
	   l-value as well as an r-value
	   @code
	   Vector<5> a = makeVector(1,2,3,4,5);
	   Vector<2> b = makeVector(8,9);
	   // replace the two elements starting from element 1 with b
	   a.slice<1, 2>() = b;       /// now a = [1 8 9 4 5]
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	template<Start, Length>
	Vector<Length,Precision>& slice();
  
	/**
	   Extract a sub-vector with runtime parameters.
	   The vector extracted will be begin at element start and will contain the next
	   length elements.
	   @code
	   Vector<5> a = makeVector(1,2,3,4,5);
	   Extract the three elements starting from element 2
	   Vector<> b = a.slice(2,3);  /// b = [3 4 5]
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	template<Start, Length>
	const Vector<Length,Precision>& slice() const;
  
	/**
	   Extract a sub-vector with runtime parameters, which can be used as an
	   l-value.
	   The vector extracted will be begin at element start and will contain the next
	   length elements.
	   @code
	   Vector<5> a = makeVector(1,2,3,4,5);
	   Extract the three elements starting from element 2
	   a.slice(2,3)[0] = 17;  /// a -> [1 2 17 4 5]
	   @endcode
	   @internal
	   This method is not defined by Vector: it is inherited.
	*/
	template<Start, Length>
	Vector<Length,Precision>& slice();
	//@}

#endif

};

}
