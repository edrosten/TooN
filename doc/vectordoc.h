/*
    Copyright (c) 2005 Paul Smith

	Permission is granted to copy, distribute and/or modify this document under
	the terms of the GNU Free Documentation License, Version 1.2 or any later
	version published by the Free Software Foundation; with no Invariant
	Sections, no Front-Cover Texts, and no Back-Cover Texts.

    You should have received a copy of the GNU Free Documentation License
    License along with this library; if not, write to the Free Software
    Foundation, Inc.
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

*/
// A proxy version of the vector class,
// cleaned up to present a comprehensible
// version of the Vector interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>

/// All classes and functions are within this namespace
namespace TooN
{

/// An object to create a Vector from a list of numbers. See the detailed documentation
/// for Vector for information on usage.
struct make_Vector
{
	/// Overloaded operator to allow this object to interpret a list of numbers after it
	/// as elements for a vector.
	VectorCreator operator,(const T& t);
};

/// Helper function to skip elements when using the <code>operator><<</code> initialisation
/// for Vector. See the detailed documentation for Vector for information on usage.
template <int N>
ComponentPlaceHolder<N> no_change();


/**
@class Vector vectordoc.h TooN/toon.h
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

The library provides classes for both statically- and dynamically-sized vectors. If you know what dimension of vector you're going to use (e.g. 3 to represent a point in 3D space), it's more efficient to statically sized vectors. The size of static vectors is determined at compile time; that of dynamically-sized vectors at run-time.

To create a 3-dimensional vector, use:
@code
Vector<3> v;
@endcode
and to create a vector of some other dimensionality just replace 3 with the positive integer of your choice, or some expression which the compiler can evaluate to an integer at compile time. Vectors can also be constructed from pointers or static arrays of doubles:
@code
void foo(double* dptr) {

  Vector<4> v1 (dptr);

  double dvals[4]={1,2,3,4};
  Vector<4> v2 (dvals);

  // ...
}
@endcode

The preferred way of defining a vector is to use make_Vector. The %make_Vector object constructs
a static vector initialised to the size and the contents of the comma-separated list following it. 
This magic is performed by overloading the comma operator, but for this to work, the expression
must be enclosed in brackets (see the examples). The %make_Vector vectors are real Vectors and 
so can be used anywhere where a vector is needed, not just in initialisations. For example
@code
// Create a vector initialised to [1 2 3];
Vector<3> v = (make_Vector, 1, 2, 3);
// Calculate the dot product with the vector [4 0 6]
double dot = v * (make_Vector, 4, 0, 6);
@endcode
Because the %make_Vector syntax creates actual vectors, compile-time checking is done to ensure
that all vectors defined in this way have the correct number of elements.

An alternative means of assigning values to vectors is to use the overloaded <code>operator=</code> function
with a comma-separated list of doubles (or anything that can be cast to a double):
@code
// Create a vector initialised to [1 2 3];
Vector<3> v;
v = 1, 2, 3;
@endcode
This is more concise than the %make_Vector syntax, but can only be used for assigning values to 
existing vectors, and only generates a compile-time error if the list is too long, not too 
short (in which case a run-time error is generated instead).

An third means of defining vectors is to use the overloaded <code>operator<<</code> function. This 
allows particular ranges of elements to be modified, since the operator sequentially replaces
elements in order, as the example below shows. The templated no_change() object allows runs of
elements to be skipped.
@code
// Create a vector initialised to [1 2 3 4];
Vector<4> v;
v << 1 << 2 << 3 << 4;
// Now modify just the first two elements
v << 5 << 6;   // v is now [5 6 3 4];
// Now modify the last two
v << no_change<2>() << 0 << 1;  // v is now [5 6 0 1]
// Now modify the last middle two
v << no_change() << 4 << 5;  // v is now [5 4 5 1]   (no_change() is a synonym for no_change<1>()
@endcode
Compile-time checking is done for attempts to over-fill a matrix. 

\par Dynamically-sized vectors

To create a dynamically sized vector, use:
@code
Vector<> v(size);
@endcode
where size is an integer which will be evaluated at run time. Dynamic sized vectors can be constructed from pointers in a similar manner to static sized vectors:
@code
void bar(int size, double* dptr) {
  Vector<> v1 (size,dptr);
  // ...
}
@endcode
Vector<> is actually a synonym for Vector<-1> being a template specialisation of Vector<N> with a special implementation.


\par Row vectors and column vectors

This library makes no distinction between row vectors and column vectors. Vectors that appear on the left of a multiplication are treated as row vectors while those that appear on the right are treated as column vectors (thus <code>v1*v2</code> means the dot product). This means that sometimes you have to be careful to include prarentheses since it is possible to write obscure stuff like <code>Vector<4> v4 = v1 * v2 * v3;</code>,
which in the absence of any extra parentheses means 'compute the dot product between <code>v1</code> and <code>v2</code> and then multiply <code>v3</code> by this scalar and assign to <code>v4</code>'.

If the row-column distinction is important, then vectors can be turned into matrices with one row or column by using as_row() or as_col():
@code
double d[3] = {1,2,3};
Vector<3> v(d);
Matrix<3,3> M = v.as_col() * v.as_row(); // creates a symmetric rank 1 matrix from v 
@endcode

@ingroup gLinAlg
**/
template<int Size>
class Vector
{
public:
  /// @name Constructors
  //@{
  
  /// Default constructor for vectors.
  /// For fixed-sized vectors, this does nothing, i.e. does not
  /// guarantee to initialise the vector to any particular values.
  /// For dynamically sized vectors, this sets the vector to have a length of 0.
  Vector();
  
  /// Constructor for dynamically-size vectors. Does nothing, i.e. does not
  /// guarantee to initialise the vector to any particular values
  Vector(int size);
  
  /**
  Construct a fixed-size vector from an array of doubles. No errors are
  generated if the array is the wrong length.
  @code
  double d[4] = {1, 2, 3, 4};
  Vector<4> a(d);  // now a = [1 2 3 4]
  @endcode
  */
  Vector(double darray[]);
  
  /**
  Construct a dynamically-sized vector from an array of doubles. No errors
  are generated if the array is the wrong length.
  @code
  double d[4] = {1, 2, 3, 4};
  Vector<> a(4, d);  // now a = [1 2 3 4]
  @endcode
  */
  Vector(int size, double* dptr);
  
  /// Copy constructor.
  Vector(const Vector<Size>& from);
  
  //@}
  
  
  /// @name Reading and writing elements, and defining vectors
  //@{
  
  /// Assignment operator. This will not cause a resize on dynamically sized
  vectors if the sizes are not matched.
  Vector<Size>& operator=(const Vector<Size>& from);
  
  
  /// Resize a dynamically sized vector.
  resize(int size);
  
  /**
  Create a (fixed-size) vector from a list of doubles.
  This operator should be followed by a comma-separated list of doubles, and the
  whole enclosed in brackets (see the example). The list of doubles most be of
  the correct length: if the wrong number of elements are provided, a compile
  error is generated.
  @code
  Vector<4> a;
  a = (make_Vector, 1, 2, 3, 4);    //  v = [ 1 2 3 4 ]
  a = make_Vector, 1, 2, 3, 4;      // compile error
  a = (make_Vector, 1, 2, 3, 4, 5); // compile error
  a = (make_Vector, 1, 2, 3);       // compile error;
  @endcode
  The <code>(make_Vector, 1, 2, 3, 4)</code> syntax generates a real vector
  object which can be used in other operations, not just for assigning to a
  Vector.
  */
  const Vector<Size>&  operator=(const make_Vector mv);
  
  /**
  Set the elements of a (fixed-size) vector from a list of doubles.
  This operator should be used together with the comma operator to provide a
  list of  doubles of the correct length. Providing too many numbers will
  generate a compile-time error, but if too few are provided, only a run-time
  error is generated. The operator=(make_Vector) syntax (see above) can check
  for both at compile time.
  @code
  Vector<4> a;
  a = 5,6,7,8;    // now a = [5 6 7 8]
  a = 5,6,7;      // run-time error (underfilling)
  a = 5,6,7,8,9;  // compile error (overfilling)
  @endcode
  */
  template <typename T>
  void operator=(const T& t);
  
  /**
  Set the elements of a (fixed-size) vector by sequential insertion.
  This operator replaces the elements of a vector starting from the beginning of
  the vector (see example). A compile error is generated if too many elements
  are provided. Single or multiple elements may be skipped and left unmodified
  by using the no_change() command.
  @code
  Vector<4> a;
  a << 9 << 10 << 11 << 12;       // now a = [9 10 11 12]
  a << 9 << 10 << 11 << 12 << 12; //compile error
  a << 13 << 14 << 15;            // now a = [13 14 15 12]
  a << 3 << no_change<2>() << 3;  // now a = [3 14 15 3];
  a << no_change() << 17 << no_change() << 17;
    // now a = [3 17 15 17], no_change() is no_change<1>()
  @endcode
  */
  template <typename T>
  void operator<<(const T& t);
  
  /**
  Access an element from the vector in the usual way.
  The index starts at zero, i.e. the first element is vect[0].
  @code
  Vector<3> a = 1,2,3;
  double d = a[1];     // now d = 2.0;
  @endcode
  */
  const double& operator[] (int i) const;
  
  /**
  Access an element from the vector in the usual way.
  This can be used as either an r-value or an l-value.
  The index starts at zero, i.e. the first element is vect[0].
  @code
  Vector<3> a = 1,2,3;
  double d = a[0]; // now d = 1.0;
  a[1] = 0;        // now a = [1 0 3];
  @endcode
  */
  double& operator[] (int i);
  
  /// Get the raw double array.
  const double* get_data_ptr() const;
  
  /// Get the raw double array.
  double* get_data_ptr();
  
  /// What is the size of this vector?
  int size() const;
  //@}
  
  /// @name Reshaping, sub-vectors and matrices
  //@{
  /**
  Convert this vector into a 1-by-Size matrix, i.e. a matrix which has this
  vector as its only row.
  @code
  Vector<3> a = 1,2,3;
  Matrix<1,3> m = a.as_row();  // now m = [1 2 3]
  @endcode
  */
  Matrix<1, Size> as_row();
  
  /**
  Convert this vector into a Size-by-1 matrix, i.e. a matrix which has this
  vector as its only column.
  @code
  Vector<3> a = 1,2,3;
  Matrix<3,1> m = a.as_col();   // now m = [1 2 3]'
  @endcode
  */
  Matrix<Size, 1> as_col();
  
  /**
  Extract a sub-vector. The vector extracted will be begin at element Start
  and will contain the next Length elements.
  @code
  Vector<5> a = 1,2,3,4,5;
  Extract the three elements starting from element 2
  Vector<3> b = a.slice<2,3>();  /// b = [3 4 5]
  @endcode
  */
  template<Start, Length>
  const Vector<Length>& slice() const;
  
  /**
  Extract a sub-vector. The vector extracted will be begin at element Start
  and will contain the next Length elements. This version can be used as an
  l-value as well as an r-value
  @code
  Vector<5> a = 1,2,3,4,5;
  Vector<2> b = 8,9;
  // replace the two elements starting from element 1 with b
  a.slice<1, 2>() = b;       /// now a = [1 8 9 4 5]
  @endcode
  */
  template<Start, Length>
  Vector<Length>& slice();
  
  /**
  Extract a sub-vector with runtime parameters.
  The vector extracted will be begin at element start and will contain the next
  length elements.
  @code
  Vector<5> a = (makeVector, 1,2,3,4,5);
  Extract the three elements starting from element 2
  Vector<3> b = a.slice(2,3);  /// b = [3 4 5]
  @endcode
  */
  template<Start, Length>
  const Vector<Length>& slice() const;
  
  /**
  Extract a sub-vector with runtime parameters, which can be used as an
  l-value.
  The vector extracted will be begin at element start and will contain the next
  length elements.
  @code
  Vector<5> a = (makeVector, 1,2,3,4,5);
  Extract the three elements starting from element 2
  a.slice(2,3)[0] = 17;  /// a -> [1 2 17 4 5]
  @endcode
  */
  template<Start, Length>
  const Vector<Length>& slice() const;
  //@}
};

/// @name Input/output
//@{

/// Write the elements of vector to a stream. The elements are space-separated and printed with a minimum
/// field width of 12 characters.
/// @relates Vector
template <int Size>
std::ostream& operator<< (std::ostream& os, const Vector<Size>& v);

/// Read a vector from a stream. The numbers can be comma or white-space separated.
/// @relates Vector
template <int Size>
std::istream& operator<< (std::istream& is, Vector<Size>& v);
//@}



}

#endif
