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
// A proxy version of the Matrix class,
// cleaned up to present a comprehensible
// version of the Mector interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

#include <iostream>

/// All classes and functions are within this namespace
namespace TooN
{

/**
@class Matrix matrixdoc.h TooN/toon.h
A matrix.
Support is provided for all the usual matrix operations: 
- the (a,b) notation can be used to access an element directly
- the [] operator can be used to yield a vector from a matrix (which can be used
as an l-value)
- they can be added and subtracted
- they can be multiplied (on either side) or divided by a scalar on the right:
- they can be multiplied by matrices or vectors
- submatrices can be extracted using the templated slice() member function
- they can be transposed (and the transpose used as an l-value)
- inverse is \e not supported. Use one of the @link gDecomps matrix
decompositions @endlink instead

See individual member function documentation for examples of usage.

\par Statically-sized matrices

The library provides classes for statically and dynamically sized matrices. As
with @link Vector Vectors@endlink, statically sized matrices are more efficient,
since their size is determined at compile-time, not run-time.
To create a \f$3\times4\f$ matrix, use:
@code
Matrix<3,4> M;
@endcode
or replace 3 and 4 with the dimensions of your choice. If the matrix is square,
it can be declared as:
@code
Matrix<3> M;
@endcode
which just is a synonym for <code>Matrix<3,3></code>. Matrices can also be
constructed from pointers or static 1D or 2D arrays of doubles:
@code
void foo(double* dptr) {

  Matrix<3> M1 (dptr);

  double dvals1[9]={1,2,3,4,5,6};
  Matrix<2,3> M2 (dvals1);

  double dvals2[2][3]={{1,2,3},{4,5,6}};
  Matrix<2,3> M3 (dvals2);

  // ...
}
@endcode

\par Dynamically-sized matrices

To create a dynamically sized matrix, use:
@code
Matrix<> M(num_rows, num_cols);
@endcode
where \a num_rows and \a num_cols are integers which will be evaluated at run
time. Dynamically-sized matrices can be constructed from pointers in a similar
manner to static sized matrices:
@code
void bar(int r, int c, double* dptr) {
  Matrix<> M (r,c,dptr);
  // ...
}
@endcode

<code>Matrix<></code> is a synonym for <code> Matrix<Dynamic, Dynamic> </code> which is
<code>%Matrix<-1,-1></code>

\par Row-major and column-major

The library supports both row major (the default - but you can change this if
you prefer) and column major layout ordering. Row major implies that the matrix
is laid out in memory in raster scan order:
\f[\begin{matrix}\text{Row major} & \text {Column major}\\
\begin{bmatrix}1&2&3\\4&5&6\\7&8&9\end{bmatrix} &
\begin{bmatrix}1&4&7\\2&5&8\\3&6&9\end{bmatrix} \end{matrix}\f]
You can override the default for a specific matrix by specifying the layout when
you construct it:
@code
Matrix<3,3,double,ColMajor> M1;
Matrix<-1,-1,double,RowMajor> M2(nrows, ncols);
@endcode
In this case the precision template argument must be given as it precedes the layout argument

@ingroup gLinAlg
**/
template<int Rows, int Cols>
class Matrix
{
public:
  /// @name Constructors
  //@{
  
  /// Default constructor.
  /// For fixed sized matrices, this does nothing, i.e. does not guarantee to
  initialise the vector to any particular values.
  /// For dynamically sized matrices, this sets up a 0x0 matrix.
  Matrix();
  
  /**
  Construct a (fixed-size) matrix from a 1D array of doubles. By default the
  layout ordering is row major (along the rows first i.e. raster scan order),
  unless <code>ColMajor</code> is specified in the template arguments (see
  detailed description). No errors are generated if the array is the wrong
  length.
  @code
  double d[4] = {1, 2, 3, 4};
  Matrix<2> m(d);  // now m = [1 2]
                    //         [3 4]
  @endcode
  */
  Matrix(double darray[Rows*Cols]);
  
  /**
  Construct a (fixed-size) matrix from a 2D array of doubles. No errors
  are generated if the array is the wrong length.
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);  // now m = [1 2 3]
                      //         [4 5 6]
  @endcode
  */
  Matrix(double darray[Rows][Cols]);
  
  /// Copy constructor.
  Matrix(const Matrix<Rows, Cols>& from);
  
  //@}
  
  
  /// @name Reading and writing elements
  //@{
  
  /// Assignment operator.
  /// For dynamically sized matrices this will not cause a resize if the sizes
are  mismatched.
  Matrix<Rows, Cols>& operator=(const Matrix<Rows, Cols>& from);
  
  /// Resize a dynamically sized matrix.
  resize(int rows, int cols);
  
  /**
  Access an element from the matrix.
  The index starts at zero, i.e. the top-left element is m(0, 0).
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  double e = m(1,2);     // now e = 6.0;
  @endcode
  */
  const double& operator() (int r, int c) const;
  
  /**
  Access an element from the matrix.
  This can be used as either an r-value or an l-value. The index starts at zero,
  i.e. the top-left element is m(0, 0).
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  m(1,2) = 8;     // now d = [1 2 3]
                  //         [4 5 8]
  @endcode
  */
  double& operator() (int r, int c);
  
  /**
  Access a row from the matrix.
  This can be used either as an r-value or an l-value. The index starts at zero,
  i.e. the first row is m[0]. To extract a column from a matrix, apply [] to the
  transpose of the matrix (see example). This can be used either as an r-value
  or an l-value. The index starts at zero, i.e. the first row (or column) is
  m[0].
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  Vector<3> v = m[1];       // now v = [4 5 6];
  Vector<2> v2 = m.T()[0];  // now v2 = [1 4];
  @endcode
  */
  const Vector& operator[] (int r) const;
  
  /**
  Access a row from the matrix.
  This can be used either as an r-value or an l-value. The index starts at zero,
  i.e. the first row is m[0]. To extract a column from a matrix, apply [] to the
  transpose of the matrix (see example). This can be used either as an r-value
  or an l-value. The index starts at zero, i.e. the first row (or column) is
  m[0].
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  Zero(m[0]);   // set the first row to zero
  Vector<2> v = 8,9;
  m.T()[1] = v; // now m = [0 8 0]
                //         [4 9 6]
  @endcode
  */
  Vector& operator[] (int r);
  
  /// How many rows does this matrix have?
  int num_rows() const;
  
  /// How many columns does this matrix have?
  int num_cols() const;
  
  /// What is the memory layout of this matrix?
  {RowMajor, ColMajor} layout const;
  //@}
  
  /// @name Transpose and sub-matrices
  //@{
  /**
  The transpose of the matrix. This is a very fast operation--it simply
  reinterprets a row-major matrix as column-major or vice-versa. This can be
  used as an l-value.
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  Zero(m[0]);   // set the first row to zero
  Vector<2> v = 8,9;
  m.T()[1] = v; // now m = [0 8 0]
                //         [4 9 6]
  @endcode
  */
  const Matrix<Cols, Rows>& T() const;
  
  /**
  The transpose of the matrix. This is a very fast operation--it simply
  reinterprets a row-major  matrix as column-major or vice-versa. The result can
  be used as an l-value.
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  Vector<2> v = 8,9;
  // Set the first column to v
  m.T()[0] = v; // now m = [8 2 3]
                //         [9 5 6]
  @endcode
  <b>This means that the semantics of <code>M=M.T()</code> are broken</b>. In
  general, it is not  necessary to say <code>M=M.T()</code>, since you can use
  M.T() for free whenever you need the transpose, but if you do need to, you
  have to use the Tranpose() function defined in <code>helpers.h</code>.
  */
  Matrix<Cols, Rows>& T();
  
  /**
  Extract a sub-matrix. The matrix extracted will be begin at element
  (Rstart, Cstart) and will contain the next Rsize by Csize elements.
  @code
  double d[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<3> m(d);
  Extract the top-left 2x2 matrix
  Matrix<2> b = m.slice<0,0,2,2>();  // b = [1 2]
                                      //     [4 5]
  @endcode
  */
  template<Rstart, Cstart, Rsize, Csize>
  const Matrix<Rsize, Csize>& slice() const;
  
  /**
  Extract a sub-matrix. The matrix extracted will be begin at element (Rstart,
  Cstart) and will contain the next Rsize by Csize elements. This can be used as
  either an r-value or an l-value.
  @code
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Matrix<2,3> m(d);
  Zero(m.slice<0,2,2,1>());  // b = [1 2 0]
                              //     [4 5 0]
  @endcode
  */
  template<Rstart, Cstart, Rsize, Csize>
  Matrix<Rsize, Csize>& slice();
  
  /**
  Extract a sub-matrix with runtime location and size. The matrix extracted will
  begin at element (rstart, cstart) and will
  contain the next rsize by csize elements.
  @code
  Matrix<> m(3,3);
  Extract the top-left 2x2 matrix
  Matrix<2> b = m.slice(0,0,2,2);
  @endcode
  */
  const Matrix<>& slice(int rstart, int cstart, int rsize, int csize) const;
  
  /**
  Extract a sub-matrix with runtime location and size, which can be used as
  an l-value. The matrix extracted will be begin at element (rstart, cstart) and
  will contain the next rsize by csize elements.
  @code
  Matrix<> m(3,3);
  Zero(m.slice(0,0,2,2));
  @endcode
  */
  Matrix<>& slice(int rstart, int cstart, int rsize, int csize);
  
  //@}
};

/// @name Input/output
//@{

/**
Write the matrix to a stream. The elements are space-separated with rows
being separated by new lines. The numbers are printed with a minimum field width
of 12 characters.
@relates Matrix
*/
template <int Rows, Cols>
std::ostream& operator<< (std::ostream& os, const Matrix<Rows, Cols>& v);

/**
Read a matrix from a stream. The numbers can be comma or white-space
separated, and are assumed to be in row-major layout (i.e. scanline order)
unless the matrix is defined as <code>ColMajor</code>
@relates Matrix
*/
template <int Rows, Cols>
std::istream& operator<< (std::istream& is, Matrix<Rows, Cols>& v);
//@}


}

#endif
