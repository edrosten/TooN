/*
    Copyright (c) 2005 Paul Smith, 2009 Edward Rosten

	Permission is granted to copy, distribute and/or modify this document under
	the terms of the GNU Free Documentation License, Version 1.2 or any later
	version published by the Free Software Foundation; with no Invariant
	Sections, no Front-Cover Texts, and no Back-Cover Texts.

    You should have received a copy of the GNU Free Documentation License
    License along with this library; if not, write to the Free Software
    Foundation, Inc.
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

*/
///////////////////////////////////////////////////////
// General Doxygen documentation
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// The main title page
/**
@mainpage

\section sIntro Introduction

The %TooN library is a set of C++ header files which provide basic numerics facilities:
	- @link TooN::Vector Vectors@endlink and @link TooN::Matrix matrices@endlink
	- @link gDecomps Matrix decompositions@endlink
	- @link gOptimize Function optimization@endlink
	- @link gTransforms Parameterized matrices (eg transformations)@endlink 
	- @link gEquations linear equations@endlink

It provides classes for statically- (known at compile time) and dynamically-
(unknown at compile time) sized vectors and matrices and it delegates advanced
functions (like SVD or multiplication of large matrices) to LAPACK and BLAS
(this means you will need libblas and liblapack).

The library makes substantial internal use of templates to achieve run-time
speed efficiency whilst retaining a clear programming syntax.

Why use this library?
 - Because it supports statically sized vectors and matrices very efficiently.
 - Because it provides extensive type safety for statically sized vectors and matrices (you can't attempt to multiply a 3x4 matrix and a 2-vector).
 - Because it supports transposition, subscripting and slicing of matrices (to obtain a vector) very efficiently.
 - Because it interfaces well to other libraries.
 - Because it exploits LAPACK and BLAS (for which optimised versions exist on many platforms).

\section sUsage How to use TooN
This section is arranged as a FAQ. Most answers include code fragments. Assume
<code>using namespace TooN;</code>.

 - \ref sDownload
 - \ref sStart
 - \ref sCreateVector
 - \ref sCreateMatrix
 - \ref sFunctionVector
 - \ref sGenericCode
 - \ref sElemOps
 - \ref sInitialize
 - \ref sScalars
 - \ref ssExamples
 - \ref sNoResize
 - \ref sDebug
 - \ref sSlices
 - \ref sPrecision
 - \ref sSolveLinear
 - \ref sOtherStuff
 - \ref sHandyFuncs
 - \ref sNoInplace
 - \ref sColMajor
 - \ref sWrap
 - \ref sWrap "How do I interface to other libraries?"
 - \ref sImplementation

 	\subsection sDownload Getting the code and installing
	
	To get the code from cvs use:

	cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/toon co TooN

	The home page for the library with a version of this documentation is at:

	http://mi.eng.cam.ac.uk/~er258/cvd/toon.html

	The code will work as-is, and comes with a default configuration, which
	should work on any system.

	On a unix system, <code>./configure && make install </code> will  install
	TooN to the correct place.  Note there is no code to be compiled, but the
	configure script performs some basic checks.

	\subsection sStart Getting started

		To begin, just in include the right file:

		@code
		#include <TooN/TooN.h>
		@endcode

		Everything lives in the <code>TooN</code> namespace.

		Then, make sure the directory containing TooN is in your compiler's
		search path. If you use any decompositions, you will need to link
		against LAPACK, BLAS and any required support libraries. On a modern
		unix system, linking against LAPACK will do this automatically.


	\subsection sCreateVector How do I create a vector?

		Vectors can be statically sized or dynamically sized.

		@code
			Vector<3> v1;    //Create a static sized vector of size 3
			Vector<>  v2(4); //Create a dynamically sized vector of size 4
			Vector<Dynamic>  v2(4); //Create a dynamically sized vector of size 4
		@endcode

		See also \ref sPrecision.


	\subsection sCreateMatrix How do I create a matrix?

		Matrices can be statically sized or dynamically sized.

		@code
			Matrix<3> m;              //A 3x3 matrix (statically sized)
			Matrix<3,2>  m;           //A 3x2 matrix (statically sized)
			Matrix<>  m(5,6);         //A 5x6 matrix (dynamically sized)
			Matrix<3,Dynamic> m(3,6); //A 3x6 matrix with a dynamic number of columns and static number of rows.
			Matrix<Dynamic,2> m(3,2); //A 2x3 matrix with a dynamic number of rows and static number of columns.
		@endcode

		See also \ref sPrecision.


	\subsection sFunctionVector How do I write a function taking a vector?
	
		To write a function taking a local copy of a vector:
		@code
			template<int Size> void func(Vector<Size> v);
		@endcode

		To write a function taking any type of vector by reference:
		@code
		template<int Size, typename Base> void func(const Vector<Size, double, Base>& v);
		@endcode
		See also \ref sPrecision, \ref sGenericCode and \ref sNoInplace


	\subsection sElemOps What elementary operations are supported?
		
		Assignments are performed using <code>=</code>. See also 
		\ref sNoResize.

		These operators apply to vectors or matrices and scalars. 
		The operator is applied to every element with the scalar.
		@code
		*=, /=, *, / 
		@endcode
		
		Vector and vectors or matrices and matrices:
		@code
		+, -, +=, -= 
		@endcode
		
		Dot product:
		@code
		Vector * Vector
		@endcode

		Matrix multiply:
		@code
		Matrix * Matrix
		@endcode

		Matrix multiplying a column vector:
		@code
		Matrix * Vector
		@endcode

		Row vector multiplying a matrix:
		@code
		Vector * Matrix
		@endcode
		
		3x3 Vector cross product:
		@code
		Vector<3> ^ Vector<3> 
		@endcode

		All the functions listed below return slices. The slices 
		are simply references to the original data and can be used as lvalues.

		Getting the transpose of a matrix:
		@code
			Matrix.T()
		@endcode
		
		Accessing elements:
		@code
		Vector[i]     //get element i
		Matrix(i,j)   //get element i,j
		Matrix[i]     //get row i as a vector
		Matrix[i][j]  //get element i,j
		@endcode
		
		Turning vectors in to matrices:
		@code
		Vector.as_row() //vector as a 1xN matrix
		Vector.as_col() //vector as a Nx1 matrix
		@endcode

		Slicing:
		@code
		Vector.slice<Start, End>();                            //Static slice
		Vector.slice<>(start, end);                            //Dynamic slice
		Matrix.slice<RowStart, ColStart, NumRows, NumCols>();  //Static slice
		Matrix.slice<>(rowstart, colstart, numrows, numcols);  //Dynamic slice
		@endcode

		See also \ref sSlices

	\subsection sInitialize How I initialize a vector/matrix?

		Vectors and matrices start off uninitialized (filled with random garbage).
		They can be easily filled with zeros, or ones (see also TooN::Ones):
		@code
			Vector<3> v = Zeros;
			Matrix<3> m = Zeros
			Vector<>  v2 = Zeros(2); //Note in they dynamic case, the size must be specified
			Matrix<>  m2 = Zeros(2,2); //Note in they dynamic case, the size must be specified
		@endcode

		Vectors can be filled with makeVector:
		@code
			Vector<> v = makeVector(2,3,4,5,6);
		@endcode
		
		Matrices can be initialized to the identity matrix:
		@code
			Matrix<2> m = Idendity;
			Matrix<> m2 = Identity(3);
		@endcode
		note that you need to specify the size in the dynamic case.

		They can also be initialized with data from another source. See also \ref  sWrap.





	\subsection sScalars How do I add a scalar to every element of a vector/matrix? 
		
		Addition to every element is not an elementary operation in the same way
		as multiplication by a scalar. It is supported throught the ::Ones
		object:
		
		@code
			Vector<3> a, b;
			...
			b = a + Ones*3;       // b_i = a_i + 3
			a+= Ones * 3;         // a_i <- a_i + 3
		@endcode

		It is supported the same way on Matrix and slices.

	\subsection sNoResize Why does assigning mismatched dynamic vectors fail?
	
	Vectors are not generic containers, and dynamic vectors have been designed
	to have the same semantics as static vectors where possible. Therefore
	trying to assign a vector of length 2 to a vector of length 3 is an error,
	so it fails. See also \ref sResize

	\subsection sResize How do I resize a dynamic vector/matrix?

	You can't. Preventing resize allow all sorts of things to be const
	which is great for optimization. If you want a genuinely resizable
	structure, you may consider using a <code>std::vector</code>, and accessing
	it as a TooN object when appropriate. See \ref sWrap. Also, the
	speed and complexity of resizable matrices depends on the memory layout, so
	you may wish to use column major matrices as opposed to the default row
	major layout.

	\subsection sDebug What debugging options are there?

	By default, everything which is checked at compile time in the static case
	is checked at run-time in the dynamic case. In other words, slices and sizes
	are checked at run-time if need be. These checks can be disabled by defining
	the macros \c TOON_NDEBUG_SLICE and \c TOON_NDEBUG_SIZE respectively. Bounds are
	not checked by default. Bounds checking can be enabled by defining the macro
	\c TOON_CHECK_BOUNDS. None of these macros change the interface, so debugging
	code can be freely mixed with optimized code.

	Errors are manifested to a call to <code>std::abort()</code>.

	TooN does not initialize data in a Vector or Matrix.  For debugging purposes
	the following macros can be defined:
	- \c TOON_INITIALIZE_QNAN Sets every element of newly defined Vectors or
	  Matrixs to quiet NaN, if it exists, and 0 otherwise. Your code will not compile
	  if you have made a Vector or Matrix of a type which cannot be constructed
	  from a number.
	- \c TOON_INITIALIZE_SNAN Sets every element of newly defined Vectors or
	  Matrixs to signalling NaN, if it exists, and 0 otherwise. 
	- \c TOON_INITIALIZE_VAL Sets every element of newly defined Vectors or
	  Matrixs to the expansion of this macro.
	- \c TOON_INITIALIZE_RANDOM Fills up newly defined Vectors and Matrixs with
	  random bytes, to trigger non repeatable behaviour. The random number
	  generator is automatically seeded with a granularity of 1 second. Your
	  code will not compile if you have a Vector or Matrix of a non-POD type.

	\subsection sSlices What are slices?

	Slices are references to data belonging to another vector or matrix. Modifying
	the data in a slice modifies the original object. Likewise, if the original 
	object changes, the change will be reflected in the slice. Slices can be
	used as lvalues. For example:

	@code
		Matrix<3> m = Identity;

		m.slice<0,0,2,2>() *= 3; //Multiply the top-left 2x2 submatrix of m by 3.

		m[2] /=10; //Divide the third row of M by 10.

		m.T()[2] +=2; //Add 2 to every element of the second column of M.

		m[1].slice<1,2>() = makeVector(3,4); //Set m_1,1 to 3 and m_1,2 to 4
		
		m[0][0]=6;
	@endcode


	\subsection sPrecision Can I have a precision other than double?

	Yes!
	@code
		Vector<3, float> v;          //Static sized vector of floats
		Vector<Dynamic, float> v(4); //Dynamic sized vector of floats
	@endcode

	Likewise for matrix.

	\subsection sSolveLinear How do I invert a matrix / solve linear equations?
	
	You use the decomposition objects (see \ref sDecompos "below"), for example to solve Ax=b:

	@code
	Matrix<3> A;
	A[0]=makeVector(1,2,3);
	A[1]=makeVector(2,3,4);
	A[2]=makeVector(3,2,1);

	Vector<3> b = makeVector (2,4,5);

	// solve Ax=b using LU
	LU<3> luA(A);
	Vector<3> x1 = luA.backsub(b);

	// solve Ax=b using SVD
	SVD<3> svdA(A);
	Vector<3> x2 = svdA.backsub(b);
	@endcode
	
	Similarly for the other \ref sDecompos "decomposition objects"

	\subsection sDecompos  Which decomposisions are there?

	For general size matrices (not necessarily square) there are:
	@link TooN::LU LU @endlink, @link TooN::SVD SVD @endlink and gauss_jordan

	For square symmetric matrices there are:
	@link TooN::SymEigen SymEigen @endlink and @link TooN::Cholesky Cholesky @endlink

	If all you want to do is solve a single Ax=b then you may want gaussian_elimination

	\subsection sOtherStuff What other stuff is there:

		Optimization: WLS, IRLS, downhill_simplex, SO2, SE2, SO3, SE3


	\subsection sHandyFuncs What handy functions are there (normalize, identity, fill, etc...)?


	\subsection sNoInplace Why don't functions work in place?

	Consider the function:
	@code
		void func(Vector<3>& v);
	@endcode
	It can accept a <code>Vector<3></code> by reference, and operate on it 
	in place. A <code>Vector<3></code> is a type which allocates memory on the
	stack. A slice merely references memory, and is a subtly different type. To
	write a function taking any kind of vector (including slices) you can write:

	@code
		template<class Base> void func(Vector<3, double, Base>& v);
	@endcode

	A slice is a
	temporary object, and according to the rules of C++, you can't pass a
	temporary to a function as a non-const reference. TooN provides the
	<code>.ref()</code> method to escape from this restriction, by returning a
	reference as a non-temporary. You would then have to write:
	@code
		Vector<4> v;
		...
		func(v.slice<0,3>().ref());
	@endcode
	to get func to accept the slice.

	Alternatively, you can observe that only TooN objects with the default base
	class own the data. All other sorts are references, so copying them only
	copies the reference, and the referred data is the same. You could therefore
	write a function to forward on TooN objects with the default base:

	@code
		template<class Base> void func(Vector<3, double, Base> v); //This will operate in-place only on slices

		void func(Vector<3>& v) //This will catch any non-slices and forward them on.
		{
			func(v.as_slice());
		}
	@endcode

	However, please consider writing functions that do not modify structures in
	place. The \c unit function of TooN computes a unit vector given an input
	vector. In the following context, the code:
	@code
		//There is some Vector, which may be a slice, etc called v;
		v = unit(v);
	@endcode
	produces exactly the same compiler output as the hypothetical
	<code>Normalize(v)</code> which operates in place (for static vectors). Consult the ChangeLog 
	entries dated ``Wed 25 Mar, 2009 20:18:16'' and ``Wed  1 Apr, 2009 16:48:45''
	for further discussion.
	

	\subsection sColMajor Can I have a column major matrix?

	Yes!
	@code
		Matrix<3, 3, double, ColMajor> m;          //3x3 Column major matrix
	@endcode

	\subsection sWrap I have a pointer to a bunch of data. How do I turn it in to a vector/matrix without copying?
	To create a vector use:
	@code
	double d[]={1,2,3,4};
	Vector<4,double,Reference> v1(d);
	Vector<Dynamic,double,Reference> v2(d,4);
	@endcode
	Or, a functional form can be used:
	@code
	double d[]={1,2,3,4};

	wrapVector<4>(d);         //Returns a Vector<4>
	wrapVector<4,double>(d);  //Returns a Vector<4>
	
	wrapVector(d,3);          //Return a Vector<Dynamic> of size 3
	wrapVector<Double>(d,3);  //Return a Vector<Dynamic> of size 3
	@endcode

	To crate a matrix use
	@code
	double d[]={1,2,3,4,5,6};
	Matrix<2,3,double,Reference::RowMajor> m1(d);
	Matrix<2,3,double,Reference::ColMajor> m2(d);
	Matrix<Dynamic, Dynamic, double, Reference::RowMajor> m3(d, 2, 3);
	Matrix<Dynamic, 3, double, Reference::RowMajor> m4(d, 2, 3); // note two size arguments are required for semi-dynamic matrices
	@endcode

	\subsection sGenericCode How do I write generic code?
	
	The constructors for TooN objects are very permissive in that they 
	accept run-time size arguments for statically sized objects, and then 
	discard the values, This allows you to easily write generic code which 
	works for both static and dynamic inputs.

	Here is a function which mixes up a vector with a random matrix:
	@code
	template<int Size, class Precision, class Base> Vector<Size, Precision> mixup(const Vector<Size, Precision, Base>& v)
	{
		//Create a square matrix, of the same size as v. If v is of dynamic
		//size, then Size == Dynamic, and so Matrix will also be dynamic. In
		//this case, TooN will use the constructor arguments to select the
		//matrix size. If Size is a real size, then TooN will simply ighore
		//the constructor values.

		Matrix<Size, Size, Precision> m(v.size(), v.size());
		
		//Fill the matrix with random values that sum up to 1.
		Precision sum=0;
		for(int i=0; i < v.size(); i++)
			for(int j=0; j < v.size(); j++)
				sum += (m[i][j] = rand());
		
		m/= sum;

		return m * v;
	}
	@endcode

	Writing functions which safely accept multiple objects requires assertions
	on the sizes since they may be either static or dynamic. TooN's built in
	size check will fail at compile time if mismatched static sizes are given,
	and at run-time if mismatched dynamic sizes are given:
	
	@code
	template<int S1, class B1, int S2, class B2> void func_of_2_vectors(const Vector<S1, double, B1>& v1, const Vector<S2, double, B2>& v2)
	{
		//Ensure that vectors are the same size
		SizeMismatch<S1, S2>::test(v1.num_rows(), v2.num_rows());


	}
	@endcode


\subsection ssExamples Are there any examples?

Create two vectors and work out their inner (dot), outer and cross products
@code
// Initialise the vectors
Vector<3> a = makeVector(3,5,0);
Vector<3> b = makeVector(4,1,3);

// Now work out the products
double dot = a*b;                            // Dot product
Matrix<3,3> outer = a.as_col() * b.as_row(); // Outer product
Vector<3> cross = a ^ b;                     // Cross product

cout << "a:" << endl << a << endl;
cout << "b:" << endl << b << endl;
cout << "Outer:" << endl << outer << endl;
cout << "Cross:" << endl << cross << endl;
@endcode

Create a vector and a matrix and multiply the two together
@code
// Initialise a vector
Vector<3> v = makeVector(1,2,3);

// Initialise a matrix
Matrix<2,3> M(d);
M[0] = makeVector(2,4,5);
M[1] = makeVector(6,8,9);

// Now perform calculations
Vector<2> v2 = M*v;  // OK - answer is a static 2D vector
Vector<> v3 = M*v;   // OK - vector is determined to be 2D at runtime
Vector<> v4 = v*M;   // Compile error - dimensions of matrix and vector incompatible
@endcode


\subsection sImplementation How is it implemented

\subsubsection ssStatic Static-sized vectors and matrices

One aspect that makes this library efficient is that when you declare a
3-vector, all you get are 3 doubles - there's no metadata. So
<code>sizeof(Vector<3>)</code> is 24. This means that when you write
<code>Vector<3> v;</code> the data for <code>v</code> is allocated on the stack
and hence <code>new</code>/<code>delete</code>
(<code>malloc</code>/<code>free</code>) overhead is avoided. However, for large
vectors and matrices, this would be a Bad Thing since <code>Vector<1000000>
v;</code> would result in an object of 8 megabytes being allocated on the stack and
potentially overflowing it. %TooN gets around
that problem by having a cutoff at which statically sized vectors are allocated
on the heap. This is completely transparent to the programmer, the objects'
behaviour is unchanged and you still get the type safety offered by statically
sized vectors and matrices. The cutoff size at which the library changes the
representation is defined in <code>TooN.h</code> as the <code>const int
TooN::Internal::max_bytes_on_stack=1000;</code>.

When you apply the subscript operator to a <code>Matrix<3,3></code> and the
function simply returns a vector which points to the the apropriate hunk of memory as a reference
(i.e. it basically does no work apart from moving around a pointer). This avoids
copying and also allows the resulting vector to be used as an l-value. Similarly
the transpose operation applied to a matrix returns a matrix which referes to the 
same memory but with the opposite layout which also means
the transpose can be used as an l-value so <code>M1 = M2.T();</code> and
<code>M1.T() = M2;</code> do exactly the same thing.

<b> Warning: This also means that <code>M = M.T();</code> does the wrong thing.</b>
However, since .T() essentially costs nothing, it should be very rare that you need to do this.

\subsubsection ssDynamic Dynamic sized vectors and matrices

These are implemented in the obvious way using metadata with the rule that the
object that allocated on the heap also deallocates. Other objects may reference
the data (e.g. when you subscript a matrix and get a vector).

\subsection ssLazy Return value optimisation vs Lazy evaluation

When you write <code>v1 = M * v2;</code> a naive implementation will compute
<code>M * v2</code> and store the result in a temporary object. It will then
copy this temporary object into <code>v1</code>. A method often advanced to
avoid this is to have <code>M * v2</code> simply return an special object
<code>O</code> which contains references to <code>M</code> and <code>v2</code>.
When the compiler then resolves <code>v1 = O</code>, the special object computes
<code>M*v2</code> directly into <code>v1</code>. This approach is often called
lazy evaluation and the special objects lazy vectors or lazy matrices.
Stroustrup (The C++ programming language Chapter 22) refers to them as
composition closure objects or compositors.


The killer is this: <b>What if v1 is just another name for v2?</b> i.e. you
write something like <code>v = M * v;</code>. In this case the semantics have
been broken because the values of <code>v</code> are being overwritten as the
computation progresses and then the remainder of the computation is using the
new values. In this library <code>v1</code> in the expression could equally well
alias part of <code>M</code>, thus you can't even solve the problem by having a
clever check for aliasing between <code>v1</code> and <code>v2</code>. This
aliasing problem means that the only time the compiler can assume it's safe to
omit the temporary is when <code>v1</code> is being constructed (and thus cannot
alias anything else) i.e. <code>Vector<3> v1 = M * v2;</code>.

%TooN provides this optimisation by providing the compiler with the opportunity
to use a return value optimisation. It does this by making <code>M * v2</code>
call a special constructor for <code>Vector<3></code> with <code>M</code> and
<code>v2</code> as arguments. Since nothing is happening between the
construction of the temporary and the copy construction of <code>v1</code> from
the temporary (which is then destroyed), the compiler is permitted to optimise
the construction of the return value directly into <code>v1</code>.

Because a naive implemenation of this strategy would result in the vector and
matrix classes having a very large number of constructors, these classes are
provided with template constructors that take a standard form. The code that
does this, declared in the header of class <code>Vector</code> is: 

@code
	template <class Op>
	inline Vector(const Operator<Op>& op)
		: Base::template VLayout<Size, Precision> (op)
	{
		op.eval(*this);
	}
@endcode

\subsubsection ssHow How it all really works

This documentation is generated from a cleaned-up version of the interface, hiding the implementation 
that allows all of the magic to work. If you want to know more and can understand idioms like:
@code

template<int, typename, int, typename> struct GenericVBase;
template<int, typename> struct VectorAlloc;

struct VBase {
	template<int Size, class Precision>
	struct VLayout : public GenericVBase<Size, Precision, 1, VectorAlloc<Size, Precision> > {
	    ...
	};
};

template <int Size, class Precision, class Base=VBase>
class Vector: public Base::template VLayout<Size, Precision> {
   ...
};
@endcode

then take a look at the source code ... 
**/

///////////////////////////////////////////////////////
// Modules classifying classes and functions

/// @defgroup gLinAlg Linear Algebra
/// %Vector and matrix classes, and helpers.

/// @defgroup gDecomps Matrix decompositions
/// Classes to perform matrix decompositions, used to solve 
/// linear equations and provide information about matrices. These are wrappers for functionality
/// provided by the LAPACK library.

/// @defgroup gTransforms Transformation matrices
/// Classes to represent particular types of transformation matrix.

/// @defgroup gEquations Linear equation solvers
/// Classes to solve linear equations.

/** 
@defgroup gOptimize Function optimization

Classes and functions to perform function optimization.

@section gOneDim One dimensional function optimization

The following functions find the minimum of a 1-D function:
 - golden_section_search()
 - brent_line_search()

@section gMultiDim Multidimensional dimensional function optimization

The following classes perform multidimensional function minimization:
 - TooN::DownhillSimplex
 - TooN::ConjugateGradient

The mode of operation is to set up a mutable class, then repeatedly call an
iterate function. This allows different sub algorithms (such as termination
conditions) to be substituted in if need be.

@defgroup gTooN Main parts of TooN

@defgroup gInternal TooN internals
*/
