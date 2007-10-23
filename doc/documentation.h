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
	- @link gDecomps Matrix decompositions@endlink (@link TooN::LU LU@endlink, @link TooN::SVD SVD@endlink, @link TooN::SymEigen symmetric eigen decomposition @endlink)
	- @link gOptimize Function optimization@endlink (@link TooN::DownhillSimplex Downhill Simplex@endlink) 
	- Some particular types of @link gTransforms transformation matrices@endlink (@link TooN::SO3 SO3@endlink and @link TooN::SE3 SE3@endlink) 
	- Solvers for systems of @link gEquations linear equations@endlink using @link TooN::WLS weighted@endlink or @link TooN::IRLS iteratively-reweighted@endlink least squares.
	
 It provides classes for statically- (known at compile time) and dynamically- (unknown at compile time) sized vectors and matrices and it delegates advanced functions (like SVD or multiplication of large matrices) to LAPACK and BLAS (this means you will need libblas and liblapack).

The library makes substantial internal use of templates to achieve run-time speed efficiency whilst retaining a (reasonably) clear programming syntax. Because of this use of templates, you will need a recent compiler (for example it does not work with g++ before version 3.0). Programs that use it can also take a \e long time to compile with -O3 under g++. Currently the library only supports double precision vectors and matrices. 

Why use this library?
 - because it supports statically sized vectors and matrices very efficiently (without introducing metadata)*
 - because it provides type safety for statically sized vectors and matrices (you can't attempt to multiply a 3x4 matrix and a 2-vector).
 - because it supports transposition and subscripting of matrices (to obtain a vector) very efficiently
 - because it exploits LAPACK and BLAS (for which optimised versions exist on many platforms).

(*) This is not always true (see Implementation, below)
 
\section sGetting How to get

The library can be obtained from savannah using anonymous cvs with the command:

<code> cvs -z3 -d:pserver:anonymous\@cvs.savannah.nongnu.org:/sources/toon co %TooN </code>

\section sUsage How to use

\subsection ssCompiler Compiler setup

	- Make sure you have a suitable compiler (g++ version < 3 is no good)
	- Make sure the directory containing TooN/ is somewhere on your search path
	- Add <code>#include <TooN/TooN.h></code>, and any other header files you might need e.g. <code>TooN/helpers.h</code> (for some more matrix and vector functions) or <code>TooN/SVD.h</code> (for the singular value decomposition) to your source code. 
	- Add <code>using namespace %TooN</code> to your code (or prefix class declarations with <code>%TooN::</code>).
	- You will need to link with <code>-llapack -lblas</code> (and <code>-lg2c</code> for g++. I'm not sure about other compilers). This means you will also need liblapack.{a,so} and libblas.{a,so} 

\subsection ssExamples Examples

Create two vectors and work out their inner (dot), outer and cross products
@code
// Initialise the vectors
Vector<3> a;
a = 3,5,0;
Vector<3> b;
b = 4,1,3;
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
Vector<3> v;
v = 1,2,3;
// Initialise a matrix
double d[2][3] = {{2,4,5},{6,8,9}};
Matrix<2,3> M(d);
// Now perform calculations
Vector<2> v2 = M*v;  // OK - answer is a static 2D vector
Vector<> v3 = M*v;   // OK - vector is determined to be 2D at runtime
Vector<> v4 = v*M;   // Compile error - dimensions of matrix and vector incompatible
@endcode

See the detailed documentation for @link TooN::Vector Vector@endlink, @link TooN::Matrix Matrix@endlink and
the various @link gDecomps matrix decompositions @endlink for more examples.



\section sImplementation Implementation

\subsection ssStatic Static-sized vectors and matrices

One aspect that makes this library efficient is that when you declare a 3-vector, all you get are 3 doubles - there's no metadata. So <code>sizeof(Vector<3>)</code> is 24. This means that when you write <code>Vector<3> v;</code> the data for <code>v</code> is allocated on the stack and hence <code>new</code>/<code>delete</code> (<code>malloc</code>/<code>free</code>) overhead is avoided. However, for large vectors and matrices, this would be a Bad Thing since <code>Vector<1000000> v;</code> would result in an object of 8 megabytes being allocated on the stack. I don't know about you, but my whole stack is only that big. %TooN gets around that problem by having a cutoff at which statically sized vectors are allocated on the heap. This is completely transparent to the programmer, the objects' behaviour is unchanged and you still get the type safety offered by statically sized vectors and matrices. The cutoff size at which the library changes the representation is defined in <code>toon.h</code> as the <code>const int MaxStackSize</code> in the class <code>NUMERICS</code>.

When you apply the subscript operator to a <code>Matrix<3,3></code> and the function simply returns the apropriate hunk of memory as a vector \e reference (i.e. it basically does no work). This avoids copying and also allows the resulting vector to be used as an l-value. Similarly the transpose operation applied to a matrix returns the memory corresponding to the matrix as a reference to a matrix with the opposite layout which also means the transpose can be used as an l-value so <code>M1 = M2.T();</code> and <code>M1.T() = M2;</code> do exactly the same thing.

\subsection ssDynamic Dynamic sized vectors and matrices

These are implemented in the obvious way using metadata with the rule that the object that allocated on the heap also deallocates. Other objects may reference the data (e.g. when you subscript a matrix and get a vector).

\subsection ssLazy Return value optimisation vs Lazy evaluation
When you write <code>v1 = M * v2;</code> a naive implementation will compute <code>M * v2</code> and store the result in a temporary object. It will then copy this temporary object into <code>v1</code>. A method often advanced to avoid this is to have <code>M * v2</code> simply return an special object <code>O</code> which contains references to <code>M</code> and <code>v2</code>. When the compiler then resolves <code>v1 = O</code>, the special object computes <code>M*v2</code> directly into <code>v1</code>. This approach is often called lazy evaluation and the special objects lazy vectors or lazy matrices. Stroustrup (The C++ programming language Chapter 22) refers to them as composition closure objects or compositors.

The killer is this: <b>What if v1 is just another name for v2?</b> i.e. you write something like <code>v = M * v;</code>. In this case the semantics have been broken because the values of <code>v</code> are being overwritten as the computation progresses and then the remainder of the computation is using the new values. In this library <code>v1</code> in the expression could equally well alias part of <code>M</code>, thus you can't even solve the problem by having a clever check for aliasing between <code>v1</code> and <code>v2</code>. This aliasing problem means that the only time the compiler can assume it's safe to omit the temporary is when <code>v1</code> is being constructed (and thus cannot alias anything else) i.e. <code>Vector<3> v1 = M * v2;</code>.

%TooN provides this optimisation by providing the compiler with the opportunity to use a return value optimisation. It does this by making <code>M * v2</code> call a special constructor for <code>Vector<3></code> with <code>M</code> and <code>v2</code> as arguments. Since nothing is happening between the construction of the temporary and the copy construction of <code>v1</code> from the temporary (which is then destroyed), the compiler is permitted to optimise the construction of the return value directly into <code>v1</code>.

Because a naive implemenation of this strategy would result in the vector and matrix classes having a very large number of constructors, these classes are provided with template constructors that take a standard form. The code that does this, declared in the header of class <code>Vector</code> is:
@code
// constructor from 2-ary operator
template <class LHS, class RHS, class Op>
inline Vector(const LHS& lhs, const RHS& rhs, const Operator<Op>&){Op::eval(*this,lhs,rhs);}
@endcode
The third argument of the constructor is a dummy, used to specify the construction method because you the standard doesn't allow you to supply template arguments when you call a constructor. Since the argument is unused, my compiler omits it (and I hope yours does too). 

\subsection ssHow How it all really works

This documentation is generated from a cleaned-up version of the interface, hiding the implementation 
that allows all of the magic to work. If you want to know more and can understand idioms like:
@code
template <int Size, class AllocZone>
class FixedVAccessor : public AllocZone {
   ...
};
@endcode
and
@code
template <int Size>
class Vector : public FixedVector<Size, FixedVAccessor<Size,typename SizeTraits<Size>::get_zone> > {
   ...
};
@endcode
Then take a look at the source code ... 
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

/// @defgroup gOptimize Function optimization
/// Classes to perform function optimization.
