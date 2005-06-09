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
// A proxy version of the SE3 class,
// cleaned up to present a comprehensible
// version of the SE3 interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

/// All classes and functions are within this namespace
namespace TooN 
{


/// @class SE3 se3doc.h TooN/SE3.h
/// Represent a three-dimensional Euclidean transformation (a rotation and a translation). 
/// This can be represented by a \f$3\times\f$4 matrix operating on a homogeneous co-ordinate, 
/// so that a vector \f$\underline{x}\f$ is transformed to a new location \f$\underline{x}'\f$
/// by
/// \f[\begin{aligned}\underline{x}' &= E\times\underline{x}\\ \begin{bmatrix}x'\\y'\\z'\end{bmatrix} &= \begin{pmatrix}r_{11} & r_{12} & r_{13} & t_1\\r_{21} & r_{22} & r_{23} & t_2\\r_{31} & r_{32} & r_{33} & t_3\end{pmatrix}\begin{bmatrix}x\\y\\z\\1\end{bmatrix}\end{aligned}\f]
/// 
/// This transformation is a member of the Special Euclidean Lie group SE3. These can be parameterised
/// six numbers (in the space of the Lie Algebra). In this class, the first three parameters are a
/// translation vector while the second three are a rotation vector, whose direction is the axis of rotation
/// and length the amount of rotation (in radians), as for SO3
/// @ingroup gTransforms
class SE3 {
  friend SE3 operator*(const SO3& lhs, const SE3& rhs);
  friend std::istream& operator>> (std::istream& is, SE3& rhs);
public:
	/// Default constructor. Initialises the the rotation to zero (the identity) and the translation to zero
  SE3();

  /// Returns the rotation part of the transformation as a SO3
  SO3& get_rotation(){return my_rotation;}
  /// @overload
  const SO3& get_rotation() const {return my_rotation;}
  /// Returns the translation part of the transformation as a Vector
  Vector<3>& get_translation() {return my_translation;}
  /// @overload
  const Vector<3>& get_translation() const {return my_translation;}

  /// Exponentiate a Vector in the Lie Algebra to generate a new SE3.
  /// See the Detailed Description for details of this vector.
  /// @param vect The Vector to exponentiate
    static SE3 exp(const Vector<6>& vect);
  /// Take the logarithm of the matrix, generating the corresponding vector in the Lie Algebra.
  /// See the Detailed Description for details of this vector.
  Vector<6> ln() const;

  /// Returns a SE3 representing the inverse transformation.
  /// An SE3 is \f$ \left[R|t \right] \f$ (where \f$R\f$ is a rotation matrix and
  /// \f$t\f$ the translation), so the inverse is \f$\left[R^T|-R^Tt\right]\f$
  SE3 inverse() const;

  /// Right-multiply by another SE3 (concatenate the two transformations)
  /// @param rhs The multipier
  SE3& operator *=(const SE3& rhs);
  /// Right-multiply by another SE3 (concatenate the two transformations)
  /// @param rhs The multipier
  SE3 operator *(const SE3& rhs) const;

  /// Returns the i-th generator multiplied by a vector. 
  /// The generators of a Lie group are the basis for the space of the Lie algebra.
  /// For %SE3, the generators are six \f$4\times4\f$ matrices representing
  /// the six possible (linearised) degrees of freedom. These matrices are usally sparse,
  /// and are usually obtained to immediately multiply them by a vector, so
  /// this function provides a fast way of doing this operation.
  /// @param i The required generator
  /// @param pos The vector to multuiply by the generator
  static Vector<4> generator_field(int i, Vector<4> pos);

  /// Transfer a vector in the Lie Algebra from one
  /// co-ordinate frame to another. This is the operation such that for a matrix 
  /// \f$ B \f$, 
  /// \f$ e^{\text{Adj}(v)} = Be^{v}B^{-1} \f$
  /// @param v The Vector to transfer
  inline void adjoint(Vector<6>& v)const;

  /// Transfer a matrix in the Lie Algebra from one
  /// co-ordinate frame to another. This is the operation such that for a matrix 
  /// \f$ B \f$, 
  /// \f$ e^{\text{Adj}(v)} = Be^{v}B^{-1} \f$
  /// @param M The Matrix to transfer
  inline void adjoint(Matrix<6,6>& M)const;

  /// Transfer covectors between frames (using the transpose of the inverse of the adjoint)
  /// so that trinvadjoint(vect1) * adjoint(vect2) = vect1 * vect2
  inline void trinvadjoint(Vector<6>& v)const;

  ///@overload
  inline void trinvadjoint(Matrix<6,6>& M)const;
};


/// Left-multiply an SE3 by an SO3 
/// @relates SE3
SE3 operator*(const SO3& lhs, const SE3& rhs);

/// Write an SE3 to a stream 
/// @relates SE3
std::ostream& operator <<(std::ostream& os, const SE3& rhs);

/// Reads an SE3 from a stream 
/// @relates SE3
std::istream& operator>>(std::istream& is, SE3& rhs);

/// Right-multiply by a Vector
/// @relates SE3
Vector<4> operator*(const SE3& lhs, const Vector<4>& rhs);


/// Left-multiply by a Vector
/// @relates SE3
Vector<4> operator*(const Vector<4>& lhs, const SE3& rhs);

/// Right-multiply by a Matrix
/// @relates SE3
template <int RHS>
Matrix<4,RHS> operator*(const SE3& lhs, const Matrix<4,RHS>& rhs);

/// Left-multiply by a Matrix
/// @relates SE3
template <int LHS>
Matrix<LHS,4> operator*(const Matrix<LHS,4>& lhs, const SE3& rhs);

}

#endif
