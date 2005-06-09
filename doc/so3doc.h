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
// A proxy version of the SO3 class,
// cleaned up to present a comprehensible
// version of the SO3 interface

#ifdef DOXYGEN_INCLUDE_ONLY_FOR_DOCS

/// All classes and functions are within this namespace
namespace TooN 
{

/// @class SO3 so3doc.h TooN/SO3.h
/// Class to represent a three-dimensional rotation matrix. Three-dimensional rotation
/// matrices are members of the Special Orthogonal Lie group SO3. This group can be parameterised
/// three numbers (a vector in the space of the Lie Algebra). In this class, the three parameters are the
/// finite rotation vector, i.e. a three-dimensional vector whose direction is the axis of rotation
/// and whose length is the angle of rotation in radians. Exponentiating this vector gives the matrix,
/// and the logarithm of the matrix gives this vector.
/// @ingroup gTransforms
class SO3 {
public:
  friend inline std::istream& operator>>(std::istream& is, SO3& rhs);
  friend inline std::istream& operator>>(std::istream& is, SE3& rhs);
  /// Default constructor. Initialises the matrix to the identity (no rotation)
  SO3();

  /// Assignment operator.
  SO3& operator=(const SO3& rhs);
  
  /// Assigment operator from a general matrix. This also calls coerce()
  /// to make sure that the matrix is a valid rotation matrix.
  SO3& operator=(const Matrix<3>& rhs);
  
  /// Modifies the matrix to make sure it is a valid rotation matrix.
  void coerce();

  /// Exponentiate a vector in the Lie algebra to generate a new SO3.
  /// See the Detailed Description for details of this vector.
  inline static SO3 exp(const Vector<3>& vect);
  
  /// Exponentiate a vector in the Lie algebra to generate a new SO3.
  /// See the Detailed Description for details of this vector.
  static SO3 exp(const double* vect);

  /// Take the logarithm of the matrix, generating the corresponding vector in the Lie Algebra.
  /// See the Detailed Description for details of this vector.
  Vector<3> ln() const;

  /// Individual member access. Access is row major i.e. <code>[4]</code> will give the 2,1 element.
  double operator[](int i);

  /// Returns the inverse of this matrix (=the transpose, so this is a fast operation)
  SO3 inverse() const;

  /// Right-multiply by another rotation matrix
  SO3& operator *=(const SO3& rhs);

  /// Right-multiply by another rotation matrix
  SO3 operator *(const SO3& rhs) const;

  /// Returns the SO3 as a Matrix<3>
  const Matrix<3>& get_matrix()const {return my_matrix;}

  /// Returns the i-th generator multiplied by a vector. 
  /// The generators of a Lie group are the basis for the space of the Lie algebra.
  /// For %SO3, the generators are three \f$3\times3\f$ matrices representing
  /// the three possible (linearised) rotations. These matrices are usally sparse,
  /// and are usually obtained to immediately multiply them by a vector, so
  /// this function provides a fast way of doing this operation.
  static Vector<3> generator_field(int i, Vector<3> pos);

  /// Transfer a vector in the Lie Algebra from one
  /// co-ordinate frame to another such that for a matrix 
  /// \f$ M \f$, the adjoint \f$Adj()\f$ obeys
  /// \f$ e^{\text{Adj}(v)} = Me^{v}M^{-1} \f$
  Vector<3> adjoint(Vector<3> v) const ;

 private:
  Matrix<3> my_matrix;
};

/// Write an SO3 to a stream 
/// @relates SO3
inline std::ostream& operator<< (std::ostream& os, const SO3& rhs);

/// Read from SO3 to a stream 
/// @relates SO3
inline std::istream& operator>>(std::istream& is, SO3& rhs);


/// Right-multiply by a Vector
/// @relates SO3
Vector<3> operator*(const SO3& lhs, const Vector<3>& rhs);

/// Left-multiply by a Vector
/// @relates SO3
Vector<3> operator*(const Vector<3>& lhs, const SO3& rhs);

/// Multiply two SO3 matrices
/// @relates SO3
template <int LHS>
Matrix<LHS,3> operator*(Matrix<LHS,3>& lhs, const SO3& rhs);


}

#endif
