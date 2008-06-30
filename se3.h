
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
#ifndef __SE3_H
#define __SE3_H

#include <TooN/so3.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

class SE3 {
    friend SE3 operator*(const SO3& lhs, const SE3& rhs);
    friend std::istream& operator>> (std::istream& is, SE3& rhs);
  
public:
    inline SE3();
    template <class A> inline SE3(const SO3& R, const FixedVector<3,A>& T) : my_rotation(R), my_translation(T) {}
      

    inline SO3& get_rotation(){return my_rotation;}
    inline const SO3& get_rotation() const {return my_rotation;}
    inline Vector<3>& get_translation() {return my_translation;}
    inline const Vector<3>& get_translation() const {return my_translation;}

    static inline SE3 exp(const Vector<6>& vect);
    static inline Vector<6> ln(const SE3& se3);
    inline Vector<6> ln() const { return SE3::ln(*this); }

    inline SE3 inverse() const;

    inline SE3& operator *=(const SE3& rhs);
    inline SE3 operator *(const SE3& rhs) const { return SE3(my_rotation*rhs.my_rotation, my_translation + my_rotation*rhs.my_translation); }
    inline SE3& left_multiply_by(const SE3& left);

    static inline Vector<4> generator_field(int i, Vector<4> pos);


    template<class Accessor>
    inline void adjoint(FixedVector<6,Accessor>& vect)const;

    template<class Accessor>
    inline void trinvadjoint(FixedVector<6,Accessor>& vect)const;

    template <class Accessor>
    inline void adjoint(FixedMatrix<6,6,Accessor>& M)const;

    template <class Accessor>
    inline void trinvadjoint(FixedMatrix<6,6,Accessor>& M)const;

private:
    SO3 my_rotation;
    Vector<3> my_translation;
};

    template <class A> inline
    Vector<3> transform(const SE3& pose, const FixedVector<3,A>& x) { return pose.get_rotation()*x + pose.get_translation(); }
    
    template <class A> inline
    Vector<3> transform_inverse_depth(const SE3& pose, const FixedVector<3,A>& uvq)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3> DqT = R.template slice<0,0,3,2>() * uvq.template slice<0,2>() + R.T()[2] + uvq[2] * pose.get_translation();
	const double inv_z = 1.0/ DqT[2];
	return makeVector(DqT[0]*inv_z, DqT[1]*inv_z, uvq[2]*inv_z);
    }

    template <class A1, class A2> inline
    Vector<3> transform(const SE3& pose, const FixedVector<3,A1>& x, FixedMatrix<3,3,A2>& J_x)
    {
	J_x = pose.get_rotation().get_matrix();
	return pose * x;
    }


    template <class A1, class A2> inline
    Vector<3> inverse_depth_point_jacobian(const SE3& pose, const double q, const FixedVector<3,A1>& DqT, const double inv_z, FixedMatrix<3,3,A2>& J_uvq)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3>& T = pose.get_translation();
	const double u1 = DqT[0] * inv_z;
	J_uvq[0][0] = inv_z * (R[0][0] - u1*R[2][0]);
	J_uvq[0][1] = inv_z * (R[0][1] - u1*R[2][1]);
	J_uvq[0][2] = inv_z * (T[0] - u1 * T[2]);
	
	const double v1 = DqT[1] * inv_z;
	J_uvq[1][0] = inv_z * (R[1][0] - v1*R[2][0]);
	J_uvq[1][1] = inv_z * (R[1][1] - v1*R[2][1]);
	J_uvq[1][2] = inv_z * (T[1] - v1 * T[2]);
	
	const double q1 = q * inv_z;
	J_uvq[2][0] = inv_z * (-q1 * R[2][0]);
	J_uvq[2][1] = inv_z * (-q1 * R[2][1]);
	J_uvq[2][2] = inv_z * (1.0 - q1 * T[2]);

	return makeVector(u1, v1, q1);
    }
    

    template <class A1, class A2> inline
    Vector<3> transform_inverse_depth(const SE3& pose, const FixedVector<3,A1>& uvq, FixedMatrix<3,3,A2>& J_uvq)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3>& T = pose.get_translation();
	const double q = uvq[2];
	const Vector<3> DqT = R.template slice<0,0,3,2>() * uvq.template slice<0,2>() + R.T()[2] + q * T;
	const double inv_z = 1.0 / DqT[2];

	return inverse_depth_point_jacobian(pose, q, DqT, inv_z, J_uvq);
    }

    template <class A1, class A2, class A3> inline
    Vector<3> transform(const SE3& pose, const FixedVector<3,A1>& x, FixedMatrix<3,3,A2>& J_x, FixedMatrix<3,6,A3>& J_pose)
    {
	J_x = pose.get_rotation().get_matrix();
	Identity(J_pose.template slice<0,0,3,3>());
	const Vector<3> cx = pose * x;
	J_pose[0][3] = J_pose[1][4] = J_pose[2][5] = 0;
	J_pose[1][3] = -(J_pose[0][4] = cx[2]);
	J_pose[0][5] = -(J_pose[2][3] = cx[1]);
	J_pose[2][4] = -(J_pose[1][5] = cx[0]);
	return cx;
    }

    template <class A1, class A2, class A3> inline
    Vector<3> transform_inverse_depth(const SE3& pose, const FixedVector<3,A1>& uvq, FixedMatrix<3,3,A2>& J_uvq, FixedMatrix<3,6,A3>& J_pose)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3>& T = pose.get_translation();
	const double q = uvq[2];
	const Vector<3> DqT = R.template slice<0,0,3,2>() * uvq.template slice<0,2>() + R.T()[2] + q * T;
	const double inv_z = 1.0 / DqT[2];
	
	const Vector<3> uvq1 = inverse_depth_point_jacobian(pose, q, DqT, inv_z, J_uvq);
	const double q1 = uvq1[2];

	J_pose[0][1] = J_pose[1][0] = J_pose[2][0] = J_pose[2][1] = 0;
	J_pose[0][0] = J_pose[1][1] = q1;
	J_pose[0][2] = -uvq1[0] * q1;
	J_pose[1][2] = -uvq1[1] * q1;
	J_pose[2][2] = -q1 * q1;
	
	J_pose[0][3] = -uvq1[1]*uvq1[0];
	J_pose[0][4] = 1 + uvq1[0]*uvq1[0];
	J_pose[0][5] = -uvq1[1];
	
	J_pose[1][3] = -1 - uvq1[1]*uvq1[1];
	J_pose[1][4] = uvq1[0] * uvq1[1];
	J_pose[1][5] = uvq1[0];

	J_pose[2][3] = -uvq1[1]*q1;
	J_pose[2][4] = uvq1[0]*q1;
	J_pose[2][5] = 0;

	return uvq1;
    }

    template <class A1, class A2, class A3> inline
    Vector<2> project_transformed_point(const SE3& pose, const FixedVector<3,A1>& in_frame, FixedMatrix<2,3,A2>& J_x, FixedMatrix<2,6,A3>& J_pose)
    {
	const double z_inv = 1.0/in_frame[2];
	const double x_z_inv = in_frame[0]*z_inv;
	const double y_z_inv = in_frame[1]*z_inv;
	const double cross = x_z_inv * y_z_inv;
	J_pose[0][0] = J_pose[1][1] = z_inv;
	J_pose[0][1] = J_pose[1][0] = 0;
	J_pose[0][2] = -x_z_inv * z_inv;
	J_pose[1][2] = -y_z_inv * z_inv;
	J_pose[0][3] = -cross;
	J_pose[0][4] = 1 + x_z_inv*x_z_inv; 
	J_pose[0][5] = -y_z_inv;  
	J_pose[1][3] = -1 - y_z_inv*y_z_inv;
	J_pose[1][4] =  cross;
	J_pose[1][5] =  x_z_inv;    
	
	const TooN::Matrix<3>& R = pose.get_rotation().get_matrix();
	J_x[0][0] = z_inv*(R[0][0] - x_z_inv * R[2][0]);
	J_x[0][1] = z_inv*(R[0][1] - x_z_inv * R[2][1]);
	J_x[0][2] = z_inv*(R[0][2] - x_z_inv * R[2][2]);
	J_x[1][0] = z_inv*(R[1][0] - y_z_inv * R[2][0]);
	J_x[1][1] = z_inv*(R[1][1] - y_z_inv * R[2][1]);
	J_x[1][2] = z_inv*(R[1][2] - y_z_inv * R[2][2]);
	
	return makeVector(x_z_inv, y_z_inv);
    }


    template <class A1> inline
    Vector<2> transform_and_project(const SE3& pose, const FixedVector<3,A1>& x)
    {
	return project(transform(pose,x));
    }

    template <class A1> inline
    Vector<2> transform_and_project_inverse_depth(const SE3& pose, const FixedVector<3,A1>& uvq)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3>& T = pose.get_translation();
	const Vector<3> DqT = R.template slice<0,0,3,2>() * uvq.template slice<0,2>() + R.T()[2] + uvq[2] * T;
	return project(DqT);
    }

    template <class A1, class A2, class A3> inline
    Vector<2> transform_and_project(const SE3& pose, const FixedVector<3,A1>& x, FixedMatrix<2,3,A2>& J_x, FixedMatrix<2,6,A3>& J_pose)
    {
	return project_transformed_point(pose, transform(pose,x), J_x, J_pose);
    }

    template <class A1, class A2, class A3> inline
    Vector<2> transform_and_project_inverse_depth(const SE3& pose, const FixedVector<3,A1>& uvq, FixedMatrix<2,3,A2>& J_uvq, FixedMatrix<2,6,A3>& J_pose)
    {
	const Matrix<3>& R = pose.get_rotation().get_matrix();	
	const Vector<3>& T = pose.get_translation();
	const double q = uvq[2];
	const Vector<3> DqT = R.template slice<0,0,3,2>() * uvq.template slice<0,2>() + R.T()[2] + q * T;
	const Vector<2> uv = project_transformed_point(pose, DqT, J_uvq, J_pose);
	J_uvq.T()[2] = J_pose.template slice<0,0,2,3>() * pose.get_translation();
	J_pose.template slice<0,0,2,3>() *= q;
	return uv;
    }




// left multiply an SE3 by an SO3 
inline SE3 operator*(const SO3& lhs, const SE3& rhs);


// transfers a vector in the Lie algebra
// from one coord frame to another
// so that exp(adjoint(vect)) = (*this) * exp(vect) * (this->inverse())
template<class Accessor>
inline void SE3::adjoint(FixedVector<6,Accessor>& vect)const{
    vect.template slice<3,3>() = my_rotation * vect.template slice<3,3>();
    vect.template slice<0,3>() = my_rotation * vect.template slice<0,3>();
    vect.template slice<0,3>() += my_translation ^ vect.template slice<3,3>();
}

// tansfers covectors between frames
// (using the transpose of the inverse of the adjoint)
// so that trinvadjoint(vect1) * adjoint(vect2) = vect1 * vect2
template<class Accessor>
inline void SE3::trinvadjoint(FixedVector<6,Accessor>& vect)const{
    vect.template slice<3,3>() = my_rotation * vect.template slice<3,3>();
    vect.template slice<0,3>() = my_rotation * vect.template slice<0,3>();
    vect.template slice<3,3>() += my_translation ^ vect.template slice<0,3>();
}

template <class Accessor>
inline void SE3::adjoint(FixedMatrix<6,6,Accessor>& M)const{
    for(int i=0; i<6; i++){
	adjoint(M.T()[i]);
    }
    for(int i=0; i<6; i++){
	adjoint(M[i]);
    }
}
  
template <class Accessor>
inline void SE3::trinvadjoint(FixedMatrix<6,6,Accessor>& M)const{
    for(int i=0; i<6; i++){
	trinvadjoint(M.T()[i]);
    }
    for(int i=0; i<6; i++){
	trinvadjoint(M[i]);
    }
}


// operator ostream& <<
inline std::ostream& operator <<(std::ostream& os, const SE3& rhs){
    for(int i=0; i<3; i++){
	os << rhs.get_rotation().get_matrix()[i] << rhs.get_translation()[i] << std::endl;
    }
    return os;
}

// operator istream& >>
inline std::istream& operator>>(std::istream& is, SE3& rhs){
    for(int i=0; i<3; i++){
	is >> rhs.get_rotation().my_matrix[i] >> rhs.get_translation()[i];
    }
    return is;
}


//////////////////
// operator *   //
// SE3 * Vector //
//////////////////

template<class VectorType>
struct SE3VMult {
    inline static void eval(Vector<4>& ret, const SE3& lhs, const VectorType& rhs){
	ret.template slice<0,3>()=lhs.get_rotation()*rhs.template slice<0,3>();
	ret.template slice<0,3>()+=lhs.get_translation() * rhs[3];
	ret[3] = rhs[3];
    }
};

template<class Accessor> inline
Vector<4> operator*(const SE3& lhs, const FixedVector<4,Accessor>& rhs){
    return Vector<4>(lhs,rhs,Operator<SE3VMult<FixedVector<4, Accessor> > >());
}

template<class Accessor> inline
Vector<4> operator*(const SE3& lhs, const DynamicVector<Accessor>& rhs){
    //FIXME: size checking
    return Vector<4>(lhs,rhs,Operator<SE3VMult<DynamicVector<Accessor> > >());
}

template <class Accessor> inline
Vector<3> operator*(const SE3& lhs, const FixedVector<3,Accessor>& rhs){
    return lhs.get_rotation()*rhs + lhs.get_translation();
}


//////////////////
// operator *   //
// Vector * SE3 //
//////////////////

template<class Accessor>
struct VSE3Mult {
    inline static void eval(Vector<4>& ret, const FixedVector<4,Accessor>& lhs, const SE3& rhs){
	ret.template slice<0,3>() = lhs.template slice<0,3>() * rhs.get_rotation();
	ret[3] = lhs[3];
	ret[3] += lhs.template slice<0,3>() * rhs.get_translation();
    }
};

template<class Accessor> inline
Vector<4> operator*(const FixedVector<4,Accessor>& lhs, const SE3& rhs){
    return Vector<4>(lhs,rhs,Operator<VSE3Mult<Accessor> >());
}



//////////////////
// operator *   //
// SE3 * Matrix //
//////////////////

template <int RHS, class Accessor>
struct SE3MMult {
    inline static void eval(Matrix<4,RHS>& ret, const SE3& lhs, const FixedMatrix<4,RHS,Accessor>& rhs){
	for(int i=0; i<RHS; i++){
	    ret.T()[i].template slice<0,3>() = lhs.get_rotation() * rhs.T()[i].template slice<0,3>();
	    ret.T()[i].template slice<0,3>() += lhs.get_translation() * rhs(3,i);
	    ret(3,i) = rhs(3,i);
	}
    }
};


template <int RHS, class Accessor> inline 
Matrix<4,RHS> operator*(const SE3& lhs, const FixedMatrix<4,RHS,Accessor>& rhs){
    return Matrix<4,RHS>(lhs,rhs,Operator<SE3MMult<RHS,Accessor> >());
}


//////////////////
// operator *   //
// Matrix * SE3 //
//////////////////

template <int LHS, class Accessor> 
struct MSE3Mult {
    inline static void eval(Matrix<LHS,4>& ret, const FixedMatrix<LHS,4,Accessor>& lhs, const SE3& rhs){
	for(int i=0; i<LHS; i++){
	    ret[i].template slice<0,3>() = lhs[i].template slice<0,3>() * rhs.get_rotation();
	    ret(i,3) = rhs.get_translation() * lhs[i].template slice<0,3>();
	    ret(i,3) += lhs(i,3);
	}
    }
};


template <int LHS, class Accessor> inline 
Matrix<LHS,4> operator*(const FixedMatrix<LHS,4,Accessor>& lhs, const SE3& rhs){
    return Matrix<LHS,4>(lhs,rhs,Operator<MSE3Mult<LHS,Accessor> >());
}


namespace SE3static 
{
    static double zero[3]={0,0,0};
}

inline SE3::SE3() :
    my_translation(SE3static::zero)
{}


inline SE3 SE3::exp(const Vector<6>& mu){
    static const double one_6th = 1.0/6.0;
    static const double one_20th = 1.0/20.0;

    SE3 result;

    const Vector<3> w = mu.slice<3,3>();
    const double theta_sq = w*w;
    const double theta = sqrt(theta_sq);
    double A, B;
    
    const Vector<3> cross = w ^ mu.slice<0,3>();
    if (theta_sq < 1e-8) {
	A = 1.0 - one_6th * theta_sq;
	B = 0.5;
	result.get_translation() = mu.slice<0,3>() + 0.5 * cross;
    } else {
	double C;
	if (theta_sq < 1e-6) {
	    C = one_6th*(1.0 - one_20th * theta_sq);
	    A = 1.0 - theta_sq * C;
	    B = 0.5 - 0.25 * one_6th * theta_sq;
	} else {
	    const double inv_theta = 1.0/theta;
	    A = sin(theta) * inv_theta;
	    B = (1 - cos(theta)) * (inv_theta * inv_theta);
	    C = (1 - A) * (inv_theta * inv_theta);
	}
	result.get_translation() = mu.slice<0,3>() + B * cross + C * (w ^ cross);
    }
    rodrigues_so3_exp(w, A, B, result.my_rotation.my_matrix);
    return result;
}

inline Vector<6> SE3::ln(const SE3& se3) {
    Vector<3> rot = se3.my_rotation.ln();
    double theta = sqrt(rot*rot);
    double shtot = 0.5;
    
    if(theta > 0.00001) {
	shtot = sin(theta/2)/theta;
    }
    
    // now do the rotation
    Vector<3> halfrot = rot * -0.5;
    SO3 halfrotator = SO3::exp(halfrot);
    
    Vector<3> rottrans = halfrotator * se3.my_translation;
    
    if(theta > 0.001){
	rottrans -= rot * ((se3.my_translation * rot) * (1-2*shtot) / (rot*rot));
    } else {
	rottrans -= rot * ((se3.my_translation * rot)/24);
    }
    
    rottrans /= (2 * shtot);
    
    Vector<6> result;
    result.slice<0,3>()=rottrans;
    result.slice<3,3>()=rot;
    return result;
}

inline SE3 SE3::inverse() const {
    const SO3& rinv = my_rotation.inverse();
    return SE3(rinv, -(rinv*my_translation));
}

inline SE3& SE3::left_multiply_by(const SE3& left) {
    my_translation = left.my_translation + left.get_rotation() * my_translation;
    my_rotation = left.my_rotation * my_rotation;
    return *this;
}

inline Vector<4> SE3::generator_field(int i, Vector<4> pos){
    double result_d[]={0,0,0,0};
    Vector<4> result(result_d);
    if(i < 3){
	result[i]=pos[3];
	return result;
    }
    result[(i+1)%3] = - pos[(i+2)%3];
    result[(i+2)%3] = pos[(i+1)%3];
    return result;
}


inline SE3 operator*(const SO3& lhs, const SE3& rhs){
    SE3 result;
    result.my_rotation = lhs*rhs.my_rotation;
    result.my_translation = lhs*rhs.my_translation;
    return result;
}

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
