/*
     Copyright (C) 2007 The Authors

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

#ifndef __SO2_H
#define __SO2_H

#include <TooN/TooN.h>
#include <TooN/helpers.h>

#ifndef TOON_NO_NAMESPACE
namespace TooN {
#endif

template<typename Precision> class SO2;
template <typename Precision> class SE2;

template<typename Precision> inline std::istream & operator>>(std::istream &, SO2<Precision> & );
template<typename Precision> inline std::istream & operator>>(std::istream &, SE2<Precision> & );

template<typename Precision = double>
class SO2 {
	friend std::istream& operator>> <Precision>(std::istream&, SO2& );
	friend std::istream& operator>> <Precision>(std::istream&, SE2<Precision>& );

public:
	SO2() : my_matrix(Identity) {} 
  
	SO2(const Matrix<2,2,Precision>& rhs) {  *this = rhs; }

	SO2(const Precision l) { *this = exp(l); }
  
	template <int R, int C, typename P, typename A> 
	inline SO2& operator=(const Matrix<R,C,P,A>& rhs){
		my_matrix = rhs;
		coerce();
		return *this;
	}

	void coerce(){
		my_matrix[0] = unit(my_matrix[0]);
		my_matrix[1] -= my_matrix[0] * (my_matrix[0]*my_matrix[1]);
		my_matrix[1] = unit(my_matrix[1]);
	}

	inline static SO2 exp(const Precision & d){
		SO2<Precision> result;
		result.my_matrix[0][0] = result.my_matrix[1][1] = cos(d);
		result.my_matrix[1][0] = sin(d);
		result.my_matrix[0][1] = -result.my_matrix[1][0];
		return result;
	}

	Precision ln() const { return atan2(my_matrix[1][0], my_matrix[0][0]); }

	SO2 inverse() const { return SO2(*this, Invert()); }

	SO2& operator *=(const SO2& rhs){
		my_matrix=my_matrix*rhs.my_matrix;
		return *this;
	}
	SO2 operator *(const SO2& rhs) const { return SO2(*this,rhs); }

	const Matrix<2,2,Precision>& get_matrix() const {return my_matrix;}

	/// returns generator matrix
	static Matrix<2,2,Precision> generator() {
		Matrix<2,2,Precision> result;
		result[0] = makeVector(0,-1);
		result[1] = makeVector(1,0);
		return result;
	}

 private:
	struct Invert {};
	inline SO2(const SO2& so2, const Invert&) : my_matrix(so2.my_matrix.T()) {}
	inline SO2(const SO2& a, const SO2& b) : my_matrix(a.my_matrix*b.my_matrix) {}

	Matrix<2,2,Precision> my_matrix;
};

template <typename Precision>
inline std::ostream& operator<< (std::ostream& os, const SO2<Precision> & rhs){
	return os << rhs.get_matrix();
}

template <typename Precision>
inline std::istream& operator>>(std::istream& is, SO2<Precision>& rhs){
	return is >> rhs.my_matrix;
	rhs.coerce();
}

template<int D, typename P1, typename PV, typename Accessor>
inline Vector<2, typename Internal::MultiplyType<P1, PV>::type> operator*(const SO2<P1> & lhs, const Vector<D, PV, Accessor> & rhs){
	return lhs.get_matrix() * rhs;
}

template<int D, typename P1, typename PV, typename Accessor>
inline Vector<2, typename Internal::MultiplyType<PV,P1>::type> operator*(const Vector<D, PV, Accessor>& lhs, const SO2<P1> & rhs){
	return lhs * rhs.get_matrix();
}

template <int R, int C, typename P1, typename P2, typename Accessor> 
inline Matrix<2,C,typename Internal::MultiplyType<P1,P2>::type> operator*(const SO2<P1> & lhs, const Matrix<R,C,P2,Accessor>& rhs){
	return lhs.get_matrix() * rhs;
}

template <int R, int C, typename P1, typename P2, typename Accessor>
inline Matrix<R,2,typename Internal::MultiplyType<P1,P2>::type> operator*(const Matrix<R,C,P1,Accessor>& lhs, const SO2<P2>& rhs){
	return lhs * rhs.get_matrix();
}

#ifndef TOON_NO_NAMESPACE
}
#endif

#endif
