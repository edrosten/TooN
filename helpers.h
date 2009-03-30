
/*                       
	 Copyright (C) 2005,2009 Tom Drummond, E. Rosten

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
#ifndef TOON_INCLUDE_HELPERS_H
#define TOON_INCLUDE_HELPERS_H

#include <TooN/TooN.h>
#include <cmath>

namespace TooN {


template<int Size, class Precision, class Base> void Fill(Vector<Size, Precision, Base>& v, const Precision& p)
{
	for(int i=0; i < v.size(); i++)
			v[i]= p;
}

template<int Rows, int Cols, class Precision, class Base> void Fill(Matrix<Rows, Cols, Precision, Base>& m, const Precision& p)
{
	for(int i=0; i < m.num_rows(); i++)
		for(int j=0; j < m.num_cols(); j++)
			m[i][j] = p;
}

template<int Size, class Precision, class Base> inline Vector<Size, Precision> unit(const Vector<Size, Precision, Base> & v)
{
	using std::sqrt;
	return v/sqrt(v*v);
}

namespace Internal{

	struct Zero {
		template<int R, int C, class P, class B>
			static void eval(Matrix<R, C, P, B>& m) {
			for(int r=0; r < m.num_rows(); r++) {
				for(int c=0; c < m.num_rows(); c++) {
					m[r][c] = 0;
				}
			}
		}

		template<int Size, class Precision, class Base>
			static void eval(Vector<Size, Precision, Base>& v) {
			for(int i=0; i < v.size(); i++) {
				v[i]= 0;
			}
		}
	};

	struct Identity {
		template<int R, int C, class P, class B> static void eval(Matrix<R, C, P, B>& m) {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());
			
			for(int r=0; r < m.num_rows(); r++) {
				for(int c=0; c < m.num_rows(); c++) {
					m[r][c] = 0;
				}
			}
			
			for(int i=0; i < m.num_rows(); i++) {
				m[i][i] = 1;
			}
		}
	};
	
	struct Copy
	{
		template<int R, int C, class P, class B, class Data> static void eval(Matrix<R, C, P, B>& m, const Data * data)
		{
			for(int r=0; r < m.num_rows(); r++)
				for(int c=0; c < m.num_rows(); c++)
					m[r][c] = *data++;
		}
	};

}


template<> class Operator<Internal::Zero> {
 public:
	template<int Size, class Precision, class Base>
		void eval(Vector<Size, Precision, Base>& v) const {
		for(int i=0; i < v.size(); i++) {
			v[i]= 0;
		}
	}
 };


static NoAliasOperator<Internal::Zero> Zero;
static Operator<Internal::Identity> Identity;


}
#endif
