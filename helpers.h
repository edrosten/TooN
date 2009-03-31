
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

		struct Zero;
		struct SizedZero;
		struct RCZero;
		struct Identity;
		struct SizedIdentity;
	
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


	template<> struct Operator<Internal::RCZero> {
	Operator(int r, int c) : my_rows(r), my_cols(c) {}

		const int my_rows;
		const int my_cols;

		int num_rows() const {return my_rows;}
		int num_cols() const {return my_cols;}

		template<int R, int C, class P, class B>
			void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}
	};


	template<> struct Operator<Internal::SizedZero> {

		// no idea why this doesn't indent properly
		Operator(int s)	: my_size(s) {}
		
		const int my_size;
		

		int size() const {return my_size;}
		int num_rows() const {return my_size;}
		int num_cols() const {return my_size;}

		template<int Size, class Precision, class Base>
			void eval(Vector<Size, Precision, Base>& v) const {
			for(int i=0; i < v.size(); i++) {
				v[i]= 0;
			}
		}

		template<int R, int C, class P, class B>
			void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}
		
	};


	template<> struct Operator<Internal::Zero> {
		template<int Size, class Precision, class Base>
			void eval(Vector<Size, Precision, Base>& v) const {
			for(int i=0; i < v.size(); i++) {
				v[i]= 0;
			}
		}

		template<int R, int C, class P, class B>
			void eval(Matrix<R,C,P,B>& m) const {
			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
		}

		Operator<Internal::SizedZero> operator()(int s){
			return Operator<Internal::SizedZero>(s);
		}

		Operator<Internal::RCZero> operator()(int r, int c){
			return Operator<Internal::RCZero>(r,c);
		}

	};

	template<> struct Operator<Internal::SizedIdentity> {
		// no idea why this doesn't indent properly
	Operator(int s)	: my_size(s) {}
		
		const int my_size;
		
		int num_rows() const {return my_size;}
		int num_cols() const {return my_size;}

		template<int R, int C, class P, class B>
			void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());

			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
						
			for(int r=0; r < m.num_rows(); r++) {
				m(r,r) = 1;
			}
		}
	};
		


	template<> struct Operator<Internal::Identity> {

		template<int R, int C, class P, class B>
			void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<R, C>::test(m.num_rows(), m.num_cols());

			for(int r=0; r<m.num_rows(); r++){
				for(int c=0; c<m.num_cols(); c++){
					m(r,c)=0;
				}
			}
						
			for(int r=0; r < m.num_rows(); r++) {
				m(r,r) = 1;
			}
		}

		Operator<Internal::SizedIdentity> operator()(int s){
			return Operator<Internal::SizedIdentity>(s);
		}
	};


	static Operator<Internal::Zero> Zero;
	static Operator<Internal::Identity> Identity;


}
#endif
