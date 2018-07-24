//-*- c++ -*-
//
// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND OTHER CONTRIBUTORS ``AS IS''
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR OTHER CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.



namespace TooN {

	namespace Internal
	{
		//Dummy struct for Diagonal operator
		template<int Size, typename Precision, typename Base>
		struct DiagMatrixOp;
	}

	template<int Size, typename Precision, typename Base>
	struct Operator<Internal::DiagMatrixOp<Size, Precision, Base> >
	{
		public:
		///@name Constructors
		///@{


		inline Operator() {}
		inline Operator(int size_in) : my_vector(size_in) {}
		inline Operator(Precision* data) : my_vector(data) {}
		inline Operator(Precision* data, int size_in) : my_vector(data,size_in) {}
		inline Operator(Precision* data_in, int size_in, int stride_in, Internal::Slicing)
			: my_vector(data_in, size_in, stride_in, Internal::Slicing() ) {}

		// constructors to allow return value optimisations
		// construction from 0-ary operator
		///my_vector constructed from a TooN::Operator 
		template <class Op>
		inline Operator(const Operator<Op>& op)
			: my_vector (op)
		{
			op.eval(my_vector);
		}
		
		// constructor from arbitrary vector
		template<int Size2, typename Precision2, typename Base2>
		inline Operator(const Vector<Size2,Precision2,Base2>& from)
			: my_vector(from.size())
		{
			my_vector=from;
		}
		///@}


		///@name Operator members
		///@{
		template<int R, int C, class P, class B>
		void eval(Matrix<R,C,P,B>& m) const {
			SizeMismatch<Size, Size>::test(m.num_rows(), m.num_cols());
			m = Zeros;
			m.diagonal_slice() = my_vector;
		}

		///The vector used to hold the leading diagonal.
		Vector<Size,Precision,Base> my_vector;
	};

/**
@class DiagonalMatrix 
A diagonal matrix

Support is limited but diagonal matrices can be multiplied by vectors, matrices
or diagonal matrices on either side.

Diagonal matrices can be created from vectors by using the <code> as_diagonal() 
</code> member function:

@code
Vector<3> v = makeVector(1,2,3);
Vector<3> v2 = v.as_diagonal() * v;   // v2 = (1,4,9)
@endcode

A vector can be obtained from the diagonal matrix by using the
<code> diagonal_slice() </code> member function.
@ingroup gLinAlg
 **/
template<int Size=Dynamic, typename Precision=DefaultPrecision, typename Base=Internal::VBase>
struct DiagonalMatrix: public Operator<Internal::DiagMatrixOp<Size, Precision, Base> > {
public:
	///@name Constructors
	///@{
	
	inline DiagonalMatrix() {}
	inline DiagonalMatrix(int size_in) : Operator<Internal::DiagMatrixOp<Size, Precision, Base> >(size_in) {}
	inline DiagonalMatrix(Precision* data) : Operator<Internal::DiagMatrixOp<Size, Precision, Base> >(data) {}
	inline DiagonalMatrix(Precision* data, int size_in) : Operator<Internal::DiagMatrixOp<Size, Precision, Base> >(data,size_in) {}
	inline DiagonalMatrix(Precision* data_in, int size_in, int stride_in, Internal::Slicing)
		: Operator<Internal::DiagMatrixOp<Size, Precision, Base> >(data_in, size_in, stride_in, Internal::Slicing() ) {}

	// constructors to allow return value optimisations
	// construction from 0-ary operator
	///my_vector constructed from a TooN::Operator 
	template <class Op>
	inline DiagonalMatrix(const Operator<Op>& op)
		: Operator<Internal::DiagMatrixOp<Size, Precision, Base> > (op)
	{
		op.eval(this->my_vector);
	}
	
	// constructor from arbitrary vector
	template<int Size2, typename Precision2, typename Base2>
	inline DiagonalMatrix(const Vector<Size2,Precision2,Base2>& from)
		: Operator<Internal::DiagMatrixOp<Size, Precision, Base> >(from.size())
	{
		this->my_vector=from;
	}
	///@}



	///Index the leading elements on the diagonal 
	Precision& operator[](int i){return this->my_vector[i];}
	///Index the leading elements on the diagonal 
	const Precision& operator[](int i) const {return this->my_vector[i];}
	
	///Return the leading diagonal as a vector.
	typename Vector<Size, Precision, Base>::as_slice_type diagonal_slice() {
		return this->my_vector.as_slice();
	}

	///Return the leading diagonal as a vector.
	const typename Vector<Size, Precision, Base>::as_slice_type diagonal_slice() const {
		return this->my_vector.as_slice();
	}

	DiagonalMatrix<Size, Precision> operator-() const
	{
		return -this->my_vector;
	}

	DiagonalMatrix<Size, Precision> inverse() const
	{
		Vector<Size, Precision> inv(this->my_vector.size());
		for(int i=0; i < inv.size(); i++)
			inv[i] = 1 / this->my_vector[i];
		return inv;
	}
};


template<int S1, typename P1, typename B1, int S2, typename P2, typename B2>
inline Vector<Internal::Sizer<S1,S2>::size, typename Internal::MultiplyType<P1,P2>::type>
operator*(const DiagonalMatrix<S1,P1,B1>& d, const Vector<S2,P2,B2>& v){
	return diagmult(d.my_vector,v);
}

template<int S1, typename P1, typename B1, int S2, typename P2, typename B2>
inline Vector<Internal::Sizer<S1,S2>::size, typename Internal::MultiplyType<P1,P2>::type>
operator*( const Vector<S1,P1,B1>& v, const DiagonalMatrix<S2,P2,B2>& d){
	return diagmult(v,d.my_vector);
}

// perhaps not the safest way to do this as we're returning the same operator used to normally make vectors
template<int S1, typename P1, typename B1, int S2, typename P2, typename B2>
inline DiagonalMatrix<Internal::Sizer<S1,S2>::size, typename Internal::MultiplyType<P1,P2>::type>
operator*( const DiagonalMatrix<S1,P1,B1>& d1, const DiagonalMatrix<S2,P2,B2>& d2){
	SizeMismatch<S1,S2>::test(d1.my_vector.size(),d2.my_vector.size());
	return Operator<Internal::VPairwise<Internal::Multiply,S1,P1,B1,S2,P2,B2> >(d1.my_vector,d2.my_vector);
}

template<int R, int C, int Size, typename P1, typename P2, typename B1, typename B2>
Matrix<R, C, typename Internal::MultiplyType<P1,P2>::type>
operator* (const Matrix<R, C, P1, B1>& m, const DiagonalMatrix<Size, P2, B2>& d){
	return diagmult(m,d.my_vector);
}

template<int R, int C, typename P1, typename B1, int Size, typename P2, typename B2> 
Matrix<R, C, typename Internal::MultiplyType<P1,P2>::type>
operator* (const DiagonalMatrix<Size,P1,B1>& d, const Matrix<R,C,P2,B2>& m)
{
	return diagmult(d.my_vector, m);
}

}
