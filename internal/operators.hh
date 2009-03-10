//-*- c++ -*-
// Operations:
//  vector+vector
//  vector-vector
//  vector dot evctor
//  vector * constant
//  vector / constant
//  vector +

//////////////////////////////////////////////////////////////////////////////////////////////
//                                       Operators
//////////////////////////////////////////////////////////////////////////////////////////////

template<class Op> struct Operator{};

namespace Internal{
	
	//Operator classes. These are evaluated in the constructor
	//of Vector = Vector+Vector, so make use of return value optimization.
	template<typename Precision, typename Op> struct Pairwise
	{
		template<int S, typename B, int S1, typename P1, typename B1, int S2, typename P2, typename B2> 
		static void eval(Vector<S, Precision, B>& res, const Vector<S1, P1, B1>& v1, const Vector<S2, P2, B2>& v2)
		{
			for(int i=0; i < res.size(); ++i)
				res[i] = Op::template op<Precision,P1, P2>(v1[i],v2[i]);
		}

		template<int R, int C, typename B, int R1, int C1, typename P1, typename B1, int R2, int C2, typename P2, typename B2> 
		static void eval(Matrix<R, C, Precision, B>& res, const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
		{
			for(int i=0; i < res.num_rows(); ++i)
				for(int j=0; j < res.num_cols(); ++j)
				res[i][j] = Op::template op<Precision,P1, P2>(m1[i][j],m2[i][j]);
		}
	};

	template<typename Precision, typename Op> struct ApplyScalar 
	{
		template<int S, typename B, int S1, typename P1, typename B1, typename P2>
		static void eval(Vector<S, Precision, B>& res, const Vector<S1, P1, B1>& v, const P2& s)
		{
			for(int i=0; i < res.size(); ++i)
				res[i] = Op::template op<Precision,P1, P2>(v[i],s);
		}

		template<int R, int C, typename B, int R1, int C1, typename P1, typename B1, typename P2>
		static void eval(Matrix<R, C, Precision, B>& res, const Matrix<R1, C1, P1, B1>& m, const P2& s)
		{		
		
			for(int i=0; i < res.num_rows(); ++i)
				for(int j=0; j < res.num_cols(); ++j)
					res[i][j] = Op::template op<Precision,P1, P2>(m[i][j],s);
		}
	};

	template<typename Precision, typename Op> struct ApplyScalarLeft
	{
		template<int S, typename B, int S2, typename P1, typename P2, typename B2>
		static void eval(Vector<S, Precision, B>& res, const P1& s, const Vector<S2, P2, B2>& v)
		{
			for(int i=0; i < res.size(); ++i)
				res[i] = Op::template op<Precision,P1, P2>(s, v[i]);
		}

		template<int R, int C, typename B, int R2, int C2, typename P1, typename B2, typename P2>
		static void eval(Matrix<R, C, Precision, B>& res, const P1& s, const Matrix<R2, C2, P2, B2>& m)
		{		
		
			for(int i=0; i < res.num_rows(); ++i)
				for(int j=0; j < res.num_cols(); ++j)
					res[i][j] = Op::template op<Precision,P1, P2>(s, m[i][j]);
		}
	};

	//FIXME what about BLAS?
	struct MatrixMultiply
	{
		template<int R, int C, typename Precision, typename B, int R1, int C1, typename P1, typename B1, int R2, int C2, typename P2, typename B2> 
		static void eval(Matrix<R, C, Precision, B>& res, const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
		{
			for(int i=0; i < res.num_rows(); ++i)
				for(int j=0; j < res.num_cols(); ++j)
					res[i][j] = m1[i] * (m2.T()[j]);
		}
	};

	struct MatrixVectorMultiply
	{
		template<int Sout, typename Pout, typename Bout, int R, int C, int Size, typename P1, typename P2, typename B1, typename B2>
		static void eval(Vector<Sout, Pout, Bout>& res, const Matrix<R, C, P1, B1>& m, const Vector<Size, P2, B2>& v)
		{
			for(int i=0; i < res.size(); ++i){
				res[i] = m[i] * v;
			}
		}
	};

	// this is distinct to cater for non communing precision types
	struct VectorMatrixMultiply
	{
		template<int Sout, typename Pout, typename Bout, int R, int C, int Size, typename P1, typename P2, typename B1, typename B2>
		static void eval(Vector<Sout, Pout, Bout>& res, const Vector<Size, P2, B2>& v, const Matrix<R, C, P1, B1>& m)
		{
			for(int i=0; i < res.size(); ++i){
				res[i] = v * m[i];
			}
		}
	};


	//Mini operators for passing to Pairwise, etc
	struct Add{ template<class A, class B, class C>      static A op(const B& b, const C& c){return b+c;} };
	struct Subtract{ template<class A, class B, class C> static A op(const B& b, const C& c){return b-c;} };
	struct Multiply{ template<class A, class B, class C> static A op(const B& b, const C& c){return b*c;} };
	struct Divide{ template<class A, class B, class C>   static A op(const B& b, const C& c){return b/c;} };
	
	//Automatic type deduction of return types
	template<class L, class R> struct AddType {      typedef TOON_TYPEOF( (*static_cast<L*>(0) + *static_cast<R*>(0))) type;};
	template<class L, class R> struct SubtractType { typedef TOON_TYPEOF( (*static_cast<L*>(0) - *static_cast<R*>(0))) type;};
	template<class L, class R> struct MultiplyType { typedef TOON_TYPEOF( (*static_cast<L*>(0) * *static_cast<R*>(0))) type;};
	template<class L, class R> struct DivideType   { typedef TOON_TYPEOF( (*static_cast<L*>(0) / *static_cast<R*>(0))) type;};
	
	//Output size, given input size. Be static if possible.
	template<int i, int j> struct Sizer{static const int size=i;};
	template<int i> struct Sizer<-1, i>{static const int size=i;};
	template<int i> struct Sizer<i, -1>{static const int size=i;};
	template<> struct Sizer<-1, -1>    {static const int size=-1;};
}

////////////////////////////////////////////////////////////////////////////////
//
// vector <op> vector
//


// Addition Vector + Vector
template<int S1, int S2, typename P1, typename P2, typename B1, typename B2> 
Vector<Internal::Sizer<S1,S2>::size, typename Internal::AddType<P1, P2>::type> operator+(const Vector<S1, P1, B1>& v1, const Vector<S2, P2, B2>& v2)
{
	typedef typename Internal::AddType<P1, P2>::type restype;
	SizeMismatch<S1, S2>:: test(v1.size(),v2.size());
	return Vector<Internal::Sizer<S1,S2>::size,restype>(v1, v2, v1.size(), Operator<Internal::Pairwise<restype, Internal::Add> >());
}

// Addition Vector - Vector
template<int S1, int S2, typename P1, typename P2, typename B1, typename B2> 
Vector<Internal::Sizer<S1,S2>::size, typename Internal::SubtractType<P1, P2>::type> operator-(const Vector<S1, P1, B1>& v1, const Vector<S2, P2, B2>& v2)
{
	typedef typename Internal::SubtractType<P1, P2>::type restype;
	SizeMismatch<S1, S2>:: test(v1.size(),v2.size());
	return Vector<Internal::Sizer<S1,S2>::size,restype>(v1, v2, v1.size(), Operator<Internal::Pairwise<restype, Internal::Subtract> >());
}

// Dot product Vector * Vector
template<int Size1, typename Precision1, typename Base1, int Size2, typename Precision2, typename Base2>
typename Internal::MultiplyType<Precision1, Precision2>::type operator*(const Vector<Size1, Precision1, Base1>& v1, const Vector<Size2, Precision2, Base2>& v2){
  SizeMismatch<Size1, Size2>:: test(v1.size(),v2.size());
  const int s=v1.size();
  typename Internal::MultiplyType<Precision1, Precision2>::type result=0;
  for(int i=0; i<s; i++){
    result+=v1[i]*v2[i];
  }
  return result;
}

template <typename P1, typename P2>
Vector<3, typename Internal::MultiplyType<P1,P2>::type> operator^(const Vector<3,P1>& v1, const Vector<3,P2>& v2){
	// assume the result of adding two restypes is also a restype
	typedef typename Internal::MultiplyType<P1,P2>::type restype;

	Vector<3, restype> result;

	result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];

	return result;
}




////////////////////////////////////////////////////////////////////////////////
//
// matrix <op> matrix
//

// Addition Matrix + Matrix
template<int R1, int C1, int R2, int C2, typename P1, typename P2, typename B1, typename B2> 
Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size, typename Internal::AddType<P1, P2>::type> operator+(const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
{
	typedef typename Internal::AddType<P1, P2>::type restype;
	SizeMismatch<R1, R2>:: test(m1.num_rows(),m2.num_rows());
	SizeMismatch<C1, C2>:: test(m1.num_cols(),m2.num_cols());
	return Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size,restype>(m1, m2, m1.num_rows(), m1.num_cols(), Operator<Internal::Pairwise<restype, Internal::Add> >());
}

// Addition Matrix - Matrix
template<int R1, int C1, int R2, int C2, typename P1, typename P2, typename B1, typename B2> 
Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size, typename Internal::SubtractType<P1, P2>::type> operator-(const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
{
	typedef typename Internal::SubtractType<P1, P2>::type restype;
	SizeMismatch<R1, R2>:: test(m1.num_rows(),m2.num_rows());
	SizeMismatch<C1, C2>:: test(m1.num_cols(),m2.num_cols());
	return Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size,restype>(m1, m2, m1.num_rows(), m1.num_cols(), Operator<Internal::Pairwise<restype, Internal::Subtract> >());
}

// Matrix multiplication Matrix * Matrix

template<int R1, int C1, int R2, int C2, typename P1, typename P2, typename B1, typename B2> 
Matrix<R1, C2, typename Internal::MultiplyType<P1, P2>::type> operator*(const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
{
	SizeMismatch<R1, C2>:: test(m1.num_rows(),m2.num_cols());
	SizeMismatch<C1, R2>:: test(m1.num_cols(),m2.num_rows());
	return Matrix<R1, C2, typename Internal::MultiplyType<P1, P2>::type>(m1, m2, m1.num_rows(), m2.num_cols(), Operator<Internal::MatrixMultiply>());
}

////////////////////////////////////////////////////////////////////////////////
//
// matrix <op> vector and vv.
//

// Matrix Vector multiplication Matrix * Vector

template<int R, int C, int Size, typename P1, typename P2, typename B1, typename B2>
Vector<R, typename Internal::MultiplyType<P1,P2>::type> operator*(const Matrix<R, C, P1, B1>& m, const Vector<Size, P2, B2>& v)
{
	SizeMismatch<C,Size>::test(m.num_cols(), v.size());
	return Vector<R, typename Internal::MultiplyType<P1,P2>::type> (m, v, m.num_rows(), Operator<Internal::MatrixVectorMultiply>() );
}
																	
// Vector Matrix multiplication Vector * Matrix

template<int Size, int R, int C, typename P1, typename P2, typename B1, typename B2>
Vector<C, typename Internal::MultiplyType<P1,P2>::type> operator*(const Vector<Size, P1, B1>& v, const Matrix<R, C, P2, B2>& m)
{
	SizeMismatch<R,Size>::test(m.num_rows(), v.size());
	return Vector<C, typename Internal::MultiplyType<P1,P2>::type> (v, m, m.num_cols(), Operator<Internal::VectorMatrixMultiply>() );
}


////////////////////////////////////////////////////////////////////////////////
//
// vector <op> scalar 
// scalar <op> vector 
// matrix <op> scalar 
// scalar <op> matrix 
//
// Except <scalar> / <matrix|vector> does not exist

// scalar on the right
#define TOON_MAKE_SCALAR_OP_RIGHT(OPNAME, OP) \
template<int R, int C, typename P1, typename B1, typename P2> \
Matrix<R, C, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const Matrix<R, C, P1, B1>& m, const P2& s)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Matrix<R, C,restype>(m, s, m.num_rows(), m.num_cols(), Operator<Internal::ApplyScalar<restype, Internal::OPNAME> >()); \
}\
\
template<int S, typename P1, typename B1, typename P2> \
Vector<S, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const Vector<S, P1, B1>& v, const P2& s)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Vector<S,restype>(v, s, v.size(), Operator<Internal::ApplyScalar<restype, Internal::OPNAME> >());\
}

// scalar on the left
#define TOON_MAKE_SCALAR_OP_LEFT(OPNAME, OP) \
template<int R, int C, typename P1, typename P2, typename B2> \
Matrix<R, C, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const P1& s, const Matrix<R, C, P2, B2>& m)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Matrix<R, C,restype>(s, m, m.num_rows(), m.num_cols(), Operator<Internal::ApplyScalarLeft<restype, Internal::OPNAME> >());\
} \
\
template<int S, typename P1, typename P2, typename B2> \
Vector<S, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const P1& s, const Vector<S, P2, B2>& v)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Vector<S,restype>(s, v, v.size(), Operator<Internal::ApplyScalarLeft<restype, Internal::OPNAME> >());\
}




#define TOON_MAKE_SCALAR_OPS(OPNAME, OP)\
TOON_MAKE_SCALAR_OP_LEFT(OPNAME, OP)\
TOON_MAKE_SCALAR_OP_RIGHT(OPNAME, OP)

TOON_MAKE_SCALAR_OPS(Add, +)
TOON_MAKE_SCALAR_OPS(Subtract, -)
TOON_MAKE_SCALAR_OPS(Multiply, *)

TOON_MAKE_SCALAR_OP_RIGHT(Divide, /)

#undef TOON_MAKE_SCALAR_OPS
#undef TOON_MAKE_SCALAR_OP_LEFT
#undef TOON_MAKE_SCALAR_OP_RIGHT





////////////////////////////////////////////////////////////////////////////////
//
// Stream I/O operators
//

// output operator <<
template <int Size, typename Precision, typename Base>
inline std::ostream& operator<< (std::ostream& os, const Vector<Size,Precision,Base>& v){
  for(int i=0; i<v.size(); i++){
    os << v[i] << " ";
  }
  return os;
}


template<int Rows, int Cols, typename Precision, class Base>
inline std::ostream& operator<< (std::ostream& os, const Matrix<Rows, Cols, Precision, Base>& m){
	for(int i=0; i < m.num_rows(); i++)
	{
		for(int j=0; j < m.num_cols(); j++)
		{
			if(j != 0)
				os << " ";
			os << m(i,j);
		}
		os << std::endl;
	}
	return os;
}



