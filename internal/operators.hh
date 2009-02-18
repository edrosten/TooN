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

	//Mini operators for passing to Pairwise, etc
	struct Add{ template<class A, class B, class C>      static A op(const B& b, const C& c){return b+c;} };
	struct Subtract{ template<class A, class B, class C> static A op(const B& b, const C& c){return b-c;} };
	struct Multiply{ template<class A, class B, class C> static A op(const B& b, const C& c){return b*c;} };
	struct Divide{ template<class A, class B, class C>   static A op(const B& b, const C& c){return b/c;} };
	
	//Automatic type deduction of return types
	template<class L, class R> struct AddType {      typedef TOON_TYPEOF((L()+R())) type; };
	template<class L, class R> struct SubtractType { typedef TOON_TYPEOF((L()-R())) type; };
	template<class L, class R> struct MultiplyType { typedef TOON_TYPEOF((L()*R())) type; };
	template<class L, class R> struct DivideType {   typedef TOON_TYPEOF((L()*R())) type; };
	
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
	return Vector<Internal::Sizer<S1,S2>::size,restype>(v1, v2, Operator<Internal::Pairwise<restype, Internal::Add> >(), v1.size());
}

// Addition Vector - Vector
template<int S1, int S2, typename P1, typename P2, typename B1, typename B2> 
Vector<Internal::Sizer<S1,S2>::size, typename Internal::SubtractType<P1, P2>::type> operator-(const Vector<S1, P1, B1>& v1, const Vector<S2, P2, B2>& v2)
{
	typedef typename Internal::SubtractType<P1, P2>::type restype;
	SizeMismatch<S1, S2>:: test(v1.size(),v2.size());
	return Vector<Internal::Sizer<S1,S2>::size,restype>(v1, v2, Operator<Internal::Pairwise<restype, Internal::Subtract> >(), v1.size());
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



// Addition Matrix + Matrix
template<int R1, int C1, int R2, int C2, typename P1, typename P2, typename B1, typename B2> 
Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size, typename Internal::AddType<P1, P2>::type> operator+(const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
{
	typedef typename Internal::AddType<P1, P2>::type restype;
	SizeMismatch<R1, R2>:: test(m1.num_rows(),m2.num_rows());
	SizeMismatch<C1, C2>:: test(m1.num_cols(),m2.num_cols());
	return Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size,restype>(m1, m2, Operator<Internal::Pairwise<restype, Internal::Add> >(), m1.num_rows(), m1.num_cols());
}

// Addition Matrix - Matrix
template<int R1, int C1, int R2, int C2, typename P1, typename P2, typename B1, typename B2> 
Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size, typename Internal::SubtractType<P1, P2>::type> operator-(const Matrix<R1, C1, P1, B1>& m1, const Matrix<R2, C2, P2, B2>& m2)
{
	typedef typename Internal::SubtractType<P1, P2>::type restype;
	SizeMismatch<R1, R2>:: test(m1.num_rows(),m2.num_rows());
	SizeMismatch<C1, C2>:: test(m1.num_cols(),m2.num_cols());
	return Matrix<Internal::Sizer<R1,R2>::size, Internal::Sizer<C1,C2>::size,restype>(m1, m2, Operator<Internal::Pairwise<restype, Internal::Subtract> >(), m1.num_rows(), m1.num_cols());
}


////////////////////////////////////////////////////////////////////////////////
//
// vector <op> scalar
//

#define TOON_MAKE_SCALAR_OP_PAIR(OPNAME, OP) \
template<int S, typename P1, typename B1, typename P2> \
Vector<S, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const Vector<S, P1, B1>& v, const P2& s)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Vector<S,restype>(v, s, Operator<Internal::ApplyScalar<restype, Internal::OPNAME> >(), v.size());\
}\
\
template<int S, typename P1, typename P2, typename B2> \
Vector<S, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const P1& s, const Vector<S, P2, B2>& v)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Vector<S,restype>(s, v, Operator<Internal::ApplyScalarLeft<restype, Internal::OPNAME> >(), v.size());\
}\
template<int R, int C, typename P1, typename B1, typename P2> \
Matrix<R, C, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const Matrix<R, C, P1, B1>& m, const P2& s)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Matrix<R, C,restype>(m, s, Operator<Internal::ApplyScalar<restype, Internal::OPNAME> >(), m.num_rows(), m.num_cols());\
}\
\
template<int R, int C, typename P1, typename P2, typename B2> \
Matrix<R, C, typename Internal::OPNAME##Type<P1, P2>::type> operator OP (const P1& s, const Matrix<R, C, P2, B2>& m)\
{	\
	typedef typename Internal::OPNAME##Type<P1, P2>::type restype;\
	return Matrix<R, C,restype>(s, m, Operator<Internal::ApplyScalarLeft<restype, Internal::OPNAME> >(), m.num_rows(), m.num_cols());\
}

TOON_MAKE_SCALAR_OP_PAIR(Add, +)
TOON_MAKE_SCALAR_OP_PAIR(Add, -)
TOON_MAKE_SCALAR_OP_PAIR(Add, *)
TOON_MAKE_SCALAR_OP_PAIR(Add, /)

#undef TOON_MAKE_SCALAR_OP_PAIR




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



