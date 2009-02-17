//-*- c++ -*-

//////////////////////////////////////////////////////////////////////////////////////////////
//                                       Operators
//////////////////////////////////////////////////////////////////////////////////////////////

template<class Op> struct Operator{};

namespace Internal{

	template<typename Precision> struct Add
	{
		template<int S, typename B, int S1, typename B1, int S2, typename B2> 
		static void eval(Vector<S, Precision, B>& res, const Vector<S1, Precision, B1>& v1, const Vector<S2, Precision, B2>& v2)
		{
			SizeMismatch<S1, S2>:: test(v1.size(),v2.size());
			for(int i=0; i < res.size(); ++i)
				res[i] = v1[i] + v2[i];
		}

		template<int S1, typename B1, int S2, typename B2> 
		static int size(const Vector<S1, Precision, B1>& v1, const Vector<S2, Precision, B2>& v2)
		{
			SizeMismatch<S1, S2>:: test(v1.size(),v2.size());
			return v1.size();
		}
	};
}


// Addition Vector + Vector
template<int Size, typename Precision, typename B1, typename B2> 
Vector<Size, Precision> operator+(const Vector<Size, Precision, B1>& v1, const Vector<Size, Precision, B2>& v2)
{
	return Vector<Size, Precision>(v1, v2, Operator<Internal::Add<double> >());
}


// Dot product Vector * Vector
template<int Size1, typename Precision1, typename Base1,
	 int Size2, typename Precision2, typename Base2>
double operator*(const Vector<Size1, Precision1, Base1>& v1,
		 const Vector<Size2, Precision2, Base2>& v2){
  SizeMismatch<Size1, Size2>:: test(v1.size(),v2.size());
  const int s=v1.size();
  double result=0;
  for(int i=0; i<s; i++){
    result+=v1[i]*v2[i];
  }
  return result;
}

// operator += Vector
template<int Size1, typename Precision1, typename Base1,
	 int Size2, typename Precision2, typename Base2>
Vector<Size1, Precision1, Base1>& operator += (Vector<Size1, Precision1, Base1>& lhs,
					       const Vector<Size2, Precision2, Base2>& rhs){
  SizeMismatch<Size1,Size2>::test(lhs.size(),rhs.size());
  const int s=lhs.size();
  for(int i=0; i<s; i++){
    lhs[i]+=rhs[i];
  }
  return lhs;
}

// operator -= Vector
template<int Size1, typename Base1, typename Precision1,
	 int Size2, typename Precision2, typename Base2>
Vector<Size1, Precision1, Base1>& operator -= (Vector<Size1, Precision1, Base1>& lhs,
					       const Vector<Size2, Precision2, Base2>& rhs){
  SizeMismatch<Size1,Size2>::test(lhs.size(),rhs.size());
  const int s=lhs.size();
  for(int i=0; i<s; i++){
    lhs[i]-=rhs[i];
  }
  return lhs;
}

// operator *= double
template<int Size, typename Precision, typename Base>
Vector<Size, Precision, Base>& operator *= (Vector<Size, Precision, Base>& lhs, double rhs){
  const int s=lhs.size();
  for(int i=0; i<s; i++){
    lhs[i]*=rhs;
  }
  return lhs;
}

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



