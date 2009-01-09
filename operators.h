//-*- c++ -*-

//////////////////////////////////////////////////////////////////////////////////////////////
//                                       Operators
//////////////////////////////////////////////////////////////////////////////////////////////


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

// operator /= double
template<int Size, typename Precision, typename Base>
Vector<Size, Precision, Base>& operator /= (Vector<Size, Precision, Base>& lhs, double rhs){
  const int s=lhs.size();
  for(int i=0; i<s; i++){
    lhs[i]/=rhs;
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


