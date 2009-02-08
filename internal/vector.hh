//-*- c++ -*-

template<int Size=-1, typename Precision=double, typename Base=VBase<Size, Precision> >
class Vector : public Base {
public:
  // sneaky hack: only one of these constructors will work with any given base
  // class but they don't generate errors unless the user tries to use one of them
  // although the error message may be less than helpful - maybe this can be changed?
  inline Vector(){}
  inline Vector(Precision* data) : Base (data) {}
  inline Vector(int size_in) : Base(size_in) {}
  inline Vector(Precision* data_in, int size_in, int stride_in, Slicing) : Base(data_in, size_in, stride_in) {}
  inline Vector(Precision* data_in, int stride_in, Slicing) : Base(data_in, stride_in) {}


  // constructors to allow return value optimisations
  // construction from 1-ary operator
  template <class T, class Op>
  inline Vector(const T& arg, const Operator<Op>& op) : Base(arg,op) {
    Op::eval(*this,arg);
  }

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline Vector(const LHS& lhs, const RHS& rhs, const Operator<Op>& op)
    : Base(lhs,rhs,op) {
    Op::eval(*this,lhs,rhs);
  }

  // copy constructor listed explicitly
  inline Vector(const Vector<Size,Precision,Base>& from)
    : Base(from) {
    (*this)=from;
  }

  // constructor from arbitrary vector
  template<int Size2, typename Precision2, typename Base2>
  inline Vector(const Vector<Size2,Precision2,Base2>& from):
    Base(from) {
    operator=(from);
  }

  // operator = from copy
  inline Vector& operator= (const Vector& from){
    SizeMismatch<Size,Size>::test(Base::size(), from.size());
    const int s=Base::size();
    for(int i=0; i<s; i++){
      (*this)[i]=from[i];
    }
    return *this;
  }

  // operator =
  template<int Size2, typename Precision2, typename Base2>
  Vector<Size,Precision,Base >& operator= (const Vector<Size2, Precision2, Base2>& from){
    SizeMismatch<Size,Size2>::test(Base::size(), from.size());
    const int s=Base::size();
    for(int i=0; i<s; i++){
      (*this)[i]=from[i];
    }
    return *this;
  }

};
