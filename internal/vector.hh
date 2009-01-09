//-*- c++ -*-
// forward declarations
template<int Size, typename Base>
class Vector;

template<int Size, int Stride>
class SVBase;

// forward declaration
template<typename T>
class Operator;

// VBase owns its data
// Size is the static number of elements in the vector
// Stride is the static gap between elements
// Type is a hack so that x?y:z expressions work
// Type=0 means the data is internal to the object
// Type=1 means the data is external to the object
template <int Size, int Type>
class VBase;

class SDVBase;


template<int Size>
class VBase<Size, 0> {
public:
  inline VBase(){}
  inline VBase(const VBase& from){}
  
  // construction from 1-ary operator
  template <class T, class Op>
  inline VBase(const T&, const Operator<Op>&){}

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline VBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&){}
  
  // constructor from arbitrary vector
  template<int Size2, class Base2>
  inline VBase(const Vector<Size2,Base2>& from){}

  double* data(){return my_data;}
  const double* data() const {return my_data;}
  
  static int size(){return Size;}
  static int stride(){return 1;}

  double& operator[](int i){
    Internal::check_index(Size, i);
    return my_data[i];
  }
  const double& operator[](int i) const {
    Internal::check_index(Size, i);
    return my_data[i];
  }

  template <int Start, int Length>
  Vector<Length, SVBase<Length,1> >
  slice(){
    Internal::CheckSlice<Size, Start, Length>::check();
    return Vector<Length, SVBase<Length,1> >(&(my_data[Start]));
  }

  Vector<-1, SDVBase>
  slice(int start, int length);

private:
  double my_data[Size];
};

template<int Size>
class VBase<Size,1>{
public:
  VBase()
    : my_data(new double[Size]){
  }

  VBase(const VBase& from)
    : my_data(new double[Size]){
    }

  VBase(int size_in)
    : my_data(new double[Size]){
  }

  // construction from 1-ary operator
  template <class T, class Op>
  inline VBase(const T&, const Operator<Op>&):
    my_data(new double[Size]){}

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline VBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&):
    my_data(new double[Size]){}

  // constructor from arbitrary vector
  template<int Size2, class Base2>
  inline VBase(const Vector<Size2,Base2>& from):
    my_data(new double[from.size()]) {}

  ~VBase(){
    delete[] my_data;
  }

  double* data(){return my_data;}
  const double* data() const {return my_data;}
  
  static int size(){return Size;}
  static int stride(){return 1;}

  double& operator[](int i){
    Internal::check_index(Size, i);
    return my_data[i];
  }
  const double& operator[](int i) const {
    Internal::check_index(Size, i);
    return my_data[i];
  }

  template <int Start, int Length>
  Vector<Length, SVBase<Length,1> >
  slice(){
    Internal::CheckSlice<Size, Start, Length>::check();
    return Vector<Length, SVBase<Length,1> >(&(my_data[Start]));
  }

  Vector<-1, SDVBase>
  slice(int start, int length);

private:
  double* const my_data;
};


// SVBase does not own its data
// and has a template stride parameter
template<int Size, int Stride>
class SVBase{
public:
  SVBase(double* data_in):
    my_data(data_in){
  }

  SVBase(const SVBase& from):
    my_data(from.my_data){
  }

  double* data(){return my_data;}
  const double* data() const {return my_data;}
  
  int size() const {return Size;}
  int stride() const {return Stride;}

  double& operator[](int i){
    return my_data[i*Stride];
  }
  const double& operator[](int i) const {
    return my_data[i*Stride];
  }

  template <int Start, int Length>
  Vector<Length, SVBase<Length,Stride> >
  slice(){
    return Vector<Length, SVBase<Length,1> >(&(my_data[Start*Stride]));
  }

private:
  double* const my_data;
};




// DVBase is for vectors whose size is determined dynamically at runtime
// They own their data
class DVBase{
public:
  DVBase(int size_in):
    my_data(new double[size_in]),
    my_size(size_in){
  }

  DVBase(const DVBase& from):
    my_data(new double[from.my_size]),
    my_size(from.my_size) {
  }

  // construction from 1-ary operator
  template <class T, class Op>
  inline DVBase(const T& arg, const Operator<Op>&):
    my_data(new double[Op::size(arg)]),
    my_size(Op::size(arg)) {
  }

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline DVBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&):
    my_data(new double[Op::size(lhs,rhs)]),
    my_size(Op::size(lhs,rhs)) {
  }

  // constructor from arbitrary vector
  template<int Size2, class Base2>
  inline DVBase(const Vector<Size2,Base2>& from):
    my_data(new double[from.size()]),
    my_size(from.size()) {
  }

  ~DVBase(){
    delete[] my_data;
  }

  double* data(){return my_data;}
  const double* data() const {return my_data;}
  
  int size() const {return my_size;}
  int stride() const {return 1;}

  double& operator[](int i){
    return my_data[i];
  }
  const double& operator[](int i) const {
    return my_data[i];
  }
private:
  double* const my_data;
  int my_size;
};

// SDVBase is for dynamically sized vectors that do not own their data
// They have an additional stride member
class SDVBase{
public:
  SDVBase(double* data_in, int size_in):
    my_data(data_in) {
    my_size=size_in;
  };

  SDVBase(const SDVBase& from)
    : my_data(from.my_data),
      my_size(from.my_size){
  }

  ~SDVBase(){
  }

  double* data(){return my_data;}
  const double* data() const {return my_data;}
  
  int size() const {return my_size;}
  int stride() const {return 1;}

  double& operator[](int i){
    return my_data[i];
  }
  const double& operator[](int i) const {
    return my_data[i];
  }
private:
  double* const my_data;
  int my_size;
};

// traits classes that help with building the vectors you actually
// construct in code
static const int MAX_SIZE=10;

template<int Size>
struct VectorSelector{
  typedef VBase<Size, (Size>MAX_SIZE)?1:0 > Type;
};

template<>
struct VectorSelector<-1>{
  typedef DVBase Type;
};


////////////////////////////////////////////////////////////////////////////////
//               The actual Vector classes
////////////////////////////////////////////////////////////////////////////////



template<int Size=-1, typename Base= typename VectorSelector<Size>::Type>
class Vector : public Base {
public:
  // sneaky hack: only one of these constructors will work with any given base
  // class but they don't generate errors unless the user tries to use one of them
  // although the error message may be less than helpful - maybe this can be changed?
  inline Vector(){}
  inline Vector(double* data) : Base (data) {}
  inline Vector(int size_in) : Base(size_in) {}
  inline Vector(double* data_in, int size_in) : Base(data_in, size_in) {}
  inline Vector(double* data_in, int size_in, int stride_in) : Base(data_in, size_in, stride_in) {}


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
  inline Vector(const Vector<Size,Base>& from)
    : Base(from) {
    (*this)=from;
  }

  // constructor from arbitrary vector
  template<int Size2, class Base2>
  inline Vector(const Vector<Size2,Base2>& from):
    Base(from) {
    operator=(from);
  }

  // operator = from copy
  inline Vector<Size,Base >& operator= (const Vector& from){
    SizeMismatch<Size,Size>::test(Base::size(), from.size());
    const int s=Base::size();
    for(int i=0; i<s; i++){
      (*this)[i]=from[i];
    }
    return *this;
  }

  // operator =
  template<int Size2, typename Base2>
  Vector<Size,Base >& operator= (const Vector<Size2, Base2>& from){
    SizeMismatch<Size,Size2>::test(Base::size(), from.size());
    const int s=Base::size();
    for(int i=0; i<s; i++){
      (*this)[i]=from[i];
    }
    return *this;
  }

};



////////////////////////////////////////////////////////////////////////////////
//
// Fill in function calls, now everything is visible

template<int Size>
Vector<-1, SDVBase> VBase<Size, 0>:: slice(int start, int length){
  Internal::CheckSlice<>::check(Size, start, length);
  return Vector<-1, SDVBase>(my_data + start, length);
}

template<int Size>
Vector<-1, SDVBase> VBase<Size, 1>:: slice(int start, int length){
  Internal::CheckSlice<>::check(Size, start, length);
  return Vector<-1, SDVBase>(my_data + start, length);
}
