//-*- c++ -*-

// All the vector base classes go in this file
// they have to provide the following capabilities:
//
// template<typename Precision> // and possibly more stuff
// class VectorBase {
//   // Constructors are only responsible for memory management
//   // not for copying any data
//   // in general the following constructors should exist:

//   VectorBase(...) // basic constructor
//   VectorBase(const VectorBase& from) // copy constructor
  
//   // construction from 1-ary operator
//   template <class T, class Op>
//   inline VectorBase(const T&, const Operator<Op>&);

//   // constructor from 2-ary operator
//   template <class LHS, class RHS, class Op>
//   inline VectorBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&);
  
//   // constructor from arbitrary vector
//   template<int Size2, class Base2>
//   inline VectorBase(const Vector<Size2, Precision, Base2>& from);

//   Precision* data();
//   const Precision* data() const;
  
//   static int size(); // returns the size of the vector
//   static int stride(); // returns the stride

//   Precision& operator[](int i); // return element i
//   const Precision& operator[](int i) const; // ditto

//   template <int Start, int Length>
//   Vector<Length,Precision,...> slice(); // static slice
//   Vector<-1, Precision, ...>  slice(int start, int length); // dynamic slice
// };


// forward declarations
template<int Size, typename Precision, typename Base>
class Vector;

template<int Size, int Stride, typename Precision>
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
template <int Size, int Type, typename Precision>
class VBase;

template<int Stride, typename Precision> 
class SDVBase; //Sliced Dynamic VBase. No ownership of data, Static stride.

template<int Size, typename Precision>
class VBase<Size, 0, Precision> {
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
  inline VBase(const Vector<Size2, Precision, Base2>&){}

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  static int size(){return Size;}
  static int stride(){return 1;}

  Precision& operator[](int i){
    Internal::check_index(Size, i);
    return my_data[i];
  }
  const Precision& operator[](int i) const {
    Internal::check_index(Size, i);
    return my_data[i];
  }

  template <int Start, int Length>
  Vector<Length, Precision, SVBase<Length,1,Precision> >
  slice(){
    Internal::CheckSlice<Size, Start, Length>::check();
    return Vector<Length, Precision, SVBase<Length,1,Precision> >(&(my_data[Start]));
  }

  Vector<-1, Precision, SDVBase<1, Precision> >
  slice(int start, int length);

private:
  Precision my_data[Size];
};

template<int Size, typename Precision>
class VBase<Size,1, Precision>{
public:

  // Constructors
  // For all these, the owning vector class is responsible
  // for filling in any data that needs to go into the vector elements

  VBase()
    : my_data(new Precision[Size]){
  }

  VBase(const VBase& from)
    : my_data(new Precision[Size]){
    }

  VBase(int size_in)
    : my_data(new Precision[Size]){
  }
  
  // construction from 1-ary operator
  template <class T, class Op>
  inline VBase(const T&, const Operator<Op>&):
    my_data(new Precision[Size]){}

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline VBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&):
    my_data(new Precision[Size]){}

  // constructor from arbitrary vector
  template<int Size2, typename Precision2, typename Base2>
  inline VBase(const Vector<Size2,Precision2,Base2>& from):
    my_data(new Precision[from.size()]) {}

  ~VBase(){
    delete[] my_data;
  }

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  static int size(){return Size;}
  static int stride(){return 1;}

  Precision& operator[](int i){
    Internal::check_index(Size, i);
    return my_data[i];
  }
  const Precision& operator[](int i) const {
    Internal::check_index(Size, i);
    return my_data[i];
  }

  template <int Start, int Length>
  Vector<Length, Precision, SVBase<Length,1, Precision> >
  slice(){
    Internal::CheckSlice<Size, Start, Length>::check();
    return Vector<Length, Precision, SVBase<Length,1,Precision> >(&(my_data[Start]));
  }

  Vector<-1, Precision, SDVBase<1, Precision> >
  slice(int start, int length);

private:
  Precision* const my_data;
};


// SVBase does not own its data
// and has a template stride parameter
template<int Size, int Stride, typename Precision>
class SVBase{
public:
  SVBase(Precision* data_in):
    my_data(data_in){
  }

  SVBase(const SVBase& from):
    my_data(from.my_data){
  }

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  int size() const {return Size;}
  int stride() const {return Stride;}

  Precision& operator[](int i){
    Internal::check_index(Size, i);
    return my_data[i*Stride];
  }
  const Precision& operator[](int i) const {
    Internal::check_index(Size, i);
    return my_data[i*Stride];
  }

  template <int Start, int Length>
  Vector<Length, Precision, SVBase<Length,Stride,Precision> >
  slice(){
    Internal::CheckSlice<Size, Start, Length>::check();
    return Vector<Length, Precision, SVBase<Length,1,Precision> >(&(my_data[Start*Stride]));
  }

  Vector<-1, Precision, SDVBase<Stride, Precision> >
  slice(int start, int length);

private:
  Precision* const my_data;
};




// DVBase is for vectors whose size is determined dynamically at runtime
// They own their data
template <typename Precision>
class DVBase{
public:
  DVBase(int size_in):
    my_data(new Precision[size_in]),
    my_size(size_in){
  }

  DVBase(const DVBase& from):
    my_data(new Precision[from.my_size]),
    my_size(from.my_size) {
  }

  // construction from 1-ary operator
  template <class T, class Op>
  inline DVBase(const T& arg, const Operator<Op>&):
    my_data(new Precision[Op::size(arg)]),
    my_size(Op::size(arg)) {
  }

  // constructor from 2-ary operator
  template <class LHS, class RHS, class Op>
  inline DVBase(const LHS& lhs, const RHS& rhs, const Operator<Op>&):
    my_data(new Precision[Op::size(lhs,rhs)]),
    my_size(Op::size(lhs,rhs)) {
  }

  // constructor from arbitrary vector
  template<int Size2, class Base2>
  inline DVBase(const Vector<Size2, Precision, Base2>& from):
    my_data(new Precision[from.size()]),
    my_size(from.size()) {
  }

  ~DVBase(){
    delete[] my_data;
  }

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  int size() const {return my_size;}
  int stride() const {return 1;}

  Precision& operator[](int i){
    Internal::check_index(my_size, i);
    return my_data[i];
  }
  const Precision& operator[](int i) const {
    Internal::check_index(my_size, i);
    return my_data[i];
  }

  template <int Start, int Length>
  Vector<Length, Precision, SVBase<Length,1,Precision> >
  slice(){
    Internal::CheckSlice<-1, Start>::check(my_size, Start, Length);
    return Vector<Length, Precision, SVBase<Length,1,Precision> >(&(my_data[Start]));
  }

  Vector<-1, Precision, SDVBase<1, Precision> >
  slice(int start, int length);

private:
  Precision* const my_data;
  int my_size;
};

// SDVBase is for dynamically sized vectors that do not own their data
// They have an additional templated stride
template <int Stride, typename Precision>
class SDVBase{
public:
  SDVBase(Precision* data_in, int size_in):
    my_data(data_in) {
    my_size=size_in;
  };

  SDVBase(const SDVBase& from)
    : my_data(from.my_data),
      my_size(from.my_size){
  }

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  int size() const {return my_size;}
  int stride() const {return Stride;}

  Precision& operator[](int i){
    Internal::check_index(my_size, i);
    return my_data[i*Stride];
  }
  const Precision& operator[](int i) const {
    Internal::check_index(my_size, i);
    return my_data[i*Stride];
  }

  template <int Start, int Length>
  Vector<Length, Precision, SVBase<Length,Stride,Precision> >
  slice(){
    Internal::CheckSlice<-1, Start>::check(my_size, Start, Length);
    return Vector<Length, Precision, SVBase<Length,Stride,Precision> >(&(my_data[Start]));
  }

  Vector<-1, Precision, SDVBase<Stride, Precision> >
  slice(int start, int length);


private:
  Precision* const my_data;
  int my_size;
};

// SDVBase is for dynamically sized vectors that do not own their data
// They have an additional stride member
template <typename Precision>
class SDDVBase{
public:
  SDDVBase(Precision* data_in, int size_in, int stride_in):
    my_data(data_in) {
    my_size=size_in;
    my_stride=stride_in;
  };

  SDDVBase(const SDDVBase& from)
    : my_data(from.my_data),
      my_size(from.my_size),
      my_stride(from.my_stride){
  }

  Precision* data(){return my_data;}
  const Precision* data() const {return my_data;}
  
  int size() const {return my_size;}
  int stride() const {return my_stride;}

  Precision& operator[](int i){
    Internal::check_index(my_size, i);
    return my_data[i*my_stride];
  }

  const Precision& operator[](int i) const {
    Internal::check_index(my_size, i);
    return my_data[i*my_stride];
  }
private:
  Precision* const my_data;
  int my_size;
  int my_stride;
};

// traits classes that help with building the vectors you actually
// construct in code
static const int MAX_SIZE=10;

template<int Size, typename Precision>
struct VectorSelector{
  typedef VBase<Size, (Size>MAX_SIZE)?1:0, Precision > Type;
};

template<typename Precision>
struct VectorSelector<-1, Precision>{
  typedef DVBase<Precision> Type;
};

