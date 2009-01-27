
#ifndef TOON_MAX_STACK_SIZE
static const int TOON_MAX_STACK_SIZE=10;
#endif

template <int Rows, int Cols>
struct MMemSize{
  static const int size=Rows*Cols;
};

template<>
struct MMemSize<-1,-1> {
  static const int size=-1;
};


template<int Size, typename Precision, int Place=(Size>TOON_MAX_STACK_SIZE?1:0)>
struct MemBase;

// stack memory
template<int Size, typename Precision>
struct MemBase<Size, Precision, 0> {



  Precision my_data[Size];
};
  
// heap memory
template<int Size, typename Precision>
struct MemBase<Size, Precision, 1>{
  MemBase() : my_data(new Precision[Size]){
  }
  ~MemBase(){
    delete[] my_data;
  }



  Precision* const my_data;
};



// dynamic memory
template<typename Precision>
struct MemBase<-1,Precision,0> {
  MemBase(const int Size) : my_data(new Precision[Size]), my_size(size){
  }
  ~MemBase(){
    delete[] my_data;
  }



  Precision* const my_data;
  int my_size;
};
  
