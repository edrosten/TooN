
#ifndef TOON_MAX_STACK_SIZE
static const int TOON_MAX_STACK_SIZE=10;
#endif

template <int Rows, int Cols, typename Precision, int Place=(Rows*Cols>TOON_MAX_STACK_SIZE?1:0)>
struct MMemBase;

// Place=0 => stack
template <int Rows, int Cols, typename Precision>
struct MMemBase<Rows, Cols, Precision, 0> {
  Precision my_data[Rows*Cols];
};


// Place=1 => heap
template <int Rows, int Cols, typename Precision>
struct MMemBase <Rows, Cols, Precision, 1>{
  MMemBase() : my_data(new Precision[Rows*Cols]) {}
  ~MMemBase() {delete[] my_data}
  Precision* const my_data;
};


// Rows=Cols=-1 => dynamic
tempate<typename Precision>
struct MMembase<-1,-1,Precision,0> {
  MMemBase(unsigned int r, unsigned int c) : my_data(new Precision[r*c]), my_rows(r), my_cols(c){}
  ~MMemBase() {delete[] my_data;}
  Precision* my_data;
  unsigned int my_rows;
  unsigned int my_cols;
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
  
