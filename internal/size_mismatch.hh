// -*- c++ -*-

// class to generate compile time error
// general case which doesn't exist
template<int Size1, int Size2>
struct SizeMismatch;

// special cases which do exist
template<int Size>
struct SizeMismatch<Size,Size>{
  static inline void test(int size1, int size2){}
};

template<int Size>
struct SizeMismatch<-1,Size>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #else
		  std::cerr << "Toon Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

template<int Size>
struct SizeMismatch<Size,-1>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #else
		  std::cerr << "Toon Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

template <>
struct SizeMismatch<-1,-1>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #else
		  std::cerr << "Toon Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};


#ifdef TOON_TEST_INTERNALS
  template<int Size1, int Size2>
  struct SizeMismatch
  {
    static inline void test(int, int)
    {
      throw Internal::StaticSizeMismatch();
    }
  };
#endif
