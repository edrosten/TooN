// -*- c++ -*-

// class to generate compile time error
// general case which doesn't exist
template<int Size1, int Size2>
struct SizeMismatch;

// special cases which do exist
template<int Size>
struct SizeMismatch<Size,Size>{
  static inline void test(int, int){}
};

template<int Size>
struct SizeMismatch<-1,Size>{
  static inline void test(int size1, int size2){
    if(size1!=size2){
	  #ifdef TOON_TEST_INTERNALS
	  	throw Internal::SizeMismatch();
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
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
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
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
	  #elif !defined TOON_NDEBUG_SIZE
		  std::cerr << "TooN Size Mismatch" << std::endl;
		  std::abort();
	  #endif
    }
  }
};

namespace Internal
{
	struct BadSize;
}

template<int Size1, int Size2>
struct SizeMismatch
{
	static inline void test(int, int)
	{
		#ifdef TOON_TEST_INTERNALS
			throw Internal::StaticSizeMismatch();
		#else
			Internal::BadSize size_mismatch;
		#endif
	}
};
