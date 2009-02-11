namespace Internal
{
	template<bool StaticBad> 
	struct BadSlice;

	template<> 
	struct BadSlice<0>{
		static void check(){}
	};

	template<int Size, int Start, int Length> 
	struct CheckStaticSlice
	{
		static void check(int/*size*/)
		{
			BadSlice<(Start < 0) || (Length < 1) || (Start+Length>Size)>::check();
		}
	};	
	
	template<int Start, int Length> struct CheckStaticSlice<-1, Start, Length>{
		static void check(int size)
		{
			BadSlice<(Start < 0) || (Length < 1)>::check();
			if(Start + Length > size)
			{
				#ifdef TOON_TEST_INTERNALS
					throw Internal::SliceError();
				#else
					std::cerr << "Toon slice out of range (static slice, synamic vector)" << std::endl;
					std::abort();
				#endif
			}
		}
	};

	struct CheckDynamicSlice{
		static void check(int size, int start, int length){
			if(start < 0 || start + length > size)
			{
				#ifdef TOON_TEST_INTERNALS
					throw Internal::SliceError();
				#else
					std::cerr << "Toon slice out of range" << std::endl;
					std::abort();
				#endif
			}
		}
	};

	#ifdef TOON_TEST_INTERNALS
		template<bool StaticBad> 
		struct BadSlice{
			static void check(){
				throw Internal::StaticSliceError();
			}
		};
	#endif



}
