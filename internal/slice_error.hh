namespace Internal
{
	template<bool StaticBad> 
	struct BadSlice;

	template<> 
	struct BadSlice<0>{
		static void check(){}
	};

	template<int Size=-1, int Start=2147483647, int Length=-1> 
	struct CheckSlice
	{
		static void check()
		{
			BadSlice<(Start < 0) || (Start+Length>=Size)>::check();
		}
	};	

	template<int Start, int Length> 
	struct CheckSlice<-1, Start, Length>
	{
		static void check(int size, int start, int length)
		{
			BadSlice<(Start != 2147483647 && Start < 0)>::check();

			if(start < 0 || start + length >= size)
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
