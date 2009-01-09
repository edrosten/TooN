namespace Internal
{
	template<int Size, bool StaticBad> 
	struct BadSlice;

	template<int Size> 
	struct BadSlice<Size, 0>{
		static void check(){}
	};

	template<int Size=-1, int Start=-1, int Length=-1> 
	struct CheckSlice
	{
		static void check()
		{
			BadSlice<Size, (Start+Length>=Size)>::check();
		}
	};	

	template<int Start, int Length> 
	struct CheckSlice<-1, Start, Length>
	{
		static void check(int size, int start, int length)
		{
			if(start + length >= size)
			{
				std::cerr << "Toon slice out of range" << std::endl;
				std::abort();
			}
		}
	};



}
