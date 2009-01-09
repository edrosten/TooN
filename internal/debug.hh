namespace Internal
{
	static inline void check_index(int i, int s)
	{
		#ifdef TOON_CHECK_BOUNDS
			if(i<0 || i >= s)
			{
				std::cerr << "Toon index out of range" << std::endl;
				std::abort();
			}
		#endif
	}
}
