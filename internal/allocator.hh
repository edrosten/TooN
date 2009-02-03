// Allocators never copy

template<int Size, class Precision, bool heap> class StackOrHeap
{
	public:
		StackOrHeap()
		{}

		StackOrHeap(const StackOrHeap&)
		{}

		Precision my_data[Size];
};

template<int Size, class Precision> class StackOrHeap<Size, Precision, 1>
{
	public:
		StackOrHeap()
		:my_data(new Precision[Size]){}


		~StackOrHeap()
		{
			delete[] my_data;
		}

		Precision *my_data;
	
	private:
		StackOrHeap(const StackOrHeap&);
};


template<class Precision> class DynamicSizedAllocator
{
	public:
		DynamicSizedAllocator(int size)
		:my_data(new Precision[size]){}


		~DynamicSizedAllocator()
		{
			delete[] my_data;
		}

		Precision *my_data;

	private:
		DynamicSizedAllocator(const DynamicSizedAllocator& d);
};



template<int Size, class Precision> class StaticSizedAllocator: public StackOrHeap<Size, Precision, (sizeof(Precision)*Size>max_bytes_on_stack) >
{
};


template<class Precision> struct SliceHolder
{
	SliceHolder(Precision* p)
	:my_data(p)
	{}

	Precision* my_data;
};


