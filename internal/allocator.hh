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










template<int R, int C, class Precision> struct MatrixAlloc: public StaticSizedAllocator<R*C, Precision>
{
	int num_rows() const {
		return R;
	}

	int num_cols() const {
		return C;
	}
};

template<class Precision> struct MatrixAlloc<-1, -1, Precision>
{
	const int my_rows;
	const int my_cols;
	Precision* const my_data;

	MatrixAlloc(int r, int c)
	:my_rows(r),my_cols(c),my_data(new Precision[r*c]) {
	}

	~MatrixAlloc() {
		delete[] my_data;
	}

	int num_rows() const {
		return num_rows;
	}

	int num_cols() const {
		return num_cols;
	}
};



template<int R, int C, class Precision> struct MatrixSlice
{
	int num_rows() const {
		return R;
	}

	int num_cols() const {
		return C;
	}
	//Optional Constructors
	
	Precision* const my_data;
	MatrixSlice(Precision* p)
	:my_data(p){}

	//MatrixSlice(Precision* p, int /*rows*/, int /*cols*/)
	//:my_data(p){}
};

template<class Precision> struct MatrixSlice<-1, -1, Precision>
{
	Precision* const my_data;
	const int my_rows;
	const int my_cols;

	MatrixSlice(Precision d, int r, int c)
	:my_data(d), my_rows(r),my_cols(c),my_data(d)
	{
	}

	int num_rows() const {
		return num_rows;
	}

	int num_cols() const {
		return num_cols;
	}
};











