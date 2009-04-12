// Allocators always copy data on copy construction.
//
// When a Vector/Matrix is constructed from a different, but compatible type
// copying is done at a much higher level: above the level that knows how the
// data is laid out in memory.
//
// At this level, copy construction is required since it is only known here 
// whether data or a reference to data should be copied.
namespace Internal
{

template<int Size, class Precision, bool heap> class StackOrHeap
{
	public:
		StackOrHeap()
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
	
		StackOrHeap(const StackOrHeap& from)
		:my_data(new Precision[Size])
		{
			for(int i=0; i < Size; i++)
				my_data[i] = from.my_data[i];
		}
};


template<int Size, class Precision> class StaticSizedAllocator: public StackOrHeap<Size, Precision, (sizeof(Precision)*Size>max_bytes_on_stack) >
{
};

template<int Size, class Precision> struct VectorAlloc : public StaticSizedAllocator<Size, Precision> {
	
	VectorAlloc() { }
	
	VectorAlloc(int /*s*/) { }

	template<class Op>
	VectorAlloc(const Operator<Op>&) {}

	int size() const {
		return Size;
	}
};

template<class Precision> struct VectorAlloc<-1, Precision> {
	Precision * const my_data;
	const int my_size;

	VectorAlloc(const VectorAlloc& v)
	:my_data(new Precision[v.my_size]), my_size(v.my_size)
	{ 
		for(int i=0; i < my_size; i++)
			my_data[i] = v.my_data[i];
	}

	VectorAlloc(int s)
	:my_data(new Precision[s]), my_size(s)
	{ }

	template <class Op>
	VectorAlloc(const Operator<Op>& op) : my_data(new Precision[op.size()]), my_size(op.size()) {}

	int size() const {
		return my_size;
	}

	~VectorAlloc(){
		delete[] my_data;
	}

};


template<int S, class Precision> struct VectorSlice
{
	int size() const {
		return S;
	}

	//Optional Constructors
	
	Precision* const my_data;
	VectorSlice(Precision* p)
	:my_data(p){}

	VectorSlice(Precision* p, int /*size*/)
	:my_data(p){}

	template<class Op>
	VectorSlice(const Operator<Op>& op) : my_data(op.data()) {}
};

template<class Precision> struct VectorSlice<-1, Precision>
{
	Precision* const my_data;
	const int my_size;

	VectorSlice(Precision* d, int s)
	:my_data(d), my_size(s)
	{ }

	template<class Op>
	VectorSlice(const Operator<Op>& op) : my_data(op.data()), my_size(op.size()) {}

	int size() const {
		return my_size;
	}
};




////////////////////////////////////////////////////////////////////////////////
//
// A class similar to StrideHolder, but to hold the size information for matrices.

template<int s> struct SizeHolder
{
	//Constructors ignore superfluous arguments
	SizeHolder(){}
	SizeHolder(int){}

	int size() const{
		return s;
	}
};

template<> struct SizeHolder<-1>
{
	SizeHolder(int s)
	:my_size(s){}

	const int my_size;
	int size() const {
		return my_size;
	}
};


template<int S> struct RowSizeHolder: private SizeHolder<S>
{
	RowSizeHolder(int i)
	:SizeHolder<S>(i){}

	RowSizeHolder()
	{}

	template<typename Op>
	RowSizeHolder(const Operator<Op>& op) : SizeHolder<S>(op.num_rows()) {}

	int num_rows() const {return SizeHolder<S>::size();}
};


template<int S> struct ColSizeHolder: public SizeHolder<S>
{
	ColSizeHolder(int i)
	:SizeHolder<S>(i){}

	ColSizeHolder()
	{}

	template<typename Op>
	ColSizeHolder(const Operator<Op>& op) : SizeHolder<S>(op.num_cols()) {}

	int num_cols() const {return SizeHolder<S>::size();}	
};



template<int R, int C, class Precision, bool FullyStatic=(R>=0 && C>=0)> 
struct MatrixAlloc: public StaticSizedAllocator<R*C, Precision>
{
	MatrixAlloc(int,int)
	{}

	MatrixAlloc()
	{}

	template <class Op>
	MatrixAlloc(const Operator<Op>&)
	{}

	int num_rows() const {
		return R;
	}

	int num_cols() const {
		return C;
	}
};


template<int R, int C, class Precision> struct MatrixAlloc<R, C, Precision, false>
	: public RowSizeHolder<R>,
	ColSizeHolder<C>
{
	Precision* const my_data;

	using RowSizeHolder<R>::num_rows;
	using ColSizeHolder<C>::num_cols;
	
	// copy constructor so guaranteed contiguous
	MatrixAlloc(const MatrixAlloc& m)
		:RowSizeHolder<R>(m.num_rows()),
		 ColSizeHolder<C>(m.num_cols()),
		 my_data(new Precision[num_rows()*num_cols()]) {
		const int size=num_rows()*num_cols();
		for(int i=0; i < size; i++) {
			my_data[i] = m.my_data[i];
		}
	}

	MatrixAlloc(int r, int c)
		:RowSizeHolder<R>(r),
		 ColSizeHolder<C>(c),
		 my_data(new Precision[num_rows()*num_cols()]) {
	}

	template <class Op>	MatrixAlloc(const Operator<Op>& op)
		:RowSizeHolder<R>(op),
		 ColSizeHolder<C>(op),
		 my_data(new Precision[num_rows()*num_cols()])
	{}

	~MatrixAlloc() {
		delete[] my_data;
	}
};


template<int R, int C, class Precision> struct MatrixSlice
	: public RowSizeHolder<R>,
	ColSizeHolder<C>
{
	Precision* const my_data;

	using RowSizeHolder<R>::num_rows;
	using ColSizeHolder<C>::num_cols;

	//Optional Constructors
	MatrixSlice(Precision* p)
	:my_data(p){}

	MatrixSlice(Precision* p, int r, int c)
		:RowSizeHolder<R>(r),
		 ColSizeHolder<C>(c),
		 my_data(p){}
	
	template<class Op>
	MatrixSlice(const Operator<Op>& op)
		:RowSizeHolder<R>(op),
		 ColSizeHolder<C>(op),
		 my_data(op.data())
	{}
};


////////////////////////////////////////////////////////////////////////////////
//
// A class similar to mem, but to hold the stride information. It is only needed
// for -1. For +int and -2, the stride is part fo teh type, or implicit.

template<int s> struct StrideHolder
{
	//Constructos ignore superfluous arguments
	StrideHolder(){}
	StrideHolder(int){}

	template<class Op>
	StrideHolder(const Operator<Op>&) {}

	int stride() const{
		return s;
	}
};

template<> struct StrideHolder<-1>
{
	StrideHolder(int s)
	:my_stride(s){}

	template<class Op>
	StrideHolder(const Operator<Op>& op) : my_stride(op.stride()) {}

	const int my_stride;
	int stride() const {
		return my_stride;
	}
};


template<int S> struct RowStrideHolder: public StrideHolder<S>
{
	RowStrideHolder(int i)
	:StrideHolder<S>(i){}

	RowStrideHolder()
	{}

	template<class Op>
	RowStrideHolder(const Operator<Op>& op)
		: StrideHolder<S>(op)
	{}

};


template<int S> struct ColStrideHolder: public StrideHolder<S>
{
	ColStrideHolder(int i)
	:StrideHolder<S>(i){}

	ColStrideHolder()
	{}

	template<class Op>
	ColStrideHolder(const Operator<Op>& op)
		: StrideHolder<S>(op)
	{}
};





}
