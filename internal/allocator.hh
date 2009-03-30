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
	VectorAlloc(const Operator<Op>& op) {}

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



template<int R, int C, class Precision> struct MatrixAlloc: public StaticSizedAllocator<R*C, Precision>
{
	MatrixAlloc(int,int)
	{}

	MatrixAlloc()
	{}
	int num_rows() const {
		return R;
	}

	int num_cols() const {
		return C;
	}
};


template<int Rows, class Precision> struct MatrixAlloc<Rows, -1, Precision>
{
	const int my_cols;
	Precision* const my_data;

	MatrixAlloc(const MatrixAlloc& m)
	:my_cols(m.my_cols),my_data(new Precision[Rows*my_cols]) {
		const int size=Rows*my_cols;
		for(int i=0; i < size; i++)
			my_data[i] = m.my_data[i];
	}

	MatrixAlloc(int, int c)
	:my_cols(c),my_data(new Precision[Rows*c]) {
	}

	~MatrixAlloc() {
		delete[] my_data;
	}

	int num_rows() const {
		return Rows;
	}

	int num_cols() const {
		return my_cols;
	}
};

template<int Cols, class Precision> struct MatrixAlloc<-1, Cols, Precision>
{
	const int my_rows;
	Precision* const my_data;

	MatrixAlloc(const MatrixAlloc& m)
	:my_rows(m.my_rows),my_data(new Precision[my_rows*Cols]) {
		const int size=Cols*my_rows;
		for(int i=0; i < size; i++)
			my_data[i] = m.my_data[i];
	}

	MatrixAlloc(int r, int)
	:my_rows(r),my_data(new Precision[r*Cols]) {
	}

	~MatrixAlloc() {
		delete[] my_data;
	}

	int num_rows() const {
		return my_rows;
	}

	int num_cols() const {
		return Cols;
	}
};

template<class Precision> struct MatrixAlloc<-1, -1, Precision>
{
	size_t elements(int r, int c)
	{
		if(r < 0 || c < 0)
			throw std::bad_alloc();
		return r*c;
	}
	const int my_rows;
	const int my_cols;
	Precision* const my_data;

	MatrixAlloc(const MatrixAlloc& m)
	:my_rows(m.my_rows),my_cols(m.my_cols),my_data(new Precision[elements(my_rows, my_cols)]){
		const int size=my_rows*my_cols;
		for(int i=0; i < size; i++)
			my_data[i] = m.my_data[i];
	}

	MatrixAlloc(int r, int c)
	:my_rows(r),my_cols(c),my_data(new Precision[elements(r, c)]){
	}

	~MatrixAlloc() {
		delete[] my_data;
	}

	int num_rows() const {
		return my_rows;
	}

	int num_cols() const {
		return my_cols;
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

	MatrixSlice(Precision* p, int /*rows*/, int /*cols*/)
	:my_data(p){}
};



template<int Rows, class Precision> struct MatrixSlice<Rows, -1, Precision>
{
	Precision* const my_data;
	const int my_cols;

	MatrixSlice(Precision* d, int r, int c) :my_data(d),my_cols(c) { } 
	MatrixSlice(Precision* d, int c, const Internal::SpecifyCols&) :my_data(d),my_cols(c) { } 

	int num_rows() const {
		return Rows;
	}

	int num_cols() const {
		return my_cols;
	}
};

template<int Cols, class Precision> struct MatrixSlice<-1, Cols, Precision>
{
	Precision* const my_data;
	const int my_rows;

	MatrixSlice(Precision* d, int r, int c) :my_data(d),my_rows(r) { }
	MatrixSlice(Precision* d, int r, const Internal::SpecifyRows&) :my_data(d),my_rows(r) { }

	int num_rows() const {
		return my_rows;
	}

	int num_cols() const {
		return Cols;
	}
};

template<class Precision> struct MatrixSlice<-1, -1, Precision>
{
	Precision* const my_data;
	const int my_rows;
	const int my_cols;

	MatrixSlice(Precision* d, int r, int c)
	:my_data(d), my_rows(r),my_cols(c)
	{
	}

	int num_rows() const {
		return my_rows;
	}

	int num_cols() const {
		return my_cols;
	}
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
	StrideHolder(const Operator<Op>& op) {}

	int stride() const{
		return s;
	}
};

template<> struct StrideHolder<-1>
{
	StrideHolder(int s)
	:my_stride(s){}
	
	//This constructor is not allowed to ignore the argument
	StrideHolder(int s, const Internal::NoIgnore&)
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

	RowStrideHolder(int i, const Internal::NoIgnore& n)
	:StrideHolder<S>(i, n){}

	RowStrideHolder()
	{}
};


template<int S> struct ColStrideHolder: public StrideHolder<S>
{
	ColStrideHolder(int i)
	:StrideHolder<S>(i){}

	ColStrideHolder(int i, const Internal::NoIgnore& n)
	:StrideHolder<S>(i, n){}

	ColStrideHolder()
	{}
};

}
