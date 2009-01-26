template <int Rows=-1, int Cols=Rows, class Precision=double, template<int R, int C, class P> class Layout = RowMajor>
class Matrix : public Layout<Rows, Cols, Precision>
{
	public:
		//Use Tom's sneaky constructor hack...
		
		Matrix(){}

		Matrix(Precision* data, Slicing)
		:Layout<Rows, Cols, Precision>(data){}







};
