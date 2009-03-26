//-*- c++ -*-
template<int Size=-1, typename Precision=double, typename Base=Internal::VBase>
class Vector : public Base::template Layout<Size, Precision> {
public:
  // sneaky hack: only one of these constructors will work with any given base
  // class but they don't generate errors unless the user tries to use one of them
  // although the error message may be less than helpful - maybe this can be changed?
	inline Vector(){}
	inline Vector(Precision* data) : Base::template Layout<Size, Precision> (data) {}
	inline Vector(int size_in) : Base::template Layout<Size, Precision>(size_in) {}
	inline Vector(Precision* data_in, int size_in, int stride_in, Internal::Slicing) : Base::template Layout<Size, Precision>(data_in, size_in, stride_in) {}
	inline Vector(Precision* data_in, int stride_in, Internal::Slicing) : Base::template Layout<Size, Precision>(data_in, stride_in) {}
	
	using Base::template Layout<Size, Precision>::size;

	// constructors to allow return value optimisations
	// construction from 0-ary operator
	template <class Op>
	inline Vector(const Operator<Op>&){
		Op::eval(*this);
	}

	// constructors to allow return value optimisations
	// construction from 1-ary operator
	template <class T, class Op>
	inline Vector(const T& arg, int size, const Operator<Op>&) : Base::template Layout<Size, Precision>(size) {
		Op::eval(*this,arg);
	}

	// constructor from 2-ary operator
	template <class LHS, class RHS, class Op>
	inline Vector(const LHS& lhs, const RHS& rhs, int size, const Operator<Op>&)
		: Base::template Layout<Size, Precision>(size) {
		Op::eval(*this,lhs,rhs);
	}

	// Copy construction is a very special case. Copy construction goes all the
	// way down to the bottom. GenericVBase has no idea how to copy itself.
	// However, the underlying allocator objects do.  In the case of static sized
	// objects, C++ automatically copies the data.  For slice objects, C++ copies
	// all parts (pointer and size), which is correct.  For dynamically sized
	// non-slice objects the copying has to be done by hand.
	
	// inline Vector(const Vector&from);

	// constructor from arbitrary vector
	template<int Size2, typename Precision2, typename Base2>
	inline Vector(const Vector<Size2,Precision2,Base2>& from):
		Base::template Layout<Size, Precision>(from.size()) {
		operator=(from);
	}

	// assignment from a 0-ary operator
	template <class Op>
	inline operator=(const Operator<Op>&){
		Op::eval(*this);
	}

	// operator = from copy
	inline Vector& operator= (const Vector& from){
		SizeMismatch<Size,Size>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}

	// operator =
	template<int Size2, typename Precision2, typename Base2>
	Vector<Size,Precision,Base >& operator= (const Vector<Size2, Precision2, Base2>& from){
		SizeMismatch<Size,Size2>::test(size(), from.size());
		const int s=size();
		for(int i=0; i<s; i++){
			(*this)[i]=from[i];
		}
		return *this;
	}


	Vector& operator+=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]+=rhs;
		return *this;
	}
	

	Vector& operator-=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]-=rhs;
		return *this;
	}
	
	Vector& operator/=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]/=rhs;
		return *this;
	}
	

	Vector& operator*=(const Precision& rhs) {
		for(int i=0; i<size(); i++)
			(*this)[i]*=rhs;
		return *this;
	}
	
	template<int Size2, class Precision2, class Base2>
	Vector& operator+=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]+=rhs[i];
		return *this;
	}

	template<int Size2, class Precision2, class Base2>
	Vector& operator-=(const Vector<Size2, Precision2, Base2>& rhs) {
		SizeMismatch<Size,Size2>::test(size(),rhs.size());
		for(int i=0; i<size(); i++)
			(*this)[i]-=rhs[i];
		return *this;
	}

};
