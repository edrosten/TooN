namespace TooN{
namespace Internal{

	template<>
	struct HaveExprTemplates<ExprTemplatePresenceDummy>
	{
		static const bool have=1;
	};


template<class B, class P> struct BaseUtilities;

template<class Expr> struct SliceExpr
{
	template<int I, class P> struct VLayout: public BaseUtilities<SliceExpr<Expr>, P>
	{
		typedef P Precision;
		const int start, length;
		const Vector<Dynamic, Precision, Expr>& vec;
		
		int size() const
		{
			return length;
		}

		Precision index(int i) const
		{
			return vec[start + i];
		}
		
		VLayout(int start_, int length_, const Vector<Dynamic, Precision, Expr>& vec_)
		:start(start_),length(length_),vec(vec_)
		{
		}

		static const int NumberOfMuls =     Vector<Dynamic, Precision, Expr>::NumberOfMuls;
		static const int NumberOfAdds =     Vector<Dynamic, Precision, Expr>::NumberOfAdds;
		static const int NumberOfDivs =     Vector<Dynamic, Precision, Expr>::NumberOfDivs;
	};
};

template<class Base, class Precision = typename Base::Precision> struct BaseUtilities
{
	//This class is used only as part of the curiously recurring template pattern.
	//This is the class that BaseUtilities muse have been derived from.
	typedef Vector<Dynamic,  Precision, Base> Derived;

	Vector<Dynamic, Precision, SliceExpr<Base> > slice(int start, int length)
	{
		const Derived& derived = static_cast<const Derived&>(*this);

		return typename SliceExpr<Base>::template VLayout<Dynamic, Precision>(start, length, derived);
	}
	
	//Magic numbers to trade off various operations
	int cost() const
	{
		const Derived& derived = static_cast<const Derived&>(*this);
		return (Derived::NumberOfMuls * 2 + Derived::NumberOfAdds + 30 * Derived::NumberOfDivs)*derived.size();
	}

	bool should_cache(int passes) const
	{
		return (cost() * passes > 1000);
	}


	//Things that Vector<> expects. They aren't meaningful, but cause a compile
	//error if they are missing.
	typedef void* PointerType;
	typedef void  try_destructive_resize;


	//The class can optionally cache values. It makes sense to do this
	//when there are multiple passes of the vector, and the cost of allocating memory
	//is less than the cost of recmoputing an complex expression many times.
	mutable Precision* cache;

	BaseUtilities()
	:cache(0)
	{}

	~BaseUtilities()
	{
		if(cache)
			delete[] cache;
	}

	//This is an advisory function which will may cause an expression template
	//to precompute and cache the values.
	void perhaps_cache(int passes) const
	{
		const Derived& derived = static_cast<const Derived&>(*this);
		if(!cache && should_cache(passes))
		{
			cache = new Precision[derived.size()];
			for(int i=0; i < derived.size(); i++)
				cache[i] = derived.index[i];
		}
	}

	Precision operator[](int i) const
	{
		const Derived& derived = static_cast<const Derived&>(*this);
		if(cache)
			return cache[i];
		else
			return derived.index(i);
	}
	
};

template<class P1, class P2, class B2> struct ScalarMulExpr
{
	typedef typename Internal::MultiplyType<P1,P2>::type Precision;

	template<int I, class P> struct VLayout: public BaseUtilities<ScalarMulExpr<P1, P2, B2> >
	{
		typedef typename Internal::MultiplyType<P1,P2>::type Precision;
		const P1& mul;
		const Vector<Dynamic, P2, B2>& vec;
		
		int size() const
		{
			return vec.size();
		}

		Precision index(int i) const
		{
			return vec[i]*mul;
		}
		
		VLayout(const P1& m, const Vector<Dynamic, P2, B2>& vec_)
		:mul(m),vec(vec_)
		{
		}

		static const int NumberOfMuls = 1 + Vector<Dynamic, P2, B2>::NumberOfMuls;
		static const int NumberOfAdds =     Vector<Dynamic, P2, B2>::NumberOfAdds;
		static const int NumberOfDivs =     Vector<Dynamic, P2, B2>::NumberOfDivs;
	};
};

template<class P1, class P2, class B2> struct ScalarDivExpr
{
	typedef typename Internal::MultiplyType<P1,P2>::type Precision;

	template<int I, class P> struct VLayout: public BaseUtilities<ScalarDivExpr<P1, P2, B2> >
	{
		typedef typename Internal::MultiplyType<P1,P2>::type Precision;
		const P1& div;
		const Vector<Dynamic, P2, B2>& vec;
		
		int size() const
		{
			return vec.size();
		}

		Precision index(int i) const
		{
			return vec[i]/div;
		}
		
		VLayout(const P1& d, const Vector<Dynamic, P2, B2>& vec_)
		:div(d),vec(vec_)
		{
		}

		static const int NumberOfMuls =     Vector<Dynamic, P2, B2>::NumberOfMuls;
		static const int NumberOfAdds =     Vector<Dynamic, P2, B2>::NumberOfAdds;
		static const int NumberOfDivs = 1+  Vector<Dynamic, P2, B2>::NumberOfDivs;
	};
};

template<class P1, class B1> struct NegExpr
{
	typedef P1 Precision;

	template<int I, class P> struct VLayout: public BaseUtilities<NegExpr<P1, B1> >
	{
		typedef P1 Precision;
		const Vector<Dynamic, P1, B1>& vec;
		
		int size() const
		{
			return vec.size();
		}

		Precision index(int i) const
		{
			return -vec[i];
		}
		
		VLayout( const Vector<Dynamic, P1, B1>& vec_)
		:vec(vec_)
		{
		}

		static const int NumberOfMuls =     Vector<Dynamic, P1, B1>::NumberOfMuls;
		static const int NumberOfAdds =     Vector<Dynamic, P1, B1>::NumberOfAdds;
		static const int NumberOfDivs =     Vector<Dynamic, P1, B1>::NumberOfDivs;
	};
};


template<class P1, class P2, class B1, class B2> 
struct AddExpr
{
	typedef typename Internal::AddType<P1,P2>::type Precision;
	typedef Vector<Dynamic, P1, B1> LHS;
	typedef Vector<Dynamic, P2, B2> RHS;

	template<int I,class P> struct VLayout: public BaseUtilities<AddExpr<P1,P2,B1,B2> >
	{
		typedef AddExpr<P1,P2,B1,B2>::Precision Precision;
		const LHS& lhs;
		const RHS& rhs;
		
		Precision index(int i) const
		{
			return lhs[i] + rhs[i];
		}

		int size() const
		{
			return lhs.size();
		}

		VLayout(const LHS& lhs_, const RHS& rhs_)
		:lhs(lhs_),rhs(rhs_)
		{
			SizeMismatch<Dynamic,Dynamic>:: test(lhs.size(),rhs.size());
		}

		static const int NumberOfMuls =     Vector<Dynamic, P1, B1>::NumberOfMuls + Vector<Dynamic, P2, B2>::NumberOfMuls;
		static const int NumberOfAdds = 1+  Vector<Dynamic, P1, B1>::NumberOfAdds + Vector<Dynamic, P2, B2>::NumberOfAdds;
		static const int NumberOfDivs =     Vector<Dynamic, P1, B1>::NumberOfDivs + Vector<Dynamic, P2, B2>::NumberOfDivs;
	};
};


template<class P1, class P2, class B1, class B2> 
struct SubExpr
{
	typedef typename Internal::SubtractType<P1,P2>::type Precision;
	typedef Vector<Dynamic, P1, B1> LHS;
	typedef Vector<Dynamic, P2, B2> RHS;
	template<int I,class P> struct VLayout: public BaseUtilities<SubExpr<P1,P2,B1,B2> >
	{
		typedef SubExpr<P1,P2,B1,B2>::Precision Precision;
		const LHS& lhs;
		const RHS& rhs;
		
		Precision index(int i) const
		{
			return lhs[i] - rhs[i];
		}

		int size() const
		{
			return lhs.size();
		}

		VLayout(const LHS& lhs_, const RHS& rhs_)
		:lhs(lhs_),rhs(rhs_)
		{
			SizeMismatch<Dynamic,Dynamic>:: test(lhs.size(),rhs.size());
		}


		static const int NumberOfMuls =     Vector<Dynamic, P1, B1>::NumberOfMuls + Vector<Dynamic, P2, P2>::NumberOfMuls;
		static const int NumberOfAdds = 1+  Vector<Dynamic, P1, B1>::NumberOfAdds + Vector<Dynamic, P2, B2>::NumberOfAdds;
		static const int NumberOfDivs =     Vector<Dynamic, P1, B1>::NumberOfDivs + Vector<Dynamic, P2, B2>::NumberOfDivs;
	};
};
}


template<class P1, class P2, class B1, class B2>
Vector<Dynamic, typename Internal::SubExpr<P1,P2,B1,B2>::Precision, Internal::SubExpr<P1,P2,B1,B2> > operator-(const Vector<Dynamic,P1,B1>& lhs, const Vector<Dynamic, P2, B2>& rhs)
{
	typedef Internal::SubExpr<P1,P2,B1,B2> Expr;
	typedef typename Expr::Precision Precision;
	return typename Expr::template VLayout<Dynamic, Precision>(lhs, rhs);
}


template<class P1, class P2, class B1, class B2>
Vector<Dynamic, typename Internal::AddExpr<P1,P2,B1,B2>::Precision, Internal::AddExpr<P1,P2,B1,B2> > operator+(const Vector<Dynamic,P1,B1>& lhs, const Vector<Dynamic, P2, B2>& rhs)
{
	typedef Internal::AddExpr<P1,P2,B1,B2> Expr;
	typedef typename Expr::Precision Precision;
	typedef typename Expr::template VLayout<Dynamic, Precision> Base;
	return Base(lhs, rhs);
}

template<class P1, class P2, class B2>
typename Internal::enable_if<
	Internal::NotVector<P1>,
    Vector<Dynamic, typename Internal::ScalarMulExpr<P1,P2,B2>::Precision, Internal::ScalarMulExpr<P1,P2,B2> >
>::type operator*(const P1& lhs, const Vector<Dynamic, P2, B2>& rhs)
{
	typedef Internal::ScalarMulExpr<P1,P2,B2> Expr;
	typedef typename Expr::Precision Precision;
	typedef typename Expr::template VLayout<Dynamic, Precision> Base;
	return Base(lhs, rhs);
}

template<class P1, class P2, class B2>
typename Internal::enable_if<
	Internal::NotVector<P1>,
	Vector<Dynamic, typename Internal::ScalarMulExpr<P1,P2,B2>::Precision, Internal::ScalarMulExpr<P1,P2,B2> > 
>::type operator*(const Vector<Dynamic, P2, B2>& lhs, const P1& rhs)
{
	typedef Internal::ScalarMulExpr<P1,P2,B2> Expr;
	typedef typename Expr::Precision Precision;
	typedef typename Expr::template VLayout<Dynamic, Precision> Base;
	return Base(rhs, lhs);
}

template<class P1, class P2, class B2>
Vector<Dynamic, typename Internal::ScalarDivExpr<P1,P2,B2>::Precision, Internal::ScalarDivExpr<P1,P2,B2> > operator/(const Vector<Dynamic, P2, B2>& lhs, const P1& rhs)
{
	typedef Internal::ScalarDivExpr<P1,P2,B2> Expr;
	typedef typename Expr::Precision Precision;
	typedef typename Expr::template VLayout<Dynamic, Precision> Base;
	return Base(rhs, lhs);
}

template<class P1, class B1>
Vector<Dynamic, P1, Internal::NegExpr<P1,B1> > operator-(const Vector<Dynamic, P1, B1>& v)
{
	return typename Internal::NegExpr<P1,B1>::template VLayout<Dynamic, P1>(v);
}

}

