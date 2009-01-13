#define TOON_TEST_INTERNALS
#include <TooN/TooN.h>
#include <string>

using namespace TooN;
using namespace std;

int lineno;

template<class A, class B> struct no_duplicates
{
	typedef B type;
};


template<class A> struct no_duplicates<A, A>
{
	typedef class IgnoreMe{} type;
};


#define TRY lineno = __LINE__; try{

#define EXPECT_CATCH(X, Y) \
catch(no_duplicates<TooN::Internal::X, TooN::Internal::Y>::type e)\
{\
    cerr << "Test FAILED on line " << lineno << " from " << func_lineno << ". Expected " << #X << ", got " << #Y << "." << endl;\
}\

#define EXPECT(X) \
	if(#X == string("NoError"))\
		cerr << "Test OK on line " << lineno << " from " << func_lineno << endl;\
	else\
		cerr << "Test FAILED on line " << lineno << " from " << func_lineno << ". Expected " << #X << ", got nothing." << endl;\
}\
catch(TooN::Internal::X e)\
{\
    cerr << "Test OK on line " << lineno << " from " << func_lineno << endl;\
}\
EXPECT_CATCH(X, BadIndex)\
EXPECT_CATCH(X, SliceError)\
EXPECT_CATCH(X, StaticSliceError)\
EXPECT_CATCH(X, SizeMismatch)\
EXPECT_CATCH(X, StaticSizeMismatch)\

#define test_bad_static_slices(...) test_bad_static_slices_(__LINE__ , __VA_ARGS__)
template<class C> void test_bad_static_slices_(int func_lineno, C v)
{
	TRY{
		Vector<2> v;
		v.slice<0, 3>();
	}
	EXPECT(StaticSliceError);

	TRY{
		Vector<2> v;
		v.slice<2, 2>();
	}
	EXPECT(StaticSliceError);

	TRY{
		Vector<2> v;
		v.slice<-1, 1>();
	}
	EXPECT(StaticSliceError);
}

#define test_bad_dynamic_slices(...) test_bad_dynamic_slices_(__LINE__, __VA_ARGS__)
template<class C> void test_bad_dynamic_slices_(int func_lineno, C v)
{
	TRY{
		v.slice(0,3);
	}
	EXPECT(SliceError);

	TRY{
		v.slice(2,2);
	}
	EXPECT(SliceError);

	TRY{
		v.slice(-1,1);
	}
	EXPECT(SliceError);
}


int main()
{
	test_bad_static_slices(Vector<2>());
	test_bad_dynamic_slices(Vector<2>());
	
	test_bad_static_slices(Vector<4>().slice<0,2>());
	test_bad_dynamic_slices(Vector<4>().slice<0,2>());

	test_bad_static_slices(Vector<4>().slice(0,2));
	test_bad_dynamic_slices(Vector<4>().slice(0,2));
	
	test_bad_static_slices(Vector<>(2));
	test_bad_dynamic_slices(Vector<>(2));

	test_bad_static_slices(Vector<>(4).slice<0,2>());
	test_bad_dynamic_slices(Vector<>(4).slice<0,2>());

	test_bad_static_slices(Vector<>(4).slice(0,2));
	test_bad_dynamic_slices(Vector<>(4).slice(0,2));
}
