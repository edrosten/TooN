#define TOON_TEST_INTERNALS
#include <TooN/TooN.h>

using namespace TooN;
using namespace std;

int lineno;

#define TRY lineno = __LINE__; try{

#define EXPECT(X) \
	cerr << "Test FAILED on line " << lineno << ". Expected " << #X << ", got nothing." << endl;\
}\
catch(TooN::Internal::X e)\
{\
    cerr << "Test OK on line " << lineno << endl;\
}\
catch(...)\
{\
    cerr << "Test FAILED on line " << lineno << ". Expected " << #X << ", got unknown error." << endl;\
}\

void test_bad_static_static_slices()
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

void test_bad_static_dynamic_slices()
{
	TRY{
		Vector<2> v;
		v.slice(0,3);
	}
	EXPECT(SliceError);

	TRY{
		Vector<2> v;
		v.slice(2,2);
	}
	EXPECT(SliceError);

	TRY{
		Vector<2> v;
		v.slice(-1,1);
	}
	EXPECT(SliceError);
}

void test_bad_dynamic_static_slices()
{
	TRY{
		Vector<> v(2);
		v.slice<0, 3>();
	}
	EXPECT(SliceError);

	TRY{
		Vector<> v(2);
		v.slice<2, 2>();
	}
	EXPECT(SliceError);

	TRY{
		Vector<> v(2);
		v.slice<-1, 1>();
	}
	EXPECT(StaticSliceError);
}




int main()
{
	test_bad_static_static_slices();
	test_bad_static_dynamic_slices();
	test_bad_dynamic_static_slices();
}
