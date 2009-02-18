#include <TooN/TooN.h>

using namespace std;
using namespace TooN;

template<class C> void type(const C&)
{
	cout << __PRETTY_FUNCTION__ << endl;
}

int main()
{
	Vector<4> v1 = makeVector(1, 2, 3, 4);
	Vector<4> v2 = makeVector(5, 6, 7, 8);

	cout << 1+(v1 + v2)+2 << endl;

	v1.slice<0, 2>() /= 2;
	cout << v1 << endl;

	type(Vector<2>() + Vector<2, int>());
}

