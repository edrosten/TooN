#include <TooN/TooN.h>

using namespace std;
using namespace TooN;

int main()
{
	Vector<4> v1 = makeVector(1, 2, 3, 4);
	Vector<4> v2 = makeVector(5, 6, 7, 8);

	cout << v1 + v2 << endl;

	v1.slice<0, 2>() /= 2;
	cout << v1 << endl;
}

