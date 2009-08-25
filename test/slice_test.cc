#include <TooN/TooN.h>
using namespace std;
using namespace TooN;

int main()
{
	Vector<6> v = makeVector(1, 2, 3, 4, 5, 6);
	Vector<>  u = v;

	cout << v.slice<0,3>() << endl;
	cout << u.slice<0,3>() << endl << endl;

	cout << v.slice<1,3>(1, 3) << endl;
	cout << u.slice<1,3>(1, 3) << endl;
	cout << v.slice<Dynamic,3>(1, 3) << endl;
	cout << u.slice<Dynamic,3>(1, 3) << endl;
	cout << v.slice<1,Dynamic>(1, 3) << endl;
	cout << u.slice<1,Dynamic>(1, 3) << endl;
	cout << v.slice<Dynamic,Dynamic>(1, 3) << endl;
	cout << u.slice<Dynamic,Dynamic>(1, 3) << endl << endl;

	cout << v.slice(2, 3) << endl;
	cout << u.slice(2, 3) << endl << endl;
	
	cout << project(v) << endl;
	
}
