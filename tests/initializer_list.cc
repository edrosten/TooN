#include <TooN/TooN.h>
using namespace TooN;
using namespace std;


int main()
{
	Vector<3> v({1, 2, 3});
	cout << v << endl;

	Vector<5> v6 = {{1,2,3,4,5}};
	cout << v6 << endl;
	
	//Not supported
	//Vector<5> v = {1,2,3,4,5};

	Vector<Dynamic> v2 = {1, 2, 3, 4};
	cout << v2 << endl;
	
	Vector<Resizable> v3 = {1, 2, 3, 4, 5};
	cout << v3 << endl;

	Vector<Resizable> v4 = {{1, 2, 3, 4, 5, 6}};
	cout << v4 << endl;

	Vector<Resizable> v5({1,2});
	cout << v5 << endl;
}
