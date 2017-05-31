#include <TooN/TooN.h>
using namespace TooN;
using namespace std;


int main()
{
	Vector<3> v({1, 2, 3});
	cout << v << endl;

	Vector<Dynamic> v2 = {1, 2, 3, 4};
	cout << v2 << endl;
	
	Vector<Resizable> v3 = {1, 2, 3, 4, 5};
	cout << v3 << endl;
}
