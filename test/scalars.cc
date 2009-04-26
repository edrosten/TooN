#include <TooN/TooN.h>
#include <TooN/helpers.h>
using namespace TooN;
using namespace std;
int main()
{
	Vector<5> v = makeVector(1,2,3,4,5);

	cout << v + Scalars(3) << endl;
	cout << v.slice(2,3) + Scalars(3) << endl;

	Matrix<> m = Identity(4);

	cout << m + Scalars(1) << endl;
	cout << m.slice<0,0,2,3>() + Scalars(2) << endl;
}
