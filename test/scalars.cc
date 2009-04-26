#include <TooN/TooN.h>
#include <TooN/helpers.h>
using namespace TooN;
using namespace std;
int main()
{
	Vector<5> v = makeVector(1,2,3,4,5);

	cout << v + Scalars(3) << endl;
	cout << v.slice(2,3) + Scalars(3) << endl;

	v+=Scalars(1);
	cout << v << endl;
	v.slice(0,3) += Scalars(3);
	cout << v << endl;

	Matrix<> m = Identity(4);
	cout << m + Scalars(1) << endl;
	cout << m.slice<0,0,2,3>() + Scalars(2) << endl;

	m+=Scalars(1);
	cout << m << endl;
	m.slice<0,0,3,2>() += Scalars(2);
	cout << m << endl;
}
