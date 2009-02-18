#include <TooN/TooN.h>
#include <TooN/helpers.h>

using namespace std;
using namespace TooN;

template<class C> void type(const C&)
{
	cout << __PRETTY_FUNCTION__ << endl;
}

int main()
{
	Matrix<3> m1;
	Matrix<3> m2;

	Zero(m1);
	Zero(m2);
	
	m1.slice<0,0,2,2>()+=3;
	m2.slice<1,1,2,2>()+=2;

	cout << m1 << endl;
	cout << m2 << endl;
	
	m1+=m2;
	cout << m1 << endl;
}

