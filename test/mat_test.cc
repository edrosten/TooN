#include <TooN/TooN.h>

using namespace TooN;
using namespace std;

int main()
{
	Matrix<3> m;	
	m[0][0] = 1;
	m[0][1] = 2;
	m[0][2] = 3;
	m[1][0] = 4;
	m[1][1] = 5;
	m[1][2] = 6;
	m[2][0] = 7;
	m[2][1] = 8;
	m[2][2] = 9;
	
	cout << m << endl;

	cout << m.slice<0,0,2,3>() << endl;

	cout << m.slice<0,0,2,3>().slice<0,1,2,2>() << endl;
}
