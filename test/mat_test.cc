#include <TooN/TooN.h>
#include <iterator>
#include <iostream>

using namespace TooN;
using namespace std;

template<int R, int C, class P, class L> void foo(Matrix<R,C,P,L>& m)
{
	cout << "In foo:\n";
	m[0][0] = 1;
	m[0][1] = 2;
	m[0][2] = 3;
	m[0][3] = 4;
	m[1][0] = 5;
	m[1][1] = 6;
	m[1][2] = 7;
	m[1][3] = 8;
	m[2][0] = 9;
	m[2][1] = 10;
	m[2][2] = 11;
	m[2][3] = 12;
	
	cout << "Slices:\n";
	cout << m << endl;
	cout << m.template slice<0,0,2,3>() << endl;
	cout << m.template slice<0,0,2,3>().template slice<0,1,2,2>() << endl;

	cout << "Transposes:\n";
	cout << m.T() << endl;
	cout << m.T().template slice<0,0,2,3>() << endl;
	cout << m.T().template slice<0,0,2,3>().template slice<0,1,2,2>() << endl;

	cout << "Subscripts:\n";
	cout << m[0] << ", " << m[1] << ", " << m[2] << endl;

	cout << m.template slice<1,1,2,2>()[0] << ", " << m.template slice<1,1,2,2>()[1] << endl;

	cout << "Layout:\n";

	copy(m.data(), m.data()+m.num_rows()*m.num_cols(), ostream_iterator<double>(cout, " "));
	cout << endl;
}


int main()
{
	Matrix<3,4> m;	
	foo(m);
	Matrix<3,4,double,ColMajor> n;	
	foo(n);
}
