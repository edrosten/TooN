#include <cstdlib>
#include <iostream>
#include <TooN/TooN.h>
#include <TooN/LU.h>

using namespace TooN;
using namespace std;


int main()
{
	Matrix<-1,5,float> m(5,5);

	for(int i=0; i< m.num_rows(); i++)
		for(int j=0; j< m.num_rows(); j++)
			m[i][j] = drand48();

	
	LU<5,float> mlu(m);

	cout << mlu.get_inverse() << endl;
	Matrix<5,5,float> inv = mlu.get_inverse();
	Matrix<5,5,float> a = m*inv;

	for(int i=0; i< m.num_rows(); i++)
		for(int j=0; j< m.num_rows(); j++)
		{
			if(round(a[i][j]) < 1e-10)
				a[i][j] = 0;
		}

	cout << a;



}
