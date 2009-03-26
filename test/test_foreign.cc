#include <TooN/TooN.h>
#include <TooN/helpers.h>

using namespace std;
using namespace TooN;

int main()
{
	double data[]={1, 2, 3, 4};

	Matrix<2> m = Identity;

	cout << m << endl;

	cout << Wrap<-1,2>::wrap(data, 1) << endl;

	cout << m + Wrap<-1,2>::wrap(data,2);

};
