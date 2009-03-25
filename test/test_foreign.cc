#include <TooN/TooN.h>
#include <TooN/helpers.h>

using namespace std;
using namespace TooN;

int main()
{
	double data[]={1, 2, 3, 4};

	Matrix<2> m = Identity;

	cout << m << endl;

	Matrix<2, 2, double, RowMajorContigRef> n(data);

	cout << n << endl;

	cout << m + n;

};
