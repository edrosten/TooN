#include <TooN/SVD.h>
#include <TooN/helpers.h>
using namespace TooN;
using namespace std;

int main()
{
	Matrix<4, 4> m = Zero;
	m[0] = makeVector(1, 2, 3, 4);
	m[1] = makeVector(1, 1, 1, 1);

	SVD<4, 4> svdm(m);

	cout << svdm.get_VT().num_rows() << endl;
	cout << svdm.get_VT().num_cols() << endl;

	cout << m[0] * svdm.get_VT()[2] << endl;
	cout << m[0] * svdm.get_VT()[3] << endl;
	cout << m[1] * svdm.get_VT()[2] << endl;
	cout << m[1] * svdm.get_VT()[3] << endl;

}
