#include "regressions/regression.h"

int main()
{
	cout << setprecision(10) 
	     << LU<2>(Matrix<2>(Data(1, 2, 3, 4))).get_inverse();
}
