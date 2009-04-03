#include <TooN/TooN.h>
#include <TooN/sl.h>

#include <iostream>
#include <iomanip>

using namespace TooN;
using namespace std;

int main(int , char ** ){
	SL<3> h(makeVector(1,0,-1,0,0,0,0,0));
	cout << h << endl;
	cout << h.inverse() << endl;
	cout << SL<3>::exp(makeVector(-1,0,1,0,0,0,0,0)) << endl;
	cout << h * h.inverse() << endl;
	
	for(int i = 0; i < SL<3>::dim; ++i)
		cout << "generator " << i << "\n" << SL<3>::generator(i) << endl;

	for(int i = 0; i < SL<5>::dim; ++i)
		cout << "generator " << i << "\n" << SL<5>::generator(i) << endl;
	
	cout << SL<2>::exp(makeVector(1,2,3)) << endl;
	
	return 0;
}
