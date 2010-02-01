#include <complex>
#include "regressions/regression.h"
#include <TooN/internal/planar_complex.hh>

int main()
{
	complex<double> i(0,1);

	Vector<3, complex<double> > v1 = makeVector<complex<double> >(1.+i, 1.+2.*i, 3);
	Vector<3, complex<double> > v2 = makeVector<complex<double> >(1.-i, 1.-2.*i, 3);

	cout << v2 * v1 << endl;
	cout << v2.as_diagonal() * v1 << endl;


	double re2[] = {1, 1, 3};	
	double im2[] = {-1, -2, 0};	

	Vector<3, complex<double>, ReferencePlanarComplex> v2ish(make_pair(re2, im2));
	cout << v1 * v2ish << endl;

	double real[] = {1,2,3,4};
	double imag[] = {5,6,7,8};


	Vector<4, complex<double>, ReferencePlanarComplex> vec(make_pair(real, imag));
	
	cout << vec << endl;
	cout << vec.slice<1,3>() << endl;
	cout << vec.slice(2,2) << endl;

	real[3] = 28;
	imag[3] = 10;
		
	cout << vec << endl;

}
