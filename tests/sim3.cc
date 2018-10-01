#include <TooN/sim3.h>
#include <TooN/functions/derivatives.h>
#include <random>

using namespace TooN;
using namespace std;

int main() {
	SIM3<> a(makeVector(1., 0, 0, 0, 0, 0, log(3)));

	cout << a << endl;
	cout << a * a.inverse() << endl;


	a = Identity;
	cout << a << endl;
	
	mt19937 rng;
	normal_distribution<> g(0, 1);
	
	double err = 0;
	
	for(int i=0; i < 10; i++)
	{
		Vector<7> v; 	
		Vector<4> p = Ones;
		for(int i=0; i < 7; i++)
			v[i] = g(rng);

		for(int i=0; i < 3; i++)
			p[i] = g(rng);
	
		SIM3<> a = SIM3<>::exp(v);

		//Compute the derivative using the generator field function
		Matrix<3,7> Jf;
		for(int i=0; i < 7; i++)
			Jf.T()[i] = SIM3<>::generator_field(i, a * p).slice<0,3>();
		
		//Compute the derivative using the generators
		Matrix<3,7> Jg;
		for(int i=0; i < 7; i++)
			Jg.T()[i] = (SIM3<>::generator(i)* a * p).slice<0,3>();
		
		//Compute the derivative numerically

		Matrix<3,7> Jn;
		for(int i=0; i < 3; i++)
			Jn[i] = numerical_gradient([&](const Vector<7>& v){
				return (SIM3<>::exp(v) * a * p)[i];	
			}, Vector<7>(Zeros));
		
		
		err = max(err, max(norm_1(Jg - Jf), norm_1(Jn - Jf)));
	}

	cout << err << endl;


	for(int i=0; i < 10; i++)
	{
		Vector<6> v1;
		Vector<7> v2;

		for(int i=0; i < 6; i++)
			v1[i] = g(rng);
		for(int i=0; i < 7; i++)
			v2[i] = g(rng);

		SE3<> se3{v1};
		SIM3<> sim3{v2};

		
		Matrix<4> r1 = (sim3 * Matrix<4>(Identity)) * (se3*Matrix<4>(Identity));
		Matrix<4> r2 = (sim3 * se3)*Matrix<4>(Identity);

		err = max(err, norm_1(r1 - r2));
	}


}

