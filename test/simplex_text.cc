#include <TooN/optimization/downhill_simplex.h>
using namespace std;
using namespace TooN;

template<class C> double length(const C& v)
{
	return sqrt(v*v);
}

double sq(double x)
{
	return x*x;
}

template<class C> double simplex_size(const C& s)
{
	return abs(length(s.get_simplex()[s.get_best()] - s.get_simplex()[s.get_worst()]) / length( s.get_simplex()[s.get_best()]));
}

double Rosenbrock(const Vector<2>& v)
{
		return sq(1 - v[0]) + 100 * sq(v[1] - sq(v[0]));
}

int main()
{
		Vector<2> starting_point = makeVector( -1, 1);

		DownhillSimplex<2> dh_fixed(Rosenbrock, starting_point, 1);
		while(simplex_size(dh_fixed) > 0.0000001)
		{
				cout << dh_fixed.get_simplex()[dh_fixed.get_best()] << endl;
				cout << dh_fixed.get_values()[dh_fixed.get_best()] << endl;
				dh_fixed.iterate(Rosenbrock);
		}

		cout << dh_fixed.get_simplex()[dh_fixed.get_best()] << endl;
}


