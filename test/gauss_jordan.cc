#include <TooN/gauss_jordan.h>
#include <TooN/helpers.h>
#include <tr1/random>
#include <fstream>

using namespace TooN;
using namespace std;
using namespace tr1;

int main()
{
	unsigned int s;
	ifstream("/dev/urandom").read((char*)&s, 4);
	
	std::tr1::mt19937 eng;
	std::tr1::uniform_real<double> rnd;
	eng.seed(s);
	Matrix<5,10> m, reduced;

	for(int i=0; i< m.num_rows(); i++)
		for(int j=0; j< m.num_rows(); j++)
			m[i][j] = rnd(eng);

	m.slice<0, 5, 5, 5>()  = Identity;

	cout << m << endl;


	reduced = m;
	gauss_jordan(reduced);

	cout << reduced.slice<0,5,5,5>() * m.slice<0,0,5,5>() << endl;

}
