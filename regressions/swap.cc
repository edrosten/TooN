#include <TooN/TooN.h>
#include <utility>
#include <cstdlib>
#include "regressions/regression.h"
using namespace TooN;
using namespace std;

template<int Size> void swap_test(int size)
{
	Vector<Size> v1 = Zeros(size);
	Vector<Size> v2 = Zeros(size);

	for(int i=0; i < v1.size(); i++)
	{
		v1[i] = xor128d();
		v2[i] = xor128d();
	}

	Vector<Size> v3 = v1, v4 = v2;

	double* v3d = v3.data(), *v4d = v4.data();
	std::swap(v3, v4);

	if(v3 == v2 && v4 == v1)
	{
		if(Size < 0)
		{
			if(v3.data() == v4d && v4.data() == v3d)
				return;
			else
			{
				cerr << "Inefficient swap on Size=" << Size << " size= " << size << endl;
				exit(1);
			}

		}
	}
	else
	{
		cerr << "Failed on Size=" << Size << " size= " << size << endl;
		exit(1);
	}
}


int main()
{
	
	for(int i=0; i < 10000; i++)
		swap_test<Dynamic>(xor128u() % 100 + 1);

	for(int i=0; i < 10000; i++)
		swap_test<Resizable>(xor128u() % 100 + 1);
	
	for(int i=0; i < 1000; i++)
	{
		swap_test<1>(1);
		swap_test<2>(1);
		swap_test<3>(1);
		swap_test<5>(1);
		swap_test<10>(1);
		swap_test<16>(1);
		swap_test<31>(1);
		swap_test<64>(1);
		swap_test<257>(1);
		swap_test<1023>(1);
	}
}
