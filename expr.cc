#include <TooN/TooN.h>

using namespace TooN;
using namespace std;
extern "C"{
void foo(const Vector<>& v1, const Vector<>& v2, Vector<>& v3, int i)
{
	v3 =-(v1+2.4*v2).slice(0,i)/2.8;
}
}

/*
int main()
{
	Vector<> v1 = makeVector(1,2,3);
	Vector<> v2 = makeVector(6,3,1);

	Vector<> v3 = -3*(v1 + 2*v2).slice(0,2) / 2;


	typedef decltype(-3*(v1 + 2*v2).slice(0,2) / 2) A;

	cout << (v1+v2)[0] << endl;
	cout << v3 << endl;
	
	
	cout <<  (-3*(v1 + 2*v2).slice(0,2) / 2).cost() << endl;

//	cout << IsExpressionTemplate<decltype(v3)>::is << endl;
	//cout << IsExpressionTemplate<decltype(v1+v2)>::is << endl;
}*/
