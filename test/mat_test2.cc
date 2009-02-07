#include <TooN/TooN.h>
using namespace TooN;
using namespace std;

void statictest()
{
	Matrix<3,4> m;

	for(int i=0; i < 12; i++)
		(&m[0][0])[i] = i;

	cout << m << endl;

	cout << "Accessing rows as vectors:\n";
	for(int i=0; i < 3; i++)
		cout << "...  " << m[i] << "...\n";
	cout << endl;

	cout << "Slice, from [1,1], size [2,3]:\n";
	cout << m.slice<1,1,2,3>() << endl;

	cout << "Accessing rows of a slice as vectors:\n";
	for(int i=0; i < 2; i++)
		cout << "...  " << m.slice<1,1,2,3>()[i] << "...\n";
	cout << endl;



	cout << "Slice, of the above slice from [0,0], size [2,1]:\n";
	cout << m.slice<1,1,2,3>().slice<0,0,2,1>() << endl;

	cout << "Accessing rows of a slice of a slice as vectors:\n";
	for(int i=0; i < 2; i++)
		cout << "...  " << m.slice<1,1,2,3>().slice<0,0,2,1>()[i] << "...\n";
	cout << endl;
}

void dynamictest()
{
	Matrix<> m(3,4);
	for(int i=0; i < 12; i++)
		(&m[0][0])[i] = i;

	cout << m << endl;

	cout << "Accessing rows as vectors:\n";
	for(int i=0; i < 3; i++)
		cout << "...  " << m[i] << "...\n";
	cout << endl;

	cout << "Slice, from [1,1], size [2,3]:\n";
	cout << m.slice<1,1,2,3>() << endl;

	cout << "Accessing rows of a slice as vectors:\n";
	for(int i=0; i < 2; i++)
		cout << "...  " << m.slice<1,1,2,3>()[i] << "...\n";
	cout << endl;



	cout << "Slice, of the above slice from [0,0], size [2,1]:\n";
	cout << m.slice<1,1,2,3>().slice<0,0,2,1>() << endl;

	cout << "Accessing rows of a slice of a slice as vectors:\n";
	for(int i=0; i < 2; i++)
		cout << "...  " << m.slice<1,1,2,3>().slice<0,0,2,1>()[i] << "...\n";
	cout << endl;
}

int main()
{
	cout << "-------------- STATIC MATRIX STATIC SLICE\n";
	statictest();

	cout << "-------------- DYNAMIC MATRIX STATIC SLICE\n";
	dynamictest();
}
