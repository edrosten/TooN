#include <TooN/TooN.h>
using namespace TooN;
using namespace std;

int main()
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
