#include <TooN/TooN.h>
#include <TooN/helpers.h>

using namespace std;
using namespace TooN;

int main()
{
	double data[]={1, 2, 3, 4, 5, 6};


	cout << Matrix<2,3,double, Reference::RowMajor> (data) << endl;
	cout << Matrix<2,-1,double, Reference::RowMajor> (data,2,3) << endl;
	cout << Matrix<-1,3,double, Reference::RowMajor> (data,2,3) << endl;
	cout << Matrix<-1,-1,double, Reference::RowMajor> (data,2,3) << endl;

	cout << Wrap<2,3>::wrap(data) << endl;
	cout << Wrap<2,Dynamic>::wrap(data,3) << endl;
	cout << Wrap<Dynamic,3>::wrap(data,2) << endl;
	cout << Wrap<Dynamic, Dynamic>::wrap(data,2,3) << endl;


	cout << Matrix<2,3,double, Reference::ColMajor> (data) << endl;
	cout << Matrix<2,-1,double, Reference::ColMajor> (data,2,3) << endl;
	cout << Matrix<-1,3,double, Reference::ColMajor> (data,2,3) << endl;
	cout << Matrix<-1,-1,double, Reference::ColMajor> (data,2,3) << endl;

	cout << Wrap<2,3,double,Reference::ColMajor>::wrap(data) << endl;
	cout << Wrap<2,Dynamic,double,Reference::ColMajor>::wrap(data,3) << endl;
	cout << Wrap<Dynamic,3,double,Reference::ColMajor>::wrap(data,2) << endl;
	cout << Wrap<Dynamic, Dynamic,double,Reference::ColMajor>::wrap(data,2,3) << endl;

};
