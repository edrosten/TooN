#include <TooN/TooN.h>
using namespace std;
using namespace TooN;

int main()
{
	Vector<6> v = makeVector(1, 2, 3, 4, 5, 6);
	Vector<>  u = v;

	cout << v.slice<0,3>() << endl;
	cout << u.slice<0,3>() << endl << endl;

	cout << v.slice<1,3>(1, 3) << endl;
	cout << u.slice<1,3>(1, 3) << endl;
	cout << v.slice<Dynamic,3>(1, 3) << endl;
	cout << u.slice<Dynamic,3>(1, 3) << endl;
	cout << v.slice<1,Dynamic>(1, 3) << endl;
	cout << u.slice<1,Dynamic>(1, 3) << endl;
	cout << v.slice<Dynamic,Dynamic>(1, 3) << endl;
	cout << u.slice<Dynamic,Dynamic>(1, 3) << endl << endl;

	cout << v.slice(2, 3) << endl;
	cout << u.slice(2, 3) << endl << endl;
	
	cout << project(v) << endl;
	cout << project(u) << endl;
	cout << project(v.slice<1,End<0> >()) << endl;
	cout << project(u.slice<1,End<0> >()) << endl;
	cout << unproject(v) << endl;
	cout << unproject(u) << endl;

	cout << v.slice<1, End<0> >() << endl;	
	cout << u.slice<1, End<0> >() << endl;	
	cout << v.slice<1, End<-1> >() << endl;	
	cout << u.slice<1, End<-1> >() << endl;	
	cout << v.slice(2, End) << endl;	
	cout << u.slice(2, End) << endl;	
	cout << v.slice(2, End(-1)) << endl << endl;	
	cout << u.slice(2, End(-1)) << endl << endl;	

	Vector<200> w = Zeros;

	cout << w.slice<100,End<-99> >() << endl;
}
