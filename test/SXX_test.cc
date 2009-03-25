#include <iostream>

using namespace std;

#include <TooN/so2.h>

int main(int, char *){
    TooN::SO2<> r(M_PI_2);
    cout << r << endl;
    cout << r.generator() << endl;

    TooN::Vector<2, float> t = TooN::makeVector(0,1);
    cout << r * t << endl;
    cout << t * r << endl;
    TooN::Matrix<2> l;
    cout << r * l << endl;
    cout << l * r << endl;
    return 0;
}
