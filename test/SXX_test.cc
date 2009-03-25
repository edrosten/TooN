#include <iostream>

using namespace std;

#include <TooN/so2.h>

int main(int, char* *){
    TooN::SO2<> r1;
    cout << r1 << endl;
    TooN::SO2<int> r2;
    cout << r2 << endl;
    TooN::SO2<> r(0.1);
    cout << r << endl;
    cout << r.generator() << endl;
    cout << r.ln() << endl;

    TooN::Vector<2, float> t = TooN::makeVector(0,1);
    cout << r * t << endl;
    cout << t * r << endl;
    TooN::Matrix<2> l;
    cout << r * l << endl;
    cout << l * r << endl;
    
    TooN::Matrix<2> m;
    m[0] = TooN::makeVector(0.5, 1);
    m[1] = TooN::makeVector(1,1);
    r = m;
    cout << r << endl;
 
    TooN::Matrix<> test(TooN::Identity);
 
 return 0;
}
