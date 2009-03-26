#include <iostream>

using namespace std;

#include <TooN/so2.h>
#include <TooN/se2.h>

void test_so2(){
    TooN::SO2<> r1;
    cout << "default constructor\n";
    cout << r1 << endl;
    cout << "default constructor <int>\n";
    TooN::SO2<int> r2;
    cout << r2 << endl;
    TooN::SO2<> r(0.1);
    cout << "constructor with 0.1\n";
    cout << r << endl;
    cout << "generator\n";
    cout << r.generator() << endl;
    cout << "ln()\n";
    cout << r.ln() << endl;
    cout << "inverse\n";
    cout << r.inverse() << endl;
    cout << "times\n";
    cout << r * r.inverse() << endl;
    cout << (r *= r.inverse()) << endl;
    r = TooN::SO2<>::exp(0.1);

    TooN::Vector<2> t = TooN::makeVector(0,1);
    cout << "right and left multiply with vector " << t << "\n";
    cout << r * t << endl;
    cout << t * r << endl;
    TooN::Matrix<2> l = TooN::Identity;
    cout << "right and left multiply with matrix\n" << l << "\n";
    cout << r * l << endl;
    cout << l * r << endl;
    TooN::Matrix<2,3> l2(TooN::Zero);
    l2[0] = TooN::makeVector(0,1,2);
    cout << "right with rectangular matrix\n";
    cout << r * l2 << endl;
    
    TooN::Matrix<2> m;
    m[0] = TooN::makeVector(0.5, 1);
    m[1] = TooN::makeVector(1,1);
    cout << "set from matrix (uses coerce) " << m << "\n";
    r = m;
    cout << r << endl;
}

void test_se2(){
    TooN::SE2<> r1;
    cout << "default constructor\n";
    cout << r1 << endl;
    cout << "default constructor <int>\n";
    TooN::SE2<int> r2;
    cout << r2 << endl;
    
    TooN::SE2<> r3(TooN::makeVector(1,1,1));
    cout << "from vector 1 1 1\n";
    cout << r3 << endl;
    cout << r3.ln() << endl;
    
    cout << "generators 0,1,2\n";
    cout << TooN::SE2<>::generator(0) ;
    cout << TooN::SE2<>::generator(1) ;
    cout << TooN::SE2<>::generator(2) << endl;
}

int main(int, char* *){
 
    test_so2();
    test_se2();
 
    return 0;
}
