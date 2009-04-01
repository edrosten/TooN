#include <iostream>

using namespace std;

#include <TooN/so2.h>
#include <TooN/se2.h>
#include <TooN/so3.h>

void test_so2(){
    cout << "---------------SO2 Tests---------------\n";
    
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
    cout << l2.T() * r << endl;

    TooN::Matrix<2> m;
    m[0] = TooN::makeVector(0.5, 1);
    m[1] = TooN::makeVector(1,1);
    cout << "set from matrix (uses coerce) " << m << "\n";
    r = m;
    cout << r << endl;
}

void test_se2(){
    cout << "---------------SE2 Tests---------------\n";

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

    TooN::Vector<2> t1 = TooN::makeVector(0,1);
    TooN::Vector<> t2(3); t2 = TooN::makeVector(1,0,1);
    cout << "se2 * vector\n";
    cout << r3 * t1 << endl;
    cout << r3 * t2 << endl;
    cout << "vector * se3\n";    
    // cout << t1 * r3 << endl; // this is not well defined, should the output be a 3 vector ?
    cout << t2 * r3 << endl;    

#if 0    
    TooN::Matrix<3> m1;
    TooN::Matrix<> m2(3,3);
    TooN::Matrix<3,10> m3;
    cout << "se2 * matrix\n";
    cout << r3 * m1 << endl;
    cout << r3 * m2 << endl;
    cout << r3 * m3 << endl;
    cout << "matrix * se2\n";
    cout << m1 * r3 << endl;
    cout << m2 * r3 << endl;
    cout << m3.T() * r3 << endl;
#endif
}

void test_so3(){
    cout << "---------------SO3 Tests---------------\n";
    
    TooN::SO3<> r1;
    cout << "default constructor\n";
    cout << r1 << endl;
    cout << "default constructor <int>\n";
    TooN::SO3<int> r2;
    cout << r2 << endl;
    TooN::SO3<> r(TooN::makeVector(0.1, 0.1, 0.1));
    cout << "constructor with 0.1\n";
    cout << r << endl;
    cout << "generator 0,1,2\n";
    cout << TooN::SO3<>::generator(0) ;
    cout << TooN::SO3<>::generator(1) ;
    cout << TooN::SO3<>::generator(2) << endl;
    cout << "ln()\n";
    cout << r.ln() << endl;
    cout << "inverse\n";
    cout << r.inverse() << endl;
    cout << "times\n";
    cout << r * r.inverse() << endl;
    cout << (r *= r.inverse()) << endl;
    r = TooN::SO3<>::exp(TooN::makeVector(0.1, 0.1, 0.1));

    TooN::Vector<3> t = TooN::makeVector(0,1,2);
    cout << "right and left multiply with vector " << t << "\n";
    cout << r * t << endl;
    cout << t * r << endl;
    TooN::Matrix<3> l = TooN::Identity;
    cout << "right and left multiply with matrix\n" << l << "\n";
    cout << r * l << endl;
    cout << l * r << endl;
    TooN::Matrix<3,6> l2(TooN::Zero);
    l2[0] = TooN::makeVector(0,1,2,3,4,5);
    cout << "right with rectangular matrix\n";
    cout << r * l2 << endl;
    cout << l2.T() * r << endl;

    TooN::Matrix<3> m;
    m[0] = TooN::makeVector(0.5, 1,2);
    m[1] = TooN::makeVector(1,1,0);
    m[2] = TooN::makeVector(0,-1,0);
    cout << "set from matrix (uses coerce)\n" << m << "\n";
    r = m;
    cout << r << endl;
}

int main(int, char* *){
 
    test_so2();
    test_so3();
 
    return 0;
}
