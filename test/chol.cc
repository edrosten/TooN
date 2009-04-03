#include <TooN/Cholesky.h>

#include <iostream>

using namespace std;

int main(int, char ** ){
    
    TooN::Matrix<3> t;
    t[0] = TooN::makeVector(1,0.5, 0.5);
    t[1] = TooN::makeVector(0.5, 2, 0.7);
    t[2] = TooN::makeVector(0.5,0.7, 3);
    
    TooN::Cholesky<3> chol(t);
    cout << chol.get_determinant() << endl;
    
    cout << t << "\n" <<  chol.get_inverse() << "\n" << t * chol.get_inverse() << endl;
    cout << chol.get_L() << endl;

    TooN::Matrix<> t2(3,3);
    t2[0] = TooN::makeVector(1,0.5, 0.5);
    t2[1] = TooN::makeVector(0.5, 2, 0.7);
    t2[2] = TooN::makeVector(0.5,0.7, 3);
    
    TooN::Cholesky<-1,float> chol2(t2);
    
    cout << t2 << "\n" <<  chol2.get_inverse() << "\n" << t2 * chol2.get_inverse() << endl;
    cout << chol2.get_L() << endl;
    
    TooN::Vector<3> bla = TooN::makeVector(1,2,3);
    
    cout << chol.backsub(bla) << endl;
    cout << chol.get_inverse() * bla << endl;
    cout << chol.mahalanobis(bla) - bla * chol.backsub(bla) << endl;
    
    return 0;
}
