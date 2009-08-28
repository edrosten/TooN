#include <iostream>
#include <TooN/TooN.h>
#include <TooN/SymEigen.h>

using namespace std;
using namespace TooN;

int main() {
  // construct M
  double d1[][3] = {{1,2,3},{2,5,6},{3,6,7}};
 Matrix<3> M(3,3);
 M[0]=makeVector(4,0,2);
 M[1]=makeVector(0,5,3);
 M[2]=makeVector(2,3,6);
 
 Vector<3> dg(makeVector(2,4,9));
 // create the eigen decomposition of M
 SymEigen<3> eigM(M);
 cout << "A=" << M << endl;
 cout << "(E,v)=eig(A)" << endl;
 // print the smallest eigenvalue
 cout << "v[0]=" << eigM.get_evalues()[0] << endl;
 // print the associated eigenvector
 cout << "E[0]=" << eigM.get_evectors()[0] << endl;
 // print the square root of the matrix.
 cout << "R=sqrtm(A)=" << eigM.get_sqrtm() << endl;
 // print the square root of the matrix squared.
 cout << "(should equal A), R^T*R="
      << eigM.get_sqrtm().T() * eigM.get_sqrtm() << endl;
 // print the inverse of the matrix.
 cout << "A^-1=" << eigM.get_pinv() << endl;
 // print the inverse square root of the matrix.
 cout << "C=isqrtm(A)=" << eigM.get_isqrtm() << endl;
 // print the inverse square root of the matrix squared.
 cout << "(should equal A^-1), C^T*C="
      << eigM.get_isqrtm().T() * eigM.get_isqrtm() << endl;

  return 0;
}
