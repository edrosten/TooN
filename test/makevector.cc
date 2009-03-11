#include <TooN/TooN.h>
using namespace TooN;

extern "C"{
double use_make_vector_double(const Vector<4>& v)
{
	return v * makeVector(0,0,2.0,0);
}

double use_make_vector_int(const Vector<4>& v)
{
	return v * makeVector<int>(0,0,2,0);
}

}

