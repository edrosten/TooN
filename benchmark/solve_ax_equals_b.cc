#include <TooN/TooN.h>
#include <TooN/LU.h>
#include <TooN/gaussian_elimination.h>
#include <cvd/timer.h>
#include <tr1/random>

using namespace TooN;
using namespace std;
using namespace tr1;
using namespace CVD;

std::tr1::mt19937 eng;
std::tr1::uniform_real<double> rnd;
double global_sum;

struct UseLU
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		LU<R> lu(a);

		x = lu.backsub(b);
	}

	static string name()
	{
		return "LU";
	}
};

struct UseLUInv
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		LU<R> lu(a);

		x = lu.get_inverse() * b;
	}

	static string name()
	{
		return "LI";
	}
};

struct UseGaussianElimination
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		x = gaussian_elimination(a, b);
	}

	static string name()
	{
		return "GE";
	}
};

template<int Size, int Cols, class Solver> void benchmark_ax_eq_b()
{
	cvd_timer t;
	double time=0;
	double sum=0;
	int n=0;

	while(time < 1)
	{
		Matrix<Size> a;
		for(int r=0; r < Size; r++)
			for(int c=0; c < Size; c++)
				a[r][c] = rnd(eng);
			

		Matrix<Size, Cols> b, x;

		for(int r=0; r < Size; r++)
			for(int c=0; c < Cols; c++)
				b[r][c] = rnd(eng);
		
		t.reset();
		Solver::template solve<Size, Cols>(a, b, x);
		time += t.get_time();

		
		for(int r=0; r < Size; r++)
			for(int c=0; c < Cols; c++)
				sum += x[r][c];

		n++;
	}

	cout << Solver::name() << " " << Size << " " << Cols << " " << n / time << endl;

	global_sum += sum;	
}

template<class A, class B> struct TypeList
{
	typedef A a;
	typedef B b;
};

struct Null{};

template<int R, int C, class L> struct benchmark_iter
{
	static void iter()
	{
		benchmark_ax_eq_b<R, C, typename L::a>();
		benchmark_iter<R, C, typename L::b>::iter();
	}
};


template<int R, int C> struct benchmark_iter<R, C, Null>
{
	static void iter()
	{
	}
};




template<int Size, int Cols, class Test> struct ColIter
{
	static void iter()
	{
		benchmark_iter<Size, Cols, Test>::iter();
		ColIter<Size, Cols-1, Test>::iter();
	}
};

template<int Size,class Test> struct ColIter<Size, 1, Test>
{
	static void iter()
	{
	}
};


template<int Size, class Test> struct SizeIter
{
	static void iter()
	{
		ColIter<Size, Size*8+1, Test>::iter();
		SizeIter<Size-1, Test>::iter();
	}
};

template<class Test> struct SizeIter<0, Test>
{
	static void iter()
	{
	}
};


int main()
{
	SizeIter<5, TypeList<UseGaussianElimination, TypeList<UseLU, Null> > >::iter();
	
	return global_sum != 123456789.0;
}
