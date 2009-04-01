#include <TooN/TooN.h>
#include <TooN/LU.h>
#include <TooN/helpers.h>
#include <TooN/gaussian_elimination.h>
#include <TooN/gauss_jordan.h>
#include <tr1/random>
#include <sys/time.h>  //gettimeofday
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iomanip>


using namespace TooN;
using namespace std;
using namespace tr1;

double get_time_of_day()
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec+tv.tv_usec * 1e-6;
}

std::tr1::mt19937 eng;
std::tr1::uniform_real<double> rnd;
double global_sum;


struct Do2x2
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		double idet = 1/(a[0][0]*a[1][1] - a[0][1] * a[1][0]);
		double i00 =  a[1][1] * idet;
		double i01 = -a[0][1] * idet;
		double i10 = -a[1][0] * idet;
		double i11 =  a[0][0] * idet;

		for(int i=0; i < x.num_cols(); i++)
		{
			x[0][i] = b[0][i] * i00 + b[1][i] * i01;
			x[1][i] = b[0][i] * i10 + b[1][i] * i11;
		}
	}

	static string name()
	{
		return "2H";
	}
};

struct Do3x3
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		double idet = a[0][0]*(a[2][2]*a[1][1]-a[2][1]*a[1][2])-a[1][0]*(a[2][2]*a[0][1]-a[2][1]*a[0][2])+a[2][0]*(a[1][2]*a[0][1]-a[1][1]*a[0][2]);

		double i00 =   (a[2][2]*a[1][1]-a[2][1]*a[1][2]) * idet;
		double i01 =  -(a[2][2]*a[0][1]-a[2][1]*a[0][2]) * idet;
		double i02 =   (a[1][2]*a[0][1]-a[1][1]*a[0][2]) * idet;
		double i10 =  -(a[2][2]*a[1][0]-a[2][0]*a[1][2]) * idet;
		double i11 =   (a[2][2]*a[0][0]-a[2][0]*a[0][2]) * idet;
		double i12 =  -(a[1][2]*a[0][0]-a[1][0]*a[0][2]) * idet;
		double i20 =   (a[2][1]*a[1][0]-a[2][0]*a[1][1]) * idet;
		double i21 =  -(a[2][1]*a[0][0]-a[2][0]*a[0][1]) * idet;
		double i22 =   (a[1][1]*a[0][0]-a[1][0]*a[0][1]) * idet;

		for(int i=0; i < x.num_cols(); i++)
		{
			x[0][i] = b[0][i] * i00 + b[1][i] * i01 + b[2][i] * i02;
			x[1][i] = b[0][i] * i10 + b[1][i] * i11 + b[2][i] * i12;
			x[2][i] = b[0][i] * i20 + b[1][i] * i21 + b[2][i] * i22;
		}
	}

	static string name()
	{
		return "3H";
	}
};


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
		//x = lu.backsub(b);
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

struct UseGaussianEliminationInverse
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		Matrix<R> i, inv;
		i = Identity;
		inv = gaussian_elimination(a, i);
		x = inv * b;
	}

	static string name()
	{
		return "GI";
	}
};

struct UseGaussJordanInverse
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		Matrix<R, 2*R> m;
		m.template slice<0,0,R,R>() = a;
		m.template slice<0,R,R,R>() = Identity;
		gauss_jordan(m);
		x = m.template slice<0,R,R,R>() * b;
	}

	static string name()
	{
		return "GJ";
	}
};

template<int Size, int Cols, class Solver> void benchmark_ax_eq_b(vector<pair<double, string> >& results)
{
	double time=0, t_tmp, start = get_time_of_day(), t_tmp2;
	double sum=0;
	int n=0;

	while(get_time_of_day() - start < 1)
	{
		Matrix<Size> a;
		for(int r=0; r < Size; r++)
			for(int c=0; c < Size; c++)
				a[r][c] = rnd(eng);
			

		Matrix<Size, Cols> b, x;

		for(int r=0; r < Size; r++)
			for(int c=0; c < Cols; c++)
				b[r][c] = rnd(eng);
		
		a[0][0] += (t_tmp=get_time_of_day()) * 1e-20;
		Solver::template solve<Size, Cols>(a, b, x);
		global_sum += (t_tmp2=get_time_of_day())*x[Size-1][Cols-1];
			
		
		time += t_tmp2 - t_tmp;


		
		for(int r=0; r < Size; r++)
			for(int c=0; c < Cols; c++)
				sum += x[r][c];

		n++;
	}

	results.push_back(make_pair(n/time, Solver::name()));

	global_sum += sum;	
}

template<int Size, int C, bool End=0> struct ColIter
{
	static void iter()
	{
		static const int Lin = Size*2;
		static const int Grow = 2;
		static const int Cols = C + (C<=Lin?0:(C-Lin)*(C-Lin)/Grow);
		vector<pair<double, string> > results;
		cout << Size << "\t" << Cols << "\t";

		benchmark_ax_eq_b<Size, Cols, UseGaussJordanInverse>(results);
		benchmark_ax_eq_b<Size, Cols, UseGaussianElimination>(results);
		benchmark_ax_eq_b<Size, Cols, UseGaussianEliminationInverse>(results);
		benchmark_ax_eq_b<Size, Cols, UseLUInv>(results);
		benchmark_ax_eq_b<Size, Cols, UseLU>(results);
		benchmark_ax_eq_b<Size, Cols, Do2x2>(results);

		sort(results.begin(), results.end());
		for(unsigned int i=0; i < results.size(); i++)
			cout << results[i].second << " " << setprecision(5) << setw(10) << results[i].first << "            ";
		cout << endl;
		ColIter<Size, C+1, (Cols> Size*20)>::iter();
	}
};

template<int Size,int C> struct ColIter<Size, C, 1>
{
	static void iter()
	{
	}
};


template<int Size, bool End=(Size<= 0)> struct SizeIter
{
	static void iter()
	{
		ColIter<Size, 1>::iter();
		SizeIter<Size-16>::iter();
	}
};

template<int S> struct SizeIter<S, 1>
{
	static void iter()
	{
	}
};


int main()
{
	SizeIter<2>::iter();
	
	return global_sum != 123456789.0;
}
