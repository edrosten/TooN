#define TOON_EXPR_TEMPLATES
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
#include <map>
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

#include "solvers.cc"

#define STATIC 1
template<int R, int C, class Precision, class Base> void gauss_jordan_new(Matrix<R, C, Precision, Base>& m)
{
	using std::swap;

	//Loop over columns to reduce.
	for(int col=0; col < m.num_rows(); col++)
	{
		//Reduce the current column to a single element


		//Search down the current column in the lower triangle for the largest
		//absolute element (pivot).  Then swap the pivot row, so that the pivot
		//element is on the diagonal. The benchmarks show that it is actually
		//faster to swap whole rows than it is to access the rows via indirection 
		//and swap the indirection element. This holds for both pointer indirection
		//and using a permutation vector over rows.
		{
		  using std::abs;
			int pivotpos = col;
			double pivotval = abs(m[pivotpos][col]);
			for(int p=col+1; p <m.num_rows(); p++)
			  if(abs(m[p][col]) > pivotval)
				{
					pivotpos = p;
					pivotval = abs(m[pivotpos][col]);
				}
			
			if(col != pivotpos)
				swap(m[col].ref(), m[pivotpos].ref());
		}

		//Reduce the current column in every row to zero, excluding elements on
		//the leading diagonal.
		for(int row = 0; row < m.num_rows(); row++)
		{
			if(row != col)
			{
				double multiple = m[row][col] / m[col][col];
		
				//We could eliminate some of the computations in the augmented
				//matrix, if the augmented half is the identity. In general, it
				//is not. 

				//Subtract the pivot row from all other rows, to make 
				//column col zero.
				m[row][col] = 0;

				int len = m.num_cols()-(col+1);
				m[row].slice(col+1, len) -= m[col].slice(col+1,len)*multiple;

				//for(int c=col+1; c < m.num_cols(); c++)
				//	m[row][c] = m[row][c] - m[col][c] * multiple;
			}
		}
	}
	
	//Final pass to make diagonal elements one. Performing this in a final
	//pass allows us to avoid any significant computations on the left-hand
	//square matrix, since it is diagonal, and ends up as the identity.
	for(int row=0;row < m.num_rows(); row++)
	{
		double mul = 1/m[row][row];

		m[row][row] = 1;

		for(int col=m.num_rows(); col < m.num_cols(); col++)
			m[row][col] *= mul;
	}
}


struct UseCompiledCramer
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{	
		#ifdef STATIC
			solve_direct(a, b, x);
		#else

		#endif
	}

	static string name()
	{
		return "CC";
	}
};

struct UseLU
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		LU<STATIC?R:Dynamic> lu(a);

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
		LU<STATIC?R:Dynamic> lu(a);

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
		#ifdef STATIC
			x = gaussian_elimination(a, b);
		#else
			x = gaussian_elimination(a.slice(0,0,a.num_rows(), a.num_cols()), b.slice(0, b.size()));
		#endif
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
		Matrix<STATIC?R:Dynamic> i=Identity(a.num_rows()), inv(a.num_rows(), a.num_rows());
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
		Matrix<STATIC?R:Dynamic, STATIC?2*R:Dynamic> m(a.num_rows(), a.num_rows()*2);
		m.template slice<0,0,STATIC?R:Dynamic,STATIC?R:Dynamic>(0,0,a.num_rows(), a.num_rows()) = a;
		m.template slice<0,STATIC?R:Dynamic,STATIC?R:Dynamic,STATIC?R:Dynamic>(0, a.num_rows(), a.num_rows(), a.num_rows()) = Identity;
		gauss_jordan(m);
		x = m.template slice<0,R,R,R>() * b;
	}

	static string name()
	{
		return "GJ";
	}
};


struct UseGaussJordanNewInverse
{
	template<int R, int C> static void solve(const Matrix<R, R>& a, const Matrix<R, C>& b, Matrix<R, C>& x)
	{
		Matrix<STATIC?R:Dynamic, STATIC?2*R:Dynamic> m(a.num_rows(), a.num_rows()*2);
		m.template slice<0,0,STATIC?R:Dynamic,STATIC?R:Dynamic>(0,0,a.num_rows(), a.num_rows()) = a;
		m.template slice<0,STATIC?R:Dynamic,STATIC?R:Dynamic,STATIC?R:Dynamic>(0, a.num_rows(), a.num_rows(), a.num_rows()) = Identity;
		gauss_jordan_new(m);
		x = m.template slice<0,R,R,R>() * b;
	}

	static string name()
	{
		return "GN";
	}
};


template<int Size, int Cols, class Solver> void benchmark_ax_eq_b(map<string, vector<double> >& results)
{
	double time=0, t_tmp, start = get_time_of_day(), t_tmp2;
	double sum=0;
	int n=0;

	while(get_time_of_day() - start < .1)
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
		n++;
	}

	results[Solver::name()].push_back(n/time);

	global_sum += sum;	
}




template<int Size, int Cols, typename Solver, bool Use> struct Optional
{
	static void solve(map<string, vector<double> >& r)
	{
		benchmark_ax_eq_b<Size, Cols, Solver>(r);
	}
};


template<int Size, int Cols, typename Solver > struct Optional<Size, Cols, Solver, 0>
{
	static void solve(map<string, vector<double> >&)
	{
	}
};

template<int Size, int C=1, bool End=0> struct ColIter
{
	static void iter()
	{
		static const int Lin = Size*2;
		static const int Grow = 1;
		static const int Cols = C + (C<=Lin?0:(C-Lin)*(C-Lin)*(C-Lin)/Grow);
		map<string, vector<double> > results;
		cout << Size << "\t" << Cols << "\t";
		
		//Run each menchmark 10 times and select the median result
		for(int i=0; i < 10; i++)
		{
			benchmark_ax_eq_b<Size, Cols, UseGaussJordanInverse>(results);
			benchmark_ax_eq_b<Size, Cols, UseGaussJordanNewInverse>(results);
			//benchmark_ax_eq_b<Size, Cols, UseGaussianElimination>(results);
			//benchmark_ax_eq_b<Size, Cols, UseGaussianEliminationInverse>(results);
			//benchmark_ax_eq_b<Size, Cols, UseLUInv>(results);
			//benchmark_ax_eq_b<Size, Cols, UseLU>(results);
			//Optional<Size, Cols, UseCompiledCramer, (Size<=highest_solver)>::solve(results);
		}
		
		vector<pair<double, string> > res;
		for(map<string, vector<double> >::iterator i=results.begin(); i != results.end(); i++)
		{
			sort(i->second.begin(), i->second.end());
			res.push_back(make_pair(i->second[i->second.size()/2], i->first));
		}



		sort(res.begin(), res.end());
		for(unsigned int i=0; i < res.size(); i++)
			cout << res[i].second << " " << setprecision(5) << setw(10) << res[i].first << "            ";
		cout << endl;
		ColIter<Size, C+1, (Cols> Size*50)>::iter();
	}
};

template<int Size, int C> struct ColIter<Size, C, 1> 
{

	static void iter()
	{
	}
};

#ifndef SIZE
	#define SIZE 5
#endif

int main()
{
	ColIter<SIZE>::iter();
	
	return global_sum != 123456789.0;
}
