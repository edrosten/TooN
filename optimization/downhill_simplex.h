#ifndef TOON_DOWNHILL_SIMPLEX_H
#define TOON_DOWNHILL_SIMPLEX_H
#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <algorithm>
#include <vector>
#include <math.h>

namespace TooN
{

template<int D> struct DSBase
{
	typedef Vector<D> Vec;
	typedef Vector<D+1> Values;
	typedef std::vector<Vector<D> > Simplex;

	static const int Dim = D;


	DSBase(int) { };

	void resize_simplex(Simplex&s) {
		s.resize(Dim+1);
	}
	void resize_values(Values&) {}
	void resize_vector(Vec&) {}

};

template<> struct DSBase<-1>
{
	typedef Vector<> Vec;
	typedef Vector<> Values;
	typedef std::vector<Vector<> > Simplex;
	int Dim;

	DSBase(int d)
	{
		Dim = d;
	};

	void resize_simplex(Simplex& s)
	{
		s.resize(Dim+1, Vector<>(Dim));
	}

	void resize_values(Values& v)
	{
		v.resize(Dim+1);
	}

	void resize_vector(Vec& v)
	{
		v.resize(Dim);
	}
};

/** This is an implementation of the Downhill Simplex (Nelder & Mead, 1965)
    algorithm. This particular instance will minimize a given function.
	
	The function maintains \f$N+1\f$ points for a $N$ dimensional function, \f$f\f$
	
	At each iteration, the following algorithm is performed:
	- Find the worst (largest) point, \f$x_w\f$.
	- Find the centroid of the remaining points, \f$x_0\f$.
	- Let \f$v = x_0 - x_w\f$
	- Compute a reflected point, \f$ x_r = x_0 + \alpha v\f$
	- If \f$f(x_r)\f$ is better than the best point
	  - Expand the simplex by extending the reflection to \f$x_e = x_0 + \rho \alpha v \f$
	  - Replace \f$x_w\f$ with the best point of \f$x_e\f$,  and \f$x_r\f$.
	- Else, if  \f$f(x_r)\f$ is between the best and second-worst point
	  - Replace \f$x_w\f$ with \f$x_r\f$.
	- Else, if  \f$f(x_r)\f$ is better than \f$x_w\f$
	  - Contract the simplex by computing \f$x_c = x_0 + \gamma v\f$
	  - If \f$f(x_c) < f(x_r)\f$
	    - Replace \f$x_w\f$ with \f$x_c\f$.
	- If \f$x_w\f$ has not been replaced, then shrink the simplex by a factor of \f$\sigma\f$ around the best point.

	This implementation uses:
	- \f$\alpha = 1\f$
	- \f$\rho = 2\f$
	- \f$\gamma = 1/2\f$
	- \f$\sigma = 1/2\f$
	
	Example usage:
	@code
	#include <utility>
	#include <TooN/optimization/downhill_simplex.h>
	using namespace std;
	using namespace TooN;

	template<class C> double length(const C& v)
	{
        return sqrt(v*v);
	}

	template<class C> double simplex_size(const C& s)
	{
		return abs(length(s.get_simplex()[s.get_best()] - s.get_simplex()[s.get_worst()]) / length( s.get_simplex()[s.get_best()]));
	}

	double Rosenbrock(Vector<2>& v)
	{
			return sq(1 - v[0]) + 100 * sq(v[1] - sq(v[0]));
	}

	int main()
	{
			Vector<2> starting_point = (make_Vector, -1, 1);

			DownhillSimplex<2> dh_fixed(RosenbrockF, starting_point, 1);
			while(simplex_size(dh_fixed) > 0.0000001)
					dh_fixed.iterate(RosenbrockF);

			cout << dh_fixed.get_simplex()[dh_fixed.get_best()] << endl;
	}
	@endcode


    @ingroup gOptimize
	@param   N The dimension of the function to optimize. As usual, the default value of <i>N</i> (-1) indicates
	         that the class is sized at run-time.


**/
template<int N=-1> class DownhillSimplex: public DSBase<N>
{
	typedef typename DSBase<N>::Vec Vec;
	typedef typename DSBase<N>::Values Values;
	typedef typename DSBase<N>::Simplex Simplex;

	using DSBase<N>::Dim;

	public:
		/// Initialize the DownhillSimplex class. The simplex is automatically
		/// generated. One point is at <i>c</i>, the remaining points are made by moving
		/// <i>c</i> by <i>spread</i> along each axis aligned unit vector.
		///
		///@param func       Functor to minimize.
		///@param c          Origin of the initial simplex. The dimension of this vector
		///                  is used to determine the dimension of the run-time sized version.
		///@param spread     Size of the initial simplex.
		template<class Function> DownhillSimplex(const Function& func, const Vec& c, double spread=1)
		:DSBase<N>(c.size())
		{
			resize_simplex(simplex);
			resize_values(values);

			for(int i=0; i < Dim+1; i++)
				simplex[i] = c;

			for(int i=0; i < Dim; i++)
				simplex[i][i] += spread;

			alpha = 1.0;
			rho = 2.0;
			gamma = 0.5;
			sigma = 0.5;

			for(int i=0; i < Dim+1; i++)
				values[i] = func(simplex[i]);
		}

		///Return the simplex
		const Simplex& get_simplex() const
		{
			return simplex;
		}
		
		///Return the score at the vertices
		const Values& get_values() const
		{
			return values;
		}
		
		///Get the index of the best vertex
		int get_best() const 
		{
			return min_element(values.begin(), values.end()) - values.begin();
		}
		
		///Get the index of the worst vertex
		int get_worst() const 
		{
			return max_element(values.begin(), values.end()) - values.begin();
		}

		///Perform one iteration of the downhill Simplex algorithm
		///@param func Functor to minimize
		template<class Function> void iterate(const Function& func)
		{
			//Find various things:
			// - The worst point
			// - The second worst point
			// - The best point
			// - The centroid of all the points but the worst
			int worst = get_worst();
			double second_worst_val=-HUGE_VAL, bestval = HUGE_VAL, worst_val = values[worst];
			int best=0;
			Vec x0;
			resize_vector(x0);
			Zero(x0);


			for(int i=0; i < Dim+1; i++)
			{
				if(values[i] < bestval)
				{
					bestval = values[i];
					best = i;
				}

				if(i != worst)
				{
					if(values[i] > second_worst_val)
						second_worst_val = values[i];

					//Compute the centroid of the non-worst points;
					x0 += simplex[i];
				}
			}
			x0 *= 1.0 / Dim;


			//Reflect the worst point about the centroid.
			Vec xr = (1 + alpha) * x0 - alpha * simplex[worst];
			double fr = func(xr);

			if(fr < bestval)
			{
				//If the new point is better than the smallest, then try expanding the simplex.
				Vec xe = rho * xr + (1-rho) * x0;
				double fe = func(xe);

				//Keep whichever is best
				if(fe < fr)
				{
					simplex[worst] = xe;
					values[worst] = fe;
				}
				else
				{
					simplex[worst] = xr;
					values[worst] = fr;
				}

				return;
			}

			//Otherwise, if the new point lies between the other points
			//then keep it and move on to the next iteration.
			if(fr < second_worst_val)
			{
				simplex[worst] = xr;
				values[worst] = fr;
				return;
			}


			//Otherwise, if the new point is a bit better than the worst point,
			//(ie, it's got just a little bit better) then contract the simplex
			//a bit.
			if(fr < worst_val)
			{
				Vec xc = (1 + gamma) * x0 - gamma * simplex[worst];
				double fc = func(xc);

				//If this helped, use it
				if(fc <= fr)
				{
					simplex[worst] = xc;
					values[worst] = fc;
					return;
				}
			}
			
			//Otherwise, fr is worse than the worst point, or the fc was worse
			//than fr. So shrink the whole simplex around the best point.
			for(int i=0; i < Dim+1; i++)
				if(i != best)
				{
					simplex[i] = simplex[best] + sigma * (simplex[i] - simplex[best]);
					values[i] = func(simplex[i]);
				}
		}





	private:
		float alpha, rho, gamma, sigma;

		//Each row is a simplex vertex
		Simplex simplex;

		//Function values for each vertex
		Values  values;


};
}
#endif
