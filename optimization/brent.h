#ifndef TOON_BRENT_H
#define TOON_BRENT_H
#include <TooN/TooN.h>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <iomanip>


#include <pstreams/pstream.h>
#include <sstream>

class Plotter
{
	std::vector<std::string> plots, pdata;
	redi::opstream plot;

	public:
		Plotter()
		:plot("gnuplot")
		{}


		std::ostream& s()
		{
			return plot;
		}

		Plotter& newline(const std::string& ps)
		{
			plots.push_back(ps);
			pdata.push_back("");
			return *this;
		}
		
		template<class C> Plotter& addpt(const C& pt)
		{
			std::ostringstream o;
			o << pt << std::endl;
			pdata.back() += o.str();
			return *this;
		}
		
		template<class C> Plotter& addpt(C p1, C p2)
		{
			std::ostringstream o;
			o << p1 << " " << p2 << std::endl;
			pdata.back() += o.str();
			return *this;
		}

		template<class C> Plotter& addpt(const std::vector<C>& pt)
		{
			std::ostringstream o;
			for(unsigned int i=0; i < pt.size(); i++)
				o << pt[i] << std::endl;
			pdata.back() += o.str();
			return *this;
		}
	
		Plotter& skip()
		{
			pdata.back() += "\n";
			return *this;
		}

		void draw()
		{
			for(unsigned int i=0; i < plots.size(); i++)
			{
				if(i == 0)
					plot << "plot \"-\"";
				else
					plot << ", \"-\"";

				if(plots[i] != "")
					plot << " with " << plots[i];
			}

			plot << std::endl;
			for(unsigned int i=0; i < plots.size(); i++)
				plot << pdata[i] << "e\n";
			plot << std::flush;
			
			plots.clear();
			pdata.clear();
		}
};

namespace TooN
{
	using std::numeric_limits;

	/// brent_line_search performs Brent's golden section/quadratic interpolation search
	/// on the functor provided. The inputs a, x, b must bracket the minimum, and
	/// must be in order, so  that \f$ a < x < b \f$ and \f$ f(a) > f(x) < f(b) \f$.
	/// @param a The most negative point along the line.
	/// @param x The central point.
	/// @param fx The value of the function at the central point (\f$b\f$).
	/// @param b The most positive point along the line.
	/// @param func The functor to minimize
	/// @param maxiterations  Maximum number of iterations
	/// @param tolerance Tolerance at which the search should be stopped (defults to sqrt machine precision)
	/// @param epsilon Minimum bracket width (defaults to machine precision)
	/// @return The minima position is returned as the first element of the vector,
	///         and the minimal value as the second element.
	template<class Functor, class Precision> Vector<2, Precision> brent_line_search(Precision a, Precision x, Precision b, Precision fx, const Functor& func, int maxiterations, Precision tolerance = sqrt(numeric_limits<Precision>::epsilon()), Precision epsilon = numeric_limits<Precision>::epsilon())
	{
		using std::min;
		using std::max;

		using std::abs;
		using std::cout;
		using std::endl;

		Plotter plot;

		plot.s() << "set nokey\n";

		std::vector<Vector<2> > curve;
		double ymin=1e99, ymax=-1e99;

		for(double xx=a-.2; xx < b+.2; xx+= (b-a)/1000)
		{
			curve.push_back(makeVector(xx, func(xx)));
			ymin = min(curve.back()[1], ymin);
			ymax = max(curve.back()[1], ymax);
		}

		//The golden ratio:
		const Precision g = (3.0 - sqrt(5))/2;
		
		//The following points are tracked by the algorithm:
		//a, b bracket the interval
		// x   is the best value so far
		// w   second best point so far
		// v   third best point so far
		// These may not be unique.
		
		//The following points are used during iteration
		// u   the point currently being evaluated
		// xm   (a+b)/2
		
		//The updates are tracked as:
		//e is the distance moved last step, or current if golden section is used
		//d is the point moved in the current step
		
		Precision w=x, v=x, fw=fx, fv=fx;
		
		Precision d=0, e=0;
		int i=0;

		while(abs(b-a) > (abs(a) + abs(b)) * tolerance + epsilon && i < maxiterations)
		{
			cout << "Starting iteration " << i << endl;
			
			//Plot the line and the brackets
			plot.newline("line lt 1").addpt(curve);
			plot.newline("line lt 1").addpt(a, ymin).addpt(a, ymax).skip().addpt(b, ymin).addpt(b,ymax);
			plot.newline("line lt 2").addpt(v, ymin).addpt(v,ymax).newline("points lt 2").addpt(v, fv); 
			plot.newline("line lt 3").addpt(w, ymin).addpt(w,ymax).newline("points lt 3").addpt(w, fw); 
			plot.newline("line lt 4").addpt(x, ymin).addpt(x,ymax).newline("points lt 4").addpt(x, fx); 

			i++;
			//The midpoint of the bracket
			const Precision xm = (a+b)/2;

			//Per-iteration tolerance 
			const Precision tol1 = abs(x)*tolerance + epsilon;

			//If we recently had an unhelpful step, then do
			//not attempt a parabolic fit. This prevents bad parabolic
			//fits spoiling the convergence. Also, do not attempt to fit if
			//there is not yet enough unique information in x, w, v.
			if(abs(e) > tol1 && w != v&& 0)
			{
				cout << "  Attempting parabolic fit\n";
				//Attempt a parabolic through the best 3 points. Imagine
				//constructing a Vandermonde matrix for polynomial fitting,
				//shifted so that x = 0, and f(x) = 0
				// xw = w-x
				// xv = v-x
				//     [  1  0   0    ]
				// V = [  1  xw  xw^2 ]
				//     [  1  xv  xv^2 ]
				//
				// det(V) = xw*xv^2-xw^2*xv
				//
				// The polynomial coefficients will be:
				//               [    0    ]
				//  C = inv(V) * [ fw - fx ]
				//               [ fv - fx ]
				//
				// If the polynmial is y=c_1 x^2 + c_2 x + c_3, then the minimum is at -c_2/2c_3
				// Which is clearly independent of the determenant. The minimum is at -c'_2/2c'_3
				// where C' = inv(V) * |V| * [ 0  fw-fx fv-fx]'
				
				const Precision fxw = fw - fx;
				const Precision fxv = fv - fx;
				const Precision xw = w-x;
				const Precision xv = v-x;

				const Precision c1 = fxw*xv-fxv*xw;
				const Precision c2 = -fxw*xv*xv+fxv*xw*xw;

				cout << "   fit parameters: " << c1  << " " << c2 << endl;

				//d is the distance offset
				//d = -0.5 * c2 / c1;

				const Precision newd = -0.5 * c2 / c1;

				if(c1 == 0 || abs(newd) > abs(e)*2 || x + newd > b || x+newd < a)
				{
					//Parabolic fit no good. Take a golden section step instead
					//and reset d and e.
					if(x > xm)
						e = a-x;
					else
						e = b-x;

					d = g*e;
				}
				else
				{
					//Parabolic fit was good. Shift d and e
					e = d;
					d = newd;
				}
			}
			else
			{
				cout << "  Going for gold\n";
				//Don't attempt a parabolic fit. Take a golden section step
				//instead and reset d and e.
				if(x > xm)
					e = a-x;
				else
					e = b-x;

				d = g*e;
			}

			const Precision u = x+d;
			//Our one function evaluation per iteration
			const Precision fu = func(u);

			if(fu < fx)
			{
				//U is the best known point.

				//Update the bracket
				if(u > x)
					a = x;
				else
					b = x;

				//Shift v, w, x
				v=w; fv = fw;
				w=x; fw = fx;
				x=u; fx = fv;
			}
			else
			{
				//u is not the best known point. However, it is within the
				//bracket.
				if(u < x)
					a = u;
				else
					b = u;

				if(fu <= fw || w == x)
				{
					//Here, u is the new second-best point
					v = w; fv = fw;
					w = u; fw = fu;
				}
				else if(fu <= fv || v==x || v == w)
				{
					//Here, u is the new third-best point.
					v = u; fv = fu;
				}
			}
			
			cout << "Iteration end: " << a << " " << b << endl;

			plot.draw();
			std::cin.get();
		}

		return makeVector(x, fx);
	}
}
#endif
