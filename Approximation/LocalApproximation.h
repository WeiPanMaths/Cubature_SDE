#pragma once

#include "Utils.h"
#include <vector>
#include <iostream>
#include <boost\math\special_functions\chebyshev.hpp>

#if dim_ == 1
// interpolation on [-1, 1]
class InterpolationMorton
{
public:

	/// CONSTRUCTOR /////////////////////////////////////////////////////////////////////////////////////////////////

	/// polynomial fitting constructor
	InterpolationMorton(unsigned int totaldegree = 1)
		: myTotalDegree(totaldegree)
		, myNoOfMonomials( 1 + totaldegree )		// for 2D only, higher dimensional it's (d+m)!/d!/m!
		, myApproximationCoefficients(myNoOfMonomials)
	{}

	~InterpolationMorton(void) {}


	/// PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////////////////////////////


	double approximateFunction(const point & x) const
	{
		// Clenshaw's evaluation algorithm using Boost
		// ! note that in Boost's implementation, the first term in the Chebyshev sum = coeff[0]/2
		
		// point needs to be remapped to [-1,1], thus we need to know the interval, or patch
		try {
			if (x[0] < -1 || x[0] > 1)
				throw (x[0]); 
		}
		catch (int myNum) {
			std::cout << "approximateFunction: value outside chebyshev range\n"
				<< myNum;
			exit(-1);
		}
		double ans = boost::math::chebyshev_clenshaw_recurrence(myApproximationCoefficients.data(), 
									myApproximationCoefficients.size(), x[0]);
		return ans;
	}

	// cross validation method //////////
	// maybe better to use Cubature double comparison //
	bool isAccurate(const std::vector<point>& pts, const std::vector<double>& values, const double tolerance)
	{
		auto rmse = [&]()->double
		{
			double err(0);
			for (auto i = 0; i < pts.size(); ++i)
			{
				auto diff = values[i] - approximateFunction(pts[i]);
				err += diff * diff;
			}
			err /= (double)pts.size();

			return std::sqrt(err);
		};

		auto std_rel = [&](const double correction)->double
		{
			double err = 0;
			for (auto i = 0; i < pts.size(); ++i)
			{
				auto abserror = std::abs(values[i] - approximateFunction(pts[i]));
				auto relerror = abserror / std::abs(values[i] + correction);
				if (relerror > err) err = relerror;
			}
			return err;
		};

		auto std_abs = [&]()->double
		{
			double err = 0;
			for (auto i = 0; i < pts.size(); ++i)
			{
				auto abserror = std::abs(values[i] - approximateFunction(pts[i]));
				if (abserror > err) err = abserror;
			}
			return err;
		};
# if RMSE		
		auto max_error = rmse();
# else
		//auto max_error = std_rel(1.e-15);// 1.e-10);
		auto max_error = std_abs();
# endif		
		if (max_error < tolerance)	return true;
		else						return false;
	}

	/// mkl least square approximation function using svd
	bool interpolate_svd(const std::vector<point>& pts, const std::vector<double>& values);

	// static size_t  NumFuncCall;			// FOR TESTING PURPOSES
	std::vector<double> getCoeffs() const
	{
		return myApproximationCoefficients;
	}


private:

	/// MEMBER DATA ////////////////////////////////////////////////////////////////////////////////////////////////////

	BLOCK						myPoint;
	const unsigned int			myTotalDegree;
	const unsigned int			myNoOfMonomials;

	std::vector<double>			myApproximationCoefficients;		// stores the result

	/// PRIVATE METHODS	////////////////////////////////////////////////////////////////////////////////////////

	/// generates matrix of monomial values for a series of points 
	void generate_monomial_matrix(std::vector<double> & a, const std::vector<point> & x)
	{
		for (unsigned int i = 0; i < x.size(); ++i)
		{
			try {
				if (x[i][0] < -1 && x[i][0] > 1)
					throw (x[i][0]);
			}
			catch (int myNum) {
				std::cout << "generate_monomial_matrix: value outside chebyshev range\n"
					<< myNum;
			}

			double T0 = 1;
			double T1 = x[i][0];
			a[i] = 1;
			a[x.size()+i] = T1;
			for(unsigned l = 2; l < myNoOfMonomials; ++l)
			{
			    std::swap(T0, T1);
			    T1 = boost::math::chebyshev_next(x[i][0], T0, T1);
				a[l * x.size() + i] = T1;
			}
		}
	}

};
#endif