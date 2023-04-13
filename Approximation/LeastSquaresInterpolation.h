#pragma once

#include "Utils.h"
#include <vector>
#include <iostream>

#if dim_ == 2
//template <typename F>
class InterpolationMorton
{
public:	
	
	/// CONSTRUCTOR /////////////////////////////////////////////////////////////////////////////////////////////////

	/// polynomial fitting constructor
	InterpolationMorton(unsigned int totaldegree = 1)
		: myTotalDegree(totaldegree)
		, myNoOfMonomials((2+totaldegree)*(1+totaldegree)/2)		// for 2D only, higher dimensional it's (d+m)!/d!/m!
		, myApproximationCoefficients(myNoOfMonomials)	
	{}

	~InterpolationMorton(void){}

	
	/// PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////////////////////////////


	double approximateFunction(const point_2d & x) const
	{		
		std::vector<double> monomial_values(myNoOfMonomials);
		generate_monomial_serie(monomial_values, x);
		double ans(0);

		// dot product to produce the answer
		for(size_t i = 0; i<myApproximationCoefficients.size(); ++i)
			ans += monomial_values[i] * myApproximationCoefficients[i];;

		return ans;
	}

	// cross validation method //////////
	bool isAccurate( const std::vector<point_2d> & pts, const std::vector<double> & values, const double tolerance )
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
		auto max_error = std_rel(1.e-10);
# endif		
		if (max_error < tolerance)	return true;
		else						return false;
	}
	
	/// mkl least square approximation function using svd
	bool interpolate_svd(const std::vector<point_2d> & pts, const std::vector<double> & values);


	//bool interpolate( const std::vector<point_2d> & pts, const std::vector<double> & values )
	//{	
	//	int no_of_rows = (int)pts.size();
	//	int no_of_cols = myNoOfMonomials;

	//	int NRHS = 1, LDA = no_of_rows, LDB = no_of_rows, N = no_of_cols;
	//	char charN('N');
	//	int m = no_of_rows, n = no_of_cols, nrhs = NRHS, lda = LDA, ldb = LDB, info, lwork;
	//	double wkopt;
	//
	//	std::vector<double> work;
	//	std::vector<double> a(LDA*N);
	//	std::vector<double> b = values;
	//	
	//	// generate matrix A
	//	this->generate_monomial_matrix(a, pts);

	//	// find the optimal lwork size
	//	lwork = -1;
	//	dgels_( &charN, &m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &wkopt, &lwork, &info );
	//	lwork = (int)wkopt;
	//
	//	//using optimal lwork value
	//	work.resize(lwork);
	//	//std::cout << "using optimal lwork value" << std::endl;
	//	dgels_( &charN,&m,&n,  &nrhs, &a[0], &m, &b[0], &m, &work[0], &lwork, &info);

	//	if( info > 0 )
	//	{
	//		std::cout << "info > 0" << std::endl;

	//		std::cout <<"The diagonal element "<< info << " of the triangular factor "<<std::endl;
	//		std::cout << "of A is zero, so that A does not have full rank;\n";
	//		std::cout << "the least squares solution could not be computed.\n";
	//		return false;
	//	} else
	//	{
	//		// success - copy result into local variable
	//		//std::cout << "mkl dgels is successful" << std::endl;
	//		for(size_t i = 0; i<myNoOfMonomials; ++i)
	//			myApproximationCoefficients[i] = b[i];
	//		return true;
	//	}
	//}

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
	
	void generate_monomial_serie( std::vector<double> & ans, const point_2d & x) const
	{
		std::vector<double>::iterator a_it = ans.begin();
		for (unsigned int d = 0; d <= myTotalDegree; ++d)
			for (int i=d; i >= 0; --i )
			{	
				*a_it = std::pow(x[0], i) * std::pow(x[1], d-i);
				++a_it;
			}
	}

	/// generates matrix of monomial values for a series of points 
	void generate_monomial_matrix( std::vector<double> & a, const std::vector<point_2d> & x )
	{
		std::vector<double>::iterator a_it = a.begin();
	
		for (unsigned int d = 0; d <= myTotalDegree; ++d)
			for ( int i=d; i >= 0; --i )
				// to store matrix values in column major format we have this as the inner for loop
				for ( std::vector<point_2d>::const_iterator x_it = x.begin(); x_it != x.end(); ++x_it )
				{
					*a_it = std::pow((*x_it)[0], i)* std::pow( (*x_it)[1], d-i);
					++a_it;
				}
	}

};

//template <typename F>  size_t Interpolation<F>::NumFuncCall = 0;
#endif