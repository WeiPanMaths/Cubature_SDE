/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once

#ifndef SharedFunctionsAndTypes_h__
#define SharedFunctionsAndTypes_h__
#include "ExactPDESolution.h"
#include "CloudWeightLocation.h"
#include "GeometricBM.h"
#include <cmath>
#include <vector>
#include <set>
#include <valarray>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
//#include "AdaptiveApproximation.h"
#include "Utils.h" // point

using boost::math::erfc;
using boost::math::erf;
using boost::math::normal;


double ErrorInCubature(double dNumberOfCubaturePoints, double dTimeToBoundaryOfCubaturePoints,  double dTimeToCubaturePoints);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double TimeStepForAcceptableError(const double& dCurrentTime, const double& dAcceptableError, const int& degree);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double BlackScholesCall(double s, double K, double r, double sigma, double timeToMat)
{
	normal n;
	double fwd = s * exp( r* timeToMat);
	double d1 = ( log(fwd/K) + .5 * sigma*sigma*timeToMat) /  ( sigma * sqrt(timeToMat));
	double nd1 (0) , nd2 (0) ;
	try {
		nd1 = cdf(n,d1);
	}
	catch(...){
		std::cout<< "exception caught nd1, value"<< nd1 << std::endl;
	}

	try {
		nd2 = cdf(n, d1-sigma*sqrt(timeToMat));
	}
	catch(...){
		std::cout << "exception caught nd2, value " << nd2 << std::endl;
	}
	
	return (fwd*nd1 - K * nd2 )* exp( - r * timeToMat );	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double Digital(double s, double K, double r, double sigma, double timeToMat)
{
	normal n;
	double fwd = s*exp(r*timeToMat);
	return .5 * exp(-r*timeToMat) * erfc( -(log(fwd/K)-.5 *sigma*sigma*timeToMat) / ( sigma*sqrt(2.*timeToMat)) );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	At Expiry Call
inline double CompareCubatureAndTrueSolution(double dProportionTimeConsumed, double dTimeToFinalHorizon, double dSpatialDisplacement, double dStrike, double dBarrier, double rate, double vol, const CVectorWeightLocation& vlwNormalApproxMean0Var1)
{
	double ans(0);
	double sd1 = sqrt(dProportionTimeConsumed * dTimeToFinalHorizon);

	for (CVectorWeightLocation::const_iterator it1(vlwNormalApproxMean0Var1.begin());
		it1 != vlwNormalApproxMean0Var1.end();
		++it1)
	{
		/*ans += (*it1).probability * 
			ExactPDESolution(
			GeometricBM(dSpatialDisplacement, rate, vol, dProportionTimeConsumed*dTimeToFinalHorizon, sd1 * (*it1).displacement),
			dStrike,dBarrier, dTimeToFinalHorizon * (1 - dProportionTimeConsumed), rate, vol
			) * exp(-rate*dTimeToFinalHorizon*dProportionTimeConsumed*dTimeToFinalHorizon);*/
		ans += (*it1).probability * 
			ExactPDESolution(dSpatialDisplacement + sd1 * (*it1).displacement,
			dStrike,dBarrier, dTimeToFinalHorizon * (1. - dProportionTimeConsumed), rate, vol);		//changed to bm
		//ans += 0;
	}
	double temp = ExactPDESolution(dSpatialDisplacement, dStrike, dBarrier, dTimeToFinalHorizon, rate, vol);
	return ans -= temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Black Scholes Call
inline double CompareCubatureAndTrueSolution(double dProportionTimeConsumed, double dTimeToFinalHorizon,
											 double dSpatialDisplacement, double dStrike, double rate, double vol, const CVectorWeightLocation&
											 vlwNormalApproxMean0Var1)
{
	double ans(0);
	double sd1 = sqrt(dProportionTimeConsumed * dTimeToFinalHorizon);
	// compute both approximations
	for (CVectorWeightLocation::const_iterator it1(vlwNormalApproxMean0Var1.begin());
		it1 != vlwNormalApproxMean0Var1.end();
		++it1)
	{
		ans += (*it1).probability * 
			ExactPDESolution(
			GeometricBM(dSpatialDisplacement, rate, vol, dProportionTimeConsumed*dTimeToFinalHorizon, sd1 * (*it1).displacement),
			dStrike, dTimeToFinalHorizon * (1 - dProportionTimeConsumed), rate, vol
			) * exp(-rate*dTimeToFinalHorizon*dProportionTimeConsumed*dTimeToFinalHorizon);
	}
	double temp = ExactPDESolution(dSpatialDisplacement, dStrike, dTimeToFinalHorizon, rate, vol);
	return ans -= temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BM with drift
// Boundary condition max(1-e(x), 0)
inline double CompareCubatureAndTrueSolution(double dProportionTimeConsumed, double dTimeToFinalHorizon,
											 double dSpatialDisplacement, double rate, double vol, const CVectorWeightLocation&
											 vlwNormalApproxMean0Var1)
{
	double ans(0);
	double drift = rate * dProportionTimeConsumed * dTimeToFinalHorizon;
	double sd1 = vol * sqrt(dProportionTimeConsumed * dTimeToFinalHorizon);
	// compute both approximations
	for (CVectorWeightLocation::const_iterator it1(vlwNormalApproxMean0Var1.begin());
		it1 != vlwNormalApproxMean0Var1.end();
		++it1)
	{
		ans += (*it1).probability * ExactPDESolution(dSpatialDisplacement + drift + sd1 * (*it1).displacement,
			dTimeToFinalHorizon * (1 - dProportionTimeConsumed), rate, vol);
	}
	double temp = ExactPDESolution(dSpatialDisplacement, dTimeToFinalHorizon, rate, vol);
	return ans -= temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BM
// boundary conditoin max(1-e(x), 0)
inline double CompareCubatureAndTrueSolution(double dProportionTimeConsumed, double dTimeToFinalHorizon,
											 double dSpatialDisplacement, const CVectorWeightLocation&
											 vlwNormalApproxMean0Var1)
{
#if 0
	double ans(0);
	double sd1 = sqrt(dProportionTimeConsumed * dTimeToFinalHorizon);
	// compute both approximations
	for (CVectorWeightLocation::const_iterator it1(vlwNormalApproxMean0Var1.begin());
		it1 != vlwNormalApproxMean0Var1.end();
		++it1)
	{
		ans += (*it1).probability * ExactPDESolution(dSpatialDisplacement + sd1 * (*it1).displacement,
			dTimeToFinalHorizon * (1 - dProportionTimeConsumed));
	}
	return ans -= ExactPDESolution(dSpatialDisplacement, dTimeToFinalHorizon);
#endif
	return CompareCubatureAndTrueSolution( dProportionTimeConsumed, dTimeToFinalHorizon, dSpatialDisplacement, 0., 1., vlwNormalApproxMean0Var1);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// terry's original terminal condition
inline double BlackScholesFromLogArgLessExp(const double arg)
{
#if 1
	return (arg < double(0)) ? double(1) - exp(arg) : double(0);
#else
	return (arg < double(0.)) ? double(0) : double(1);	//indicator function
#endif
}
struct FunctionData {};

struct LogArgLessExp
{
	double operator() (const point& pt) const
	{
		return BlackScholesFromLogArgLessExp(pt[0]);
	}

	double operator() (double pt) const
	{
		return BlackScholesFromLogArgLessExp(pt);
	}

	double operator() (const point& pt, FunctionData& data) const
	{
		return BlackScholesFromLogArgLessExp(pt[0]);
	}
};

struct BenchmarkExact
{
	double operator() (const point& pt) const
	{
		return std::max( ExactPDESolution(pt[0], 1) - 0.05, 0. );
	}

	double operator() (double pt) const
	{
		return std::max(ExactPDESolution(pt, 1) - 0.05, 0.);
	}

	double operator() (const point& pt, FunctionData& data) const
	{
		return std::max(ExactPDESolution(pt[0], 1) - 0.05, 0.);
	}
};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double BlackScholesFromCall(double s, double k)
{
	return (s-k < double(0))? double(0) : s - k;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double BlackScholesFromAtExpCall(double s, double k, double b)
{
#if 1
	return ((s > k) && (s<b+k))? (s-k) : double(0);		// at exp payoff
#else
	return (s>b+k)? double(1) : double(0);	// digital payoff. Find out the sizes of the patch and interval.
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double AtExpCall(double s, double K, double B, double r, double sigma, double timeToMat)
{
	normal n;
	double d1 = ( log((B+K)/s) - ( r + .5*sigma*sigma)*timeToMat )  /  ( sigma * sqrt(timeToMat));
	double d2 = ( log((B+K)/s) - ( r - .5*sigma*sigma)*timeToMat )  /  ( sigma * sqrt(timeToMat));
	double nd1 (0) , nd2 (0) ;
	
	nd1 = cdf(n,d1);
	nd2 = cdf(n,d2);
	
	return (0 < B)? s*nd1 - K*nd2*exp(-r*timeToMat) : 0.;	
}


#endif //__SharedFunctionsAndTypes__




//template <class T>  std::ostream& operator << (std::ostream & os, const std::vector<T> & in)
//{
//	os << "{" ; 
//	for (typename std::vector<T>::const_iterator it(in.begin()); it != in.end() ; ++it)
//	{
//		os << *it;
//		if (it+1 != in.end()) os << ", "; 
//	}
//	os << "}";
//	return os;
//}