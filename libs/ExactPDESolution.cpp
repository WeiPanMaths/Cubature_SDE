/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
//ExactPDESolution
#include "stdafx.h"
#include "ExactPDESolution.h"
#include "SharedFunctionsAndTypes.h"
#include <boost/math/special_functions/erf.hpp>
using boost::math::erfc;
using boost::math::erf;
#include <boost/math/distributions/normal.hpp>	// normal cdf


// the closed form solution to B&S pde with BlackScholesFromLogArgLessExp boundary
double ExactPDESolution(double y, double t)
{
	// boundary condition   max(1-exp(x), 0)
	return (erfc(y / (sqrt(2 * t))) - exp(t / 2. + y) * erfc((t + y) / (sqrt(2 * t)))) / 2.;
	//return ExactPDESolution(y,t,0.,1.);
}

// closed form solution
// boundary condition   max( 1 - exp(x) , 0 )
// dX=rdt+sigma dW
double ExactPDESolution( double x, double t, double r, double sigma )
{	 
#if	0
	/*return 0.5*(	erfc(  (r*t + x)/ (sigma*sqrt(2.*t)) ) 
			 -	exp( r*t + x + .5*sigma*sigma*t ) *  erfc( (r*t + x + sigma*sigma*t)/ (sigma*sqrt(2.*t))  ));*/

	return (erfc(x / (sqrt(2 * t))) - exp(t / 2. + x) * erfc((t + x) / (sqrt(2 * t)))) / 2.;	//terry's original
#else

	// binary
	return (1.-boost::math::cdf(boost::math::normal(x,sqrt(t)),double(0.)));

#endif
}

// closed form solution BS Call
// dS = rS dt + vol S dW
double ExactPDESolution( double spot, double strike, double timeToMat, double rate, double vol )
{
	return BlackScholesCall( spot, strike, rate, vol, timeToMat );
}

// closed form solution BS at expiry call option
double ExactPDESolution(double spot, double strike, double barrier, double timeToMat, double rate, double vol)
{
	return AtExpCall(spot, strike, barrier, rate, vol, timeToMat);
}