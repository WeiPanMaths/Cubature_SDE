#include "BlackScholesGreeks.h"
#include "CubaturePDESolver.h"
//#include "ConstantValues.h"
#include <boost/math/distributions/normal.hpp>

using boost::math::normal;

// five point
double BSDeltaFivePointCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon)
{
	double spot(ipIntParams.dLocation), strike(ipIntParams.dStrike), rate(ipIntParams.dRate), vol(ipIntParams.dVol), timeToMaturity(ipIntParams.dTime);	
	double soln1 = CubaturePDESolver(spot-2*epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	double soln2 = CubaturePDESolver(spot-epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	double soln3 = CubaturePDESolver(spot+epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	double soln4 = CubaturePDESolver(spot+2*epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);

	return (soln1 - soln2*8.+8.*soln3 - soln4) / (12.*epsilon);
}

// first order forward difference
double BSDeltaFwdDifferenceCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon)
{
	double spot(ipIntParams.dLocation), strike(ipIntParams.dStrike), rate(ipIntParams.dRate), vol(ipIntParams.dVol), timeToMaturity(ipIntParams.dTime);
	double solution = CubaturePDESolver(ipIntParams, vlwNormalApproxMean0Var1, 0.);
	double solutionFwd = CubaturePDESolver(spot+epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	return (solutionFwd - solution)/epsilon;
}

// central difference
double BSDeltaCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon)
{
#if 1
	double spot(ipIntParams.dLocation), strike(ipIntParams.dStrike), rate(ipIntParams.dRate), vol(ipIntParams.dVol), timeToMaturity(ipIntParams.dTime);	
	double solution = CubaturePDESolver(ipIntParams, vlwNormalApproxMean0Var1, 0.);
	double solutionMinusEpsilon = CubaturePDESolver(spot-epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	double solutionPlusEpsilon  = CubaturePDESolver(spot+epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	
	return .5*(solutionPlusEpsilon-solutionMinusEpsilon)/epsilon;
#else
	double spot(ipIntParams.dLocation), strike(ipIntParams.dStrike), rate(ipIntParams.dRate), vol(ipIntParams.dVol), timeToMaturity(ipIntParams.dTime);	
	double solution = CubaturePDESolver(ipIntParams, vlwNormalApproxMean0Var1, 0.);
	double solutionPlusTwoEpsilon = CubaturePDESolver(spot+2.*epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);
	double solutionPlusEpsilon  = CubaturePDESolver(spot+epsilon,strike,timeToMaturity,rate,vol,vlwNormalApproxMean0Var1,ipIntParams,0.);

	return .5*(4.*solutionPlusEpsilon - 3.*solution - solutionPlusTwoEpsilon)/epsilon;
#endif
}

double BSDeltaExactSolution(CMainProgramParameters& ipIntParams)
{
	double spot(ipIntParams.dLocation), strike(ipIntParams.dStrike), rate(ipIntParams.dRate), vol(ipIntParams.dVol), timeToMaturity(ipIntParams.dTime);
	normal n;
	double d1 = (log(spot/strike) + (rate + .5*vol*vol)*timeToMaturity)/(vol*sqrt(timeToMaturity));
	return cdf(n,d1);
}