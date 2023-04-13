/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
// my approximate solution to B&S	pde with BlackScholesFromLogArgLessExp boundary
//#include "stdafx.h"
#include "SharedFunctionsAndTypes.h"
#include "CubaturePDESolver.h"
#include "ExactPDESolution.h"
#include "CloudWeightLocation.h"
#include "DebugTools.h"
#include "GeometricBM.h"
#include <iostream>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/hermite.hpp>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// black scholes solver AtExpCall
///*double CubaturePDESolver(double x, double strike, double barrier, double t, double r, double sigma,
//						 const CVectorWeightLocation & vlwNormalApproxMean0Var1,
//						 CMainProgramParameters& ipIntParams, double e)
//						 */
//double CubaturePDESolver(double x, double strike, double barrier, double t, double r, double sigma,
//						 const CVectorWeightLocation & vlwNormalApproxMean0Var1,
//						 CMainProgramParameters& ipIntParams, double e, unsigned int & DCount)
//{
//	//unsigned int DCount(0);
//
//	DCount = 0;
//
//	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
//		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
//		dProportionOfRemainingTimeConsumedThisStep, t, x,strike, barrier, r, sigma, vlwNormalApproxMean0Var1))
//		: (e != 0)
//		? e
//		: ipIntParams.tolerance;
//
//	CSetWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation arg;
//	arg.displacement = x;
//	arg.probability = 1;
//	slwLocationsWeightsLeft.insert(arg);
//	double dTimeLeft(t);
//	double dAccumulatedSolution(0);
//	double df = exp(-r*dTimeLeft);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		double drift  = (r-.5*sigma*sigma)*dTimeLeft;
//		double drift1 = (r-.5*sigma*sigma)*dTimeLeft*ipIntParams.dProportionOfRemainingTimeConsumedThisStep;
//		double drift2 = (r-.5*sigma*sigma)*dTimeLeft*(1-ipIntParams.dProportionOfRemainingTimeConsumedThisStep);
//		double sd	= sigma* sqrt(dTimeLeft);
//		double sd1	= sigma* sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd2	= sigma* sqrt(dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
//
//		CCloudWeightLocation slwCarryForward;
//
//		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			itExistingLocation != slwLocationsWeightsLeft.end();
//			++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//				itOuterNormal != vlwNormalApproxMean0Var1.end();
//				++itOuterNormal)
//			{
//				double temp = df * BlackScholesFromAtExpCall( GeometricBM(itExistingLocation->displacement, drift, sd*(*itOuterNormal).displacement), strike, barrier);
//				ans += (*itOuterNormal).probability * temp;	
//
//				for (
//					CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
//					itInnerNormal != vlwNormalApproxMean0Var1.end();
//				++itInnerNormal
//					)
//				{
//					ipIntParams.iCountBoundaryFunctionEvaluations ++;
//					/*ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
//						BlackScholesFromAtExpCall( 
//						GeometricBM(itExistingLocation->displacement,drift1, sd1*(*itOuterNormal).displacement) *
//						GeometricBM(1, drift2, sd2*(*itInnerNormal).displacement)
//						,strike,barrier) * df;*/	// gbm version
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
//						BlackScholesFromAtExpCall( 
//						GeometricBM(itExistingLocation->displacement,drift1, sd1*(*itOuterNormal).displacement) +
//						GeometricBM(0, drift2, sd2*(*itInnerNormal).displacement)
//						,strike,barrier) * df;
//				}
//			}
//			double dError;
//			if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//					itOuterNormal != vlwNormalApproxMean0Var1.end();
//					++itOuterNormal)
//				{
//					CWeightLocation lwNextWeightLocation;
//					lwNextWeightLocation.displacement = GeometricBM(itExistingLocation->displacement,drift1, sd1*(*itOuterNormal).displacement );
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					slwCarryForward.insert(lwNextWeightLocation);
//				}
//			}
//		}
//		
//		if (ipIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//		double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
//			rbegin()->displacement);
//#if 1		
//		
//		slwCarryForward.ReducePointSet
//			(dAccuracySpread * ipIntParams.dAdjustClusterDiameter
//			, (unsigned int)(ipIntParams.INumerator() * vlwNormalApproxMean0Var1
//			.size()) / ipIntParams.IDenominator()
//			, ipIntParams);
//
//		std::cout << slwCarryForward.size() << std::endl;
//
//#else		
//		///// changed for test/////////////////////////////
//		
//		// patch determination
//
//		double patchEpsilon(1.0e-10);
//		int iNoCubatureToDimension = (ipIntParams.INumerator() * vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
//		
//		int n = iNoCubatureToDimension;	//rename to coincide with experiment
//		double z(barrier/std::sqrt(dTimeLeft));
//		double y(0.8*z);
//		double dAtExpCallnPlusOne = 1./std::sqrt(boost::math::constants::pi<double>()) * std::pow(1./std::sqrt(2.),n) * 
//			(  std::exp(-(y-z)*(y-z)/2.) * ( z/std::sqrt(2.) * boost::math::hermite(n,(y-z)/std::sqrt(2.)) * std::pow(-1.,n+1)  -   boost::math::hermite(n-1, (z-y)/std::sqrt(2.))  )
//												+  boost::math::hermite(n-1,y/std::sqrt(2.))* std::pow(-1.,n-1) * std::exp(-y*y/2.)	);
//
//		double spread = std::pow(patchEpsilon * boost::math::factorial<double>(n) / std::sqrt(dTimeLeft) / std::abs(dAtExpCallnPlusOne), 1./n);
//		slwCarryForward.ReducePointSet( spread, iNoCubatureToDimension, ipIntParams);
//		
//		std::cout << dAccuracySpread << "," << dAtExpCallnPlusOne << "," << spread << "," << slwCarryForward.size() << std::endl;
//#endif		
//		///       end test ///////////////////////////////
//
//		slwLocationsWeightsLeft.swap(slwCarryForward);
//		slwCarryForward.clear();
//		DCount++;
//	}
//
//	return dAccumulatedSolution; 
//
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// black scholes solver CALL
//double CubaturePDESolver(double x, double strike, double t, double r, double sigma,
//						 const CVectorWeightLocation & vlwNormalApproxMean0Var1,
//						 CMainProgramParameters& ipIntParams, double e)
//{
//	unsigned int DCount(0);
//
//	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
//		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
//		dProportionOfRemainingTimeConsumedThisStep, t, exp(x), strike, r, sigma, vlwNormalApproxMean0Var1))
//		: (e != 0)
//		? e
//		: ipIntParams.tolerance;
//	
//	CSetWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation arg;
//	arg.displacement = x;
//	arg.probability = 1;
//	slwLocationsWeightsLeft.insert(arg);
//	double dTimeLeft(t);
//	double dAccumulatedSolution(0);
//	double df = exp(-r*dTimeLeft);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		double drift  = (r-.5*sigma*sigma)*dTimeLeft;
//		double drift1 = (r-.5*sigma*sigma)*dTimeLeft*ipIntParams.dProportionOfRemainingTimeConsumedThisStep;
//		double drift2 = (r-.5*sigma*sigma)*dTimeLeft*(1-ipIntParams.dProportionOfRemainingTimeConsumedThisStep);
//		double sd	= sigma* sqrt(dTimeLeft);
//		double sd1	= sigma* sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd2	= sigma* sqrt(dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
//		
//		CCloudWeightLocation slwCarryForward;
//
//		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			 itExistingLocation != slwLocationsWeightsLeft.end();
//			 ++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//				 itOuterNormal != vlwNormalApproxMean0Var1.end();
//				 ++itOuterNormal)
//			{
// 				// double temp = df * BlackScholesFromCall( GeometricBM(itExistingLocation->displacement, drift, sd*(*itOuterNormal).displacement), strike);
//				double temp = df * BlackScholesFromCall(exp(drift + itExistingLocation->displacement + sd*(*itOuterNormal).displacement), strike);
//				ans += (*itOuterNormal).probability * temp;	
//
//				for (
//					CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
//					itInnerNormal != vlwNormalApproxMean0Var1.end();
//					++itInnerNormal
//					)
//				{
//					ipIntParams.iCountBoundaryFunctionEvaluations ++;
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
//						BlackScholesFromCall( exp(drift1 + itExistingLocation->displacement + sd1 * (*itOuterNormal).
//						displacement + drift2 + sd2 * (*itInnerNormal).displacement)
//						//GeometricBM(itExistingLocation->displacement,drift1, sd1*(*itOuterNormal).displacement) *
//						//GeometricBM(1., drift2, sd2*(*itInnerNormal).displacement)  // was GeometricBM(1, drift2, sd2*(*itInnerNormal).displacement)
//						,strike) * df;
//				}
//			}
//			double dError;
//			if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//					 itOuterNormal != vlwNormalApproxMean0Var1.end();
//					 ++itOuterNormal)
//				{
//					CWeightLocation lwNextWeightLocation;
//					lwNextWeightLocation.displacement = drift1 + itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement;
//						//GeometricBM(itExistingLocation->displacement,drift1, sd1*(*itOuterNormal).displacement );
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					slwCarryForward.insert(lwNextWeightLocation);
//				}
//			}
//		}
//		double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
//										rbegin()->displacement);
//
//		if (ipIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//		slwCarryForward.ReducePointSet
//			(dAccuracySpread * ipIntParams.dAdjustClusterDiameter
//			, (unsigned int)(ipIntParams.INumerator() * vlwNormalApproxMean0Var1
//							 .size()) / ipIntParams.IDenominator()
//			, ipIntParams);
//
//		slwLocationsWeightsLeft.swap(slwCarryForward);
//		slwCarryForward.clear();
//		DCount++;
//	}
//
//	return dAccumulatedSolution;
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bm with drift
double CubaturePDESolver(double x, double t, double r, double sigma, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
						 CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount )
{
	//CDebugTools& dbgLog = ipIntParams.dbgLogFile;


#ifdef TJLVERBOSE
	//	dbgLog.m_out << "TJLVERBOSE defined in libs, Producing Debug Data" << std::endl;
#endif

	//unsigned int DCount(0);

	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
		dProportionOfRemainingTimeConsumedThisStep, t, x, r, sigma, vlwNormalApproxMean0Var1))
		: (e != 0)
		? e
		: ipIntParams.tolerance;

	//e=1.4521717162097048e-012;

	//double patchEpsilon(10.e-10);
	
	CSetWeightLocation slwLocationsWeightsLeft;
	CWeightLocation arg;
	arg.displacement = x;
	arg.probability = 1;
	slwLocationsWeightsLeft.insert(arg);
	double dTimeLeft(t);
	double dAccumulatedSolution(0);

	/*double spreadConst = std::pow(patchEpsilon * 
								(double)(iNoCubatureToDimension + 1) * 
								std::sqrt(2. * boost::math::constants::pi<double>()) * 
								boost::math::factorial<double>(iNoCubatureToDimension / 2) * 
								std::pow(2.,(double)(iNoCubatureToDimension/2))      ,    1./(double)(iNoCubatureToDimension+1));*/

	while (slwLocationsWeightsLeft.size() != 0)
	{

//		std::cout<< slwLocationsWeightsLeft.size() << std::endl;
		// drift
		double drift	= r * dTimeLeft;
		double drift1	= r * ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft;
		double drift2	= r * (1.-ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dTimeLeft; 
		// std dev
		double sd	= sigma * sqrt(dTimeLeft);
		double sd1	= sigma * sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
		double sd2	= sigma * sqrt(dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
			
		CCloudWeightLocation slwCarryForward;

		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
			 itExistingLocation != slwLocationsWeightsLeft.end();
			 ++itExistingLocation)
		{
			double ans(0), ans1(0);
			// compute both approximations
			for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
				 itOuterNormal != vlwNormalApproxMean0Var1.end();
				 ++itOuterNormal)
			{
				double temp = BlackScholesFromLogArgLessExp(drift+itExistingLocation->displacement+sd*(*itOuterNormal).displacement);
				ans += (*itOuterNormal).probability*temp;	

				for (
					CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
					itInnerNormal != vlwNormalApproxMean0Var1.end();
					++itInnerNormal
					)
				{
					ipIntParams.iCountBoundaryFunctionEvaluations ++;
					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
						BlackScholesFromLogArgLessExp(drift1 + itExistingLocation->displacement + sd1 * (*itOuterNormal).
						displacement + drift2 + sd2 * (*itInnerNormal).displacement); // path concatenation
				}
			}
			double dError;
			if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
			{
				dAccumulatedSolution += itExistingLocation->probability * ans1;
			}
			else
			{
				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
					 itOuterNormal != vlwNormalApproxMean0Var1.end();
					 ++itOuterNormal)
				{
					CWeightLocation lwNextWeightLocation;
					lwNextWeightLocation.displacement = drift1 + itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement;
					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
					slwCarryForward.insert(lwNextWeightLocation);
				}
			}
		}
		double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
										rbegin()->displacement);

		if (ipIntParams.bVerbose)
			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;

#if 0
		double spread = spreadConst * std::pow(dTimeLeft,0.5);
		slwCarryForward.ReducePointSet
			( spread
			, iNoCubatureToDimension
			, ipIntParams);
		std::cout << spread << std::endl;
#else
		unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int) vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
		double spread = dAccuracySpread * ipIntParams.dAdjustClusterDiameter;	
		slwCarryForward.ReducePointSet //Multithreaded
			( spread
			, iNoCubatureToDimension
			, ipIntParams);
#endif
		slwLocationsWeightsLeft.swap(slwCarryForward);
		slwCarryForward.clear();
		DCount++;		
	}

	return dAccumulatedSolution;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//bm
double CubaturePDESolver(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
						 CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount )
{
	return CubaturePDESolver(x,t, 0., 1., vlwNormalApproxMean0Var1, ipIntParams, e, DCount);
}


// heat equation with adaptive approximation
//template <typename F>
//double CubaturePDESolverWithAdaptiveApprox(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
//	CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int& DCount, F & funcwrapper)
//{
//	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
//		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
//			dProportionOfRemainingTimeConsumedThisStep, t, x, vlwNormalApproxMean0Var1))
//		: (e != 0)
//		? e
//		: ipIntParams.tolerance;
//
//	//LogArgLessExp myfunc;
//	//FunctionData data;
//	//ApproximationMorton<LogArgLessExp, FunctionData > funcwrapper(myfunc, 6, std::pow(10., -12));	// wrapper function algorithm
//
//	CSetWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation arg;
//	arg.displacement = x;
//	arg.probability = 1;
//	slwLocationsWeightsLeft.insert(arg);
//	double dTimeLeft(t);
//	double dAccumulatedSolution(0);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		double sd = sqrt(dTimeLeft);
//		double sd1 = sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd2 = sqrt(dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
//
//		CCloudWeightLocation slwCarryForward;
//
//		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			itExistingLocation != slwLocationsWeightsLeft.end();
//			++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//				itOuterNormal != vlwNormalApproxMean0Var1.end();
//				++itOuterNormal)
//			{
//				point pt_ = { itExistingLocation->displacement + sd * (*itOuterNormal).displacement };
//				double temp = funcwrapper(pt_, data);
//				ans += (*itOuterNormal).probability * temp;
//
//				for (
//					CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
//					itInnerNormal != vlwNormalApproxMean0Var1.end();
//					++itInnerNormal
//					)
//				{
//					ipIntParams.iCountBoundaryFunctionEvaluations++;
//					point pt_outer = { itExistingLocation->displacement
//							+ sd1 * (*itOuterNormal).displacement
//							+ sd2 * (*itInnerNormal).displacement };
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
//						funcwrapper(pt_outer, data); // path concatenation
//				}
//			}
//			double dError;
//			if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//					itOuterNormal != vlwNormalApproxMean0Var1.end();
//					++itOuterNormal)
//				{
//					CWeightLocation lwNextWeightLocation;
//					lwNextWeightLocation.displacement = itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement;
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					slwCarryForward.insert(lwNextWeightLocation);
//				}
//			}
//		}
//		double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
//			rbegin()->displacement);
//
//		if (ipIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//		unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int)vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
//		double spread = dAccuracySpread * ipIntParams.dAdjustClusterDiameter;
//		slwCarryForward.ReducePointSet
//		(spread
//			, iNoCubatureToDimension
//			, ipIntParams);
//
//		slwLocationsWeightsLeft.swap(slwCarryForward);
//		slwCarryForward.clear();
//		DCount++;
//	}
//
//	return dAccumulatedSolution;
//}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
double CubaturePDESolverThreeStepOneStep(double x, double t, double r, double sigma, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
											 CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount )
{
	//CDebugTools& dbgLog = ipIntParams.dbgLogFile;

#ifdef TJLVERBOSE
	//	dbgLog.m_out << "TJLVERBOSE defined in libs, Producing Debug Data" << std::endl;
#endif

	//unsigned int DCount(0);

	/*
	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
		dProportionOfRemainingTimeConsumedThisStep, t, x, r, sigma, vlwNormalApproxMean0Var1))
		: (e != 0)
		? e
		: ipIntParams.tolerance;
	*/

	CSetWeightLocation slwLocationsWeightsLeft;
	CWeightLocation arg;
	arg.displacement = x;
	arg.probability = 1;
	slwLocationsWeightsLeft.insert(arg);
	double dTimeLeft(t);
	double dAccumulatedSolution(0);

	while (slwLocationsWeightsLeft.size() != 0)
	{
		double dAfterStepOneTimeLeft = dTimeLeft * (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep);
		// drift
		double stepOneDrift = r * dTimeLeft * ipIntParams.dProportionOfRemainingTimeConsumedThisStep ;
		double stepTwoDrift	= r * dAfterStepOneTimeLeft;
		double stepTwoDrift1	= r * ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dAfterStepOneTimeLeft;
		double stepTwoDrift2	= r * (1.-ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dAfterStepOneTimeLeft; 
		// std dev
		double stepOneSD = sigma * sqrt (dTimeLeft * ipIntParams.dProportionOfRemainingTimeConsumedThisStep );
		double stepTwoSD	= sigma * sqrt(dAfterStepOneTimeLeft);
		double stepTwoSD1	= sigma * sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dAfterStepOneTimeLeft);
		double stepTwoSD2	= sigma * sqrt(dAfterStepOneTimeLeft * (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
		dTimeLeft = dAfterStepOneTimeLeft;

		CSetWeightLocation slwStepOneLocationsWeights;
		CCloudWeightLocation slwCarryForward;

		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
			itExistingLocation != slwLocationsWeightsLeft.end();
			++itExistingLocation)
		{
			// go fwd one step
			for (CVectorWeightLocation::const_iterator itStepOneNormal(vlwNormalApproxMean0Var1.begin());
				itStepOneNormal != vlwNormalApproxMean0Var1.end();
				++itStepOneNormal)
			{
				CWeightLocation lwNextStepOneWeightLocation;
				lwNextStepOneWeightLocation.displacement = itExistingLocation->displacement + stepOneDrift + stepOneSD*(*itStepOneNormal).displacement;
				lwNextStepOneWeightLocation.probability = itExistingLocation->probability * (*itStepOneNormal).probability;
				slwStepOneLocationsWeights.insert(lwNextStepOneWeightLocation);
			}

			bool stepFwdFlag(0);
			double val(0);

			// for each intermediate one step...do one-step-two-step comparison
			for (CSetWeightLocation::iterator itExistingLocation = slwStepOneLocationsWeights.begin();
				itExistingLocation != slwStepOneLocationsWeights.end() && stepFwdFlag==0;
				++itExistingLocation)
			{
				double ans(0), ans1(0);
				// compute both approximations
				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
					itOuterNormal != vlwNormalApproxMean0Var1.end();
					++itOuterNormal)
				{
					double temp = BlackScholesFromLogArgLessExp(stepTwoDrift+itExistingLocation->displacement+stepTwoSD*(*itOuterNormal).displacement);
					ans += (*itOuterNormal).probability*temp;

					for (CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
						itInnerNormal != vlwNormalApproxMean0Var1.end();
						++itInnerNormal)
					{
						ipIntParams.iCountBoundaryFunctionEvaluations ++;
						ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
							BlackScholesFromLogArgLessExp(stepTwoDrift1 + itExistingLocation->displacement + stepTwoSD1 * (*itOuterNormal).
							displacement + stepTwoDrift2 + stepTwoSD2 * (*itInnerNormal).displacement); // path concatenation
					}
				}
				double dError;
				if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
					val += itExistingLocation->probability * ans1;
				else
					stepFwdFlag = 1;
			}	// end of two steps.

			if (stepFwdFlag)
			{
				for (CSetWeightLocation::iterator itExistingLocation = slwStepOneLocationsWeights.begin();
					itExistingLocation != slwStepOneLocationsWeights.end();
					++itExistingLocation)
				{
					CWeightLocation lwNextWeightLocation;
					lwNextWeightLocation.displacement = itExistingLocation->displacement;
					lwNextWeightLocation.probability = itExistingLocation->probability;
					slwCarryForward.insert(lwNextWeightLocation);
				}
			} else
				dAccumulatedSolution+= val;	// don't need to further expand so add to solution
			
			slwStepOneLocationsWeights.clear();
		
		} // end for each to be expanded

		if (ipIntParams.bVerbose)
			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;

		unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int) vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();	

		//std::cout<< "number of particles is " << slwLocationsWeightsLeft.size() << std::endl;
#if 1
		// patch determination
		double patchEpsilon(1.0e-10);
		int n = 1+iNoCubatureToDimension;
		double spreadConst = std::sqrt(2.) * std::pow(patchEpsilon * std::sqrt(boost::math::constants::pi<double>()) * (double)(n) * boost::math::factorial<double>((n-1)/2) ,1./(double)(n));
		double spread = spreadConst * std::pow(dTimeLeft,0.5);

		//spread = std::pow(dTimeLeft,0.5); // test

		slwCarryForward.ReducePointSet 
			( spread
			, iNoCubatureToDimension
			, ipIntParams);

#else
		double dAccuracySpread = stepOneSD * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.rbegin()->displacement);
		double spread = dAccuracySpread * ipIntParams.dAdjustClusterDiameter;	
		slwCarryForward.ReducePointSet
			( spread
			, iNoCubatureToDimension
			, ipIntParams);

#endif

		//std::cout << "spreadConst is " << spreadConst << std::endl;
		//std::cout<< "spread is " << spread << std::endl;

		slwLocationsWeightsLeft.swap(slwCarryForward);
		slwCarryForward.clear();
		DCount++;		

	}	// end while loop

	return dAccumulatedSolution;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CubaturePDESolverThreeStepOneStep(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
											 CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount )
{
	return CubaturePDESolverThreeStepOneStep(x,t, 0., 1., vlwNormalApproxMean0Var1, ipIntParams, e, DCount);
}

// a 15 point cubature produces 
//the value of ExactPDESolution  at (0, 1) is 0.238422
//CubaturePDESolver was evaluated 22291 times
//with 15 points the difference between cubature and answer is -1.11022e-016
//with 15 points the difference between cubature and answer is -2.77556e-017
//with 15 points the difference between cubature and answer is 0
//with 15 points the difference between cubature and answer is -2.77556e-017
//with 15 points the difference between cubature and answer is -2.77556e-017
//with 15 points the difference between cubature and answer is -5.55112e-017
//with 15 points the difference between cubature and answer is -1.52656e-016
//with 15 points the difference between cubature and answer is 1.38778e-017
//with 15 points the difference between cubature and answer is 0
//Timer factor was 0.8
//with 15 points the difference between cubature and answer is -1.61013e-012
//with 15 points the difference between cubature and answer is -1.38131e-012
//with 15 points the difference between cubature and answer is -1.16709e-012
//with 15 points the difference between cubature and answer is -9.66449e-013
//with 15 points the difference between cubature and answer is -7.78461e-013
//with 15 points the difference between cubature and answer is -6.02185e-013
//with 15 points the difference between cubature and answer is -4.37261e-013
//with 15 points the difference between cubature and answer is -2.83357e-013
//with 15 points the difference between cubature and answer is -1.41484e-013
//Timer factor was 0.64
//with 15 points the difference between cubature and answer is -5.13909e-010
//with 15 points the difference between cubature and answer is -4.38583e-010
//with 15 points the difference between cubature and answer is -3.6805e-010
//with 15 points the difference between cubature and answer is -3.0196e-010
//with 15 points the difference between cubature and answer is -2.39998e-010
//with 15 points the difference between cubature and answer is -1.81935e-010
//with 15 points the difference between cubature and answer is -1.2766e-010
//with 15 points the difference between cubature and answer is -7.72365e-011
//with 15 points the difference between cubature and answer is -3.09422e-011
//Timer factor was 0.512
//with 15 points the difference between cubature and answer is -2.52465e-008
//with 15 points the difference between cubature and answer is -2.14399e-008
//with 15 points the difference between cubature and answer is -1.78748e-008
//with 15 points the difference between cubature and answer is -1.45334e-008
//with 15 points the difference between cubature and answer is -1.14008e-008
//with 15 points the difference between cubature and answer is -8.46677e-009
//with 15 points the difference between cubature and answer is -5.7283e-009
//with 15 points the difference between cubature and answer is -3.1919e-009
//with 15 points the difference between cubature and answer is -8.76071e-010
//Timer factor was 0.4096
//with 15 points the difference between cubature and answer is -4.29416e-007
//with 15 points the difference between cubature and answer is -3.63011e-007
//with 15 points the difference between cubature and answer is -3.00807e-007
//with 15 points the difference between cubature and answer is -2.42498e-007
//with 15 points the difference between cubature and answer is -1.87836e-007
//with 15 points the difference between cubature and answer is -1.36673e-007
//with 15 points the difference between cubature and answer is -8.89964e-008
//with 15 points the difference between cubature and answer is -4.49763e-008
//with 15 points the difference between cubature and answer is -5.00755e-009
//Timer factor was 0.32768
//with 15 points the difference between cubature and answer is -3.65915e-006
//with 15 points the difference between cubature and answer is -3.08056e-006
//with 15 points the difference between cubature and answer is -2.53847e-006
//with 15 points the difference between cubature and answer is -2.03028e-006
//with 15 points the difference between cubature and answer is -1.55396e-006
//with 15 points the difference between cubature and answer is -1.10844e-006
//with 15 points the difference between cubature and answer is -6.93956e-007
//with 15 points the difference between cubature and answer is -3.12441e-007
//with 15 points the difference between cubature and answer is 3.21011e-008
//Timer factor was 0.262144
//with 15 points the difference between cubature and answer is -1.91581e-005
//with 15 points the difference between cubature and answer is -1.60696e-005
//with 15 points the difference between cubature and answer is -1.31757e-005
//with 15 points the difference between cubature and answer is -1.04625e-005
//with 15 points the difference between cubature and answer is -7.91999e-006
//with 15 points the difference between cubature and answer is -5.54367e-006
//with 15 points the difference between cubature and answer is -3.33637e-006
//with 15 points the difference between cubature and answer is -1.31056e-006
//with 15 points the difference between cubature and answer is 5.09691e-007
//Timer factor was 0.209715
//with 15 points the difference between cubature and answer is -7.00006e-005
//with 15 points the difference between cubature and answer is -5.85275e-005
//with 15 points the difference between cubature and answer is -4.77754e-005
//with 15 points the difference between cubature and answer is -3.7695e-005
//with 15 points the difference between cubature and answer is -2.82511e-005
//with 15 points the difference between cubature and answer is -1.94306e-005
//with 15 points the difference between cubature and answer is -1.12495e-005
//with 15 points the difference between cubature and answer is -3.7613e-006
//with 15 points the difference between cubature and answer is 2.93615e-006
//Timer factor was 0.167772
//with 15 points the difference between cubature and answer is -0.000194572
//with 15 points the difference between cubature and answer is -0.000162231
//with 15 points the difference between cubature and answer is -0.000131919
//with 15 points the difference between cubature and answer is -0.000103502
//with 15 points the difference between cubature and answer is -7.68851e-005
//with 15 points the difference between cubature and answer is -5.20418e-005
//with 15 points the difference between cubature and answer is -2.90305e-005
//with 15 points the difference between cubature and answer is -8.01837e-006
//with 15 points the difference between cubature and answer is 1.06975e-005
//Timer factor was 0.134218
//with 15 points the difference between cubature and answer is -0.000437618
//with 15 points the difference between cubature and answer is -0.000364016
//with 15 points the difference between cubature and answer is -0.00029503
//with 15 points the difference between cubature and answer is -0.000230354
//with 15 points the difference between cubature and answer is -0.000169792
//with 15 points the difference between cubature and answer is -0.000113298
//with 15 points the difference between cubature and answer is -6.10326e-005
//with 15 points the difference between cubature and answer is -1.34088e-005
//with 15 points the difference between cubature and answer is 2.88575e-005
//Timer factor was 0.107374
//Press any key to continue