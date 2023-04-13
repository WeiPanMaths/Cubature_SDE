/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
// my approximate solution to B&S	pde with BlackScholesFromLogArgLessExp boundary
//#include "stdafx.h"
//#include "SharedFunctionsAndTypes.h"
#include "CubaturePDESolver2D.h"
//#include "ExactPDESolution.h"
#include "CloudWeightLocation2D.h"
//#include "DebugTools.h"
//#include "GeometricBM.h"
#include <iostream>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <random>
#include <ql/errors.hpp>
//#include "Utils.h"
//#include <wpUtilities.h>
//#include <float.h>


/// exact solution of E(X1(t) * X2(t)) for a 2D geometric BM (X1(t), X2(t)), with correlation rho
double ExactSolution(double x1, double x2, double r, double sigma1, double sigma2, double rho, double t)
{
	//auto a = (r - 0.5 * sigma1 * sigma1);
	//auto b = (r - 0.5 * sigma2 * sigma2);
	//auto ans = (x1 + a * t) * (x2 + b * t) + sigma1 * sigma2 * rho * t;
	//return ans * std::exp(-r*t);
	return std::exp(x1 + x2 + (2. * r + sigma1 * sigma2 * rho) * t) * std::exp(-r * t);
}

void QuadratureDoublePrecisionDeg5(CWeightLocation2D::CVectorWeightLocation& vlwNormalApproxMean0Var1)
{
	vlwNormalApproxMean0Var1.clear();
	std::vector<wpUtilities::Point<dim__>> vlwGaussianQuadratureDeg5Locations;
	//std::vector<wpUtilities2D::Point> vlwGaussianQuadratureDeg5Locations;
	std::vector<double> vlwGaussianQuadratureDeg5Weights;

	/*vlwGaussianQuadratureDeg5Locations = { {0., 0.}, {1., std::sqrt(3.)}, {-1., std::sqrt(3.)}, {1., -std::sqrt(3.)},
							{-1., -std::sqrt(3.)}, {1., std::sqrt(3.)}, {-1., std::sqrt(3.)},  {1., -std::sqrt(3.)},
							{-1., -std::sqrt(3.)},{2., 0.}, {-2., 0.}, {2., 0.}, {-2., 0.} };

	vlwGaussianQuadratureDeg5Weights = { 0.5, 1./24., 1./24., 1./24., 1./24., 1./24., 1./24.,
											  1./24., 1./24., 1./24., 1./24., 1./24., 1./24. };*/

	vlwGaussianQuadratureDeg5Locations = { {0., 0.}, {1., std::sqrt(3.)}, {-1., std::sqrt(3.)}, {1., -std::sqrt(3.)},
			{-1., -std::sqrt(3.)}, {2., 0.},  {-2., 0.} };  
	vlwGaussianQuadratureDeg5Weights = { 0.5, 1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12.};

	for (int i = 0; i < vlwGaussianQuadratureDeg5Locations.size(); ++i)
	{
		CWeightLocation2D::CWeightLocation arg;
		arg.displacement = vlwGaussianQuadratureDeg5Locations[i];
		arg.probability = vlwGaussianQuadratureDeg5Weights[i];
		vlwNormalApproxMean0Var1.push_back(arg);
	}
}

// compute 1 step cubature spread
void ComputeCubatureSpread(CWeightLocation2D::CVectorWeightLocation& vlwNormalApproxMean0Var1, const double rho, wpUtilities::Point<dim__> & spread)
{
	//Point max = vlwNormalApproxMean0Var1[0].displacement;
	//Point min = vlwNormalApproxMean0Var1[0].displacement;
	wpUtilities::Point<dim__> max = vlwNormalApproxMean0Var1[0].displacement;
	wpUtilities::Point<dim__> min = vlwNormalApproxMean0Var1[0].displacement;

	for (int i = 0; i < vlwNormalApproxMean0Var1.size(); ++i)
	{
		auto w1 = vlwNormalApproxMean0Var1[i].displacement[0];
		auto w2 = vlwNormalApproxMean0Var1[i].displacement[1];

		auto x = w1;
		auto y = rho * w1 + std::sqrt(1. - rho * rho) * w2;
		max[0] = std::max(x, max[0]);
		max[1] = std::max(y, max[1]);
		min[0] = std::min(x, min[0]);
		min[1] = std::min(y, min[1]);
	}
	spread[0] = max[0] - min[0];
	spread[1] = max[1] - min[1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// boundary data
//inline double BoundaryData(const double arg1, const double arg2)
//{
//	return std::exp(arg1 + arg2);
//}
//
//// asset spread option payoff log space
//inline double BoundaryData(const double arg1, const double arg2, const double dStrike)
//{
//	//return arg1 * arg2;
//	//return std::exp(arg1 + arg2);
//
//	double underlying = std::exp(arg1) - std::exp(arg2);
//	double ans = (dStrike - underlying);
//	return (ans < 0) ? 0 : ans;
//}
//
//// Put option spread
//inline double BoundaryData(const double arg1, const double arg2, const wpUtilities::Point<dim__> & dStrike)
//{
//	QL_REQUIRE(dStrike[1] > dStrike[0], "strike values inconsistent with assumption");
//	double underlying = std::exp(arg1) - std::exp(arg2);
//	//double long_call_payoff = ((underlying - dStrike[0]) < 0) ? 0 : underlying - dStrike[0];
//	//double short_call_payoff = ((underlying - dStrike[1]) < 0) ? 0 : underlying - dStrike[1];
//	//return long_call_payoff - short_call_payoff;
//	double short_put_payoff = ((dStrike[0] - underlying) < 0) ? 0 : dStrike[0] - underlying;
//	double long_put_payoff  = ((dStrike[1] - underlying) < 0) ? 0 : dStrike[1] - underlying;
//	return long_put_payoff - short_put_payoff;
//}


inline void makeSmaller(double& arg)
{
	arg -= 1e-15; // DBL_MIN; // std::numeric_limits<double>::min();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// log gbm 
// ln(x1) + (r - 0.5*sigma1*sigma1) t + sigma1 * B(t)
// ln(x1) + (r - 0.5*sigma2*sigma2) t + sigma2 * (rho * B(t) + sqrt(1 - rho*rho) * W(t))
// B(t), W(t) independent
//double CubaturePDESolver2D(double x1, double x2, wpUtilities::Point<dim__> strike, double t, double r, double sigma1, double sigma2, double rho,
//	CMainProgramParameters& ipIntParams, 
//	double e, 
//	unsigned int& DCount, unsigned int  stCubatureDegree)
//{
//	// fill quadrature weights
//	CWeightLocation2D::CVectorWeightLocation vlwNormalApproxMean0Var1;
//	QuadratureDoublePrecisionDeg5(vlwNormalApproxMean0Var1);
//	wpUtilities::Point<dim__> _spread = { 0, 0 };		// 1 step cubature (BM) spread
//	// compute cubature points fixed spread
//	ComputeCubatureSpread(vlwNormalApproxMean0Var1, rho, _spread);
//
//	CWeightLocation2D::CVectorWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation2D::CVectorWeightLocation vlwLocationsWeightsCarryForward;
//	CWeightLocation2D::CWeightLocation arg;
//	arg.displacement = { x1, x2 };
//	arg.probability = 1;
//	slwLocationsWeightsLeft.push_back(arg);
//	double dTimeLeft(t);
//	double dAccumulatedSolution(0);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		// drift
//		double drift1 = (r - 0.5 * sigma1 * sigma1) * dTimeLeft;
//		double drift11 = (r - 0.5 * sigma1 * sigma1) * ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft;
//		double drift12 = (r - 0.5 * sigma1 * sigma1) * (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dTimeLeft;
//
//		double drift2 = (r - 0.5 * sigma2 * sigma2) * dTimeLeft;
//		double drift21 = (r - 0.5 * sigma2 * sigma2) * ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft;
//		double drift22 = (r - 0.5 * sigma2 * sigma2) * (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dTimeLeft;
//
//		// std dev
//		double sd1 = sigma1 * sqrt(dTimeLeft);
//		double sd11 = sigma1 * sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd12 = sigma1 * sqrt((1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dTimeLeft);
//
//		double sd2 = sigma2 * sqrt(dTimeLeft);
//		double sd21 = sigma2 * sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd22 = sigma2 * sqrt((1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep) * dTimeLeft);
//		dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep);
//
//		for (CWeightLocation2D::CVectorWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			itExistingLocation != slwLocationsWeightsLeft.end();
//			++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CWeightLocation2D::CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//				itOuterNormal != vlwNormalApproxMean0Var1.end();
//				++itOuterNormal)
//			{
//				double _temp1 = itExistingLocation->displacement[0] + drift1 + sd1 * (*itOuterNormal).displacement[0];
//				double _temp2 = itExistingLocation->displacement[1] + drift2 + sd2 * (rho * (*itOuterNormal).displacement[0]
//					+ std::sqrt(1. - rho * rho) * (*itOuterNormal).displacement[1]);
//
//				ipIntParams.iCountBoundaryFunctionEvaluations++;
//				double temp = BoundaryData(_temp1, _temp2, strike);
//				ans += (*itOuterNormal).probability * temp;
//
//				for (
//					CWeightLocation2D::CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
//					itInnerNormal != vlwNormalApproxMean0Var1.end();
//					++itInnerNormal
//					)
//				{
//					ipIntParams.iCountBoundaryFunctionEvaluations++;
//					double __temp1 = itExistingLocation->displacement[0] + drift11 + sd11 * (*itOuterNormal).displacement[0]
//						+ drift12 + sd12 * (*itInnerNormal).displacement[0];
//					double __temp2 = itExistingLocation->displacement[1] + drift21
//						+ sd21 * (rho * (*itOuterNormal).displacement[0] + std::sqrt(1. - rho * rho) * (*itOuterNormal).displacement[1])
//						+ drift22
//						+ sd22 * (rho * (*itInnerNormal).displacement[0] + std::sqrt(1. - rho * rho) * (*itInnerNormal).displacement[1]);
//
//					double temp_ = BoundaryData(__temp1, __temp2, strike);
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * temp_;
//				}
//			}
//			double dError(0);
//			if (((dError = abs(ans - ans1)) <= e) || (DCount >= ipIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CWeightLocation2D::CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
//					itOuterNormal != vlwNormalApproxMean0Var1.end();
//					++itOuterNormal)
//				{
//					CWeightLocation2D::CWeightLocation lwNextWeightLocation;
//
//					lwNextWeightLocation.displacement[0] = itExistingLocation->displacement[0] + drift11 + sd11 * (*itOuterNormal).displacement[0];
//					lwNextWeightLocation.displacement[1] = itExistingLocation->displacement[1] + drift21
//						+ sd21 * (rho * (*itOuterNormal).displacement[0] + std::sqrt(1. - rho * rho) * (*itOuterNormal).displacement[1]);
//
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					//slwCarryForward.insert(lwNextWeightLocation);
//					vlwLocationsWeightsCarryForward.push_back(lwNextWeightLocation);
//				}
//			}
//		}
//
//		// rescale all carry forward values to [1,2)^2
//		// get the max and min, which are needed for the rescaling, also compute the maximum diameter of the points
//		std::vector<double> MAX = { 0, 0 };
//		if (vlwLocationsWeightsCarryForward.size() > 0)
//			MAX = { vlwLocationsWeightsCarryForward.begin()->displacement[0], vlwLocationsWeightsCarryForward.begin()->displacement[1]
//		};
//		std::vector<double> MIN(MAX);
//
//		for (size_t i = 0; i < dim__; ++i)
//		{
//			for (CWeightLocation2D::CVectorWeightLocation::const_iterator itCarryForwardLocation(vlwLocationsWeightsCarryForward.begin());
//				itCarryForwardLocation != vlwLocationsWeightsCarryForward.end();
//				++itCarryForwardLocation)
//			{
//				MAX[i] = std::max(itCarryForwardLocation->displacement[i], MAX[i]);
//				MIN[i] = std::min(itCarryForwardLocation->displacement[i], MIN[i]);
//			}
//		}
//
//		// make max slightly bigger than 2 -- for fixing exponent in bit representation
//		std::vector<double> MAX_enlarged(MAX);
//		makeLarger(MAX_enlarged[0]);
//		makeLarger(MAX_enlarged[1]);
//		wpUtilities::Point<dim__> _scaling = { sd11, sd21 };
//		double dNonzeroComponents = 0;
//
//		// dAccuracySpread - recombination patch size, set to rescale as sqrt(time remaining)
//		// compute dAccuracySpread consistent with rescaling.
//		// double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
//		//		rbegin()->displacement);     // original 1D dAccuracySpread
//		double dAccuracySpread(0.);
//		for (size_t i = 0; i < dim__; ++i)
//		{
//			dAccuracySpread += ((MAX[i] - MIN[i]) == 0.)
//				? 0. : std::pow(_scaling[i] * _spread[i] / (MAX_enlarged[i] - MIN[i]), 2);
//			dNonzeroComponents += ((MAX[i] - MIN[i]) == 0.) ? 0. : 1;
//		}
//		dAccuracySpread = std::sqrt(dAccuracySpread);
//
//		// get the spread in each component, since we scale all points to [1, 2)^2, spread must be < 1
//		double spread(0);
//		try
//		{
//			spread = (dNonzeroComponents > 0) ? dAccuracySpread * ipIntParams.dAdjustClusterDiameter / std::sqrt(dNonzeroComponents) : 0;
//			if (spread > 1)
//				throw(2);
//		}
//		catch (int i)
//		{
//			std::cout << "spread >= 1";
//			throw i;
//		}
//		// Obtain the bit mask value that gives the spread:  argmax_n {2^(-n) =< spread}
//		size_t uiBitMask = (spread > 0.) ? 52 + std::floor(std::log2(spread)) : 0;
//		//uiBitMask = 52;
//
//#if _VERBOSE_TEST_		
//		std::cout << "\n itOriginal" << std::endl;
//		for (CWeightLocation2D::CVectorWeightLocation::const_iterator itOriginal(vlwLocationsWeightsCarryForward.begin());
//			itOriginal != vlwLocationsWeightsCarryForward.end();
//			++itOriginal)
//			std::cout << "[" << itOriginal->displacement[0] << ", " << itOriginal->displacement[1] << "], ";
//		std::cout << std::endl;
//#endif
//		CWeightLocation2D::CCloudWeightLocation slwCarryForward;
//
//		// we make the choice to refactor the points to the dyadic set [1, 2)^d and add the points to slwCarryForward
//		// this is so that we can fix the exponent bits in order to use morton ordering to do partitioning in the ReducePointSet member function
//		// for (ptrdiff_t j = 0; j < ptrdiff_t(vlwLocationsWeightsCarryForward.size()); ++j)
//		for (CWeightLocation2D::CVectorWeightLocation::const_iterator itCarryForwardLocation(vlwLocationsWeightsCarryForward.begin());
//			itCarryForwardLocation != vlwLocationsWeightsCarryForward.end();
//			++itCarryForwardLocation)
//		{
//			CWeightLocation2D::CWeightLocation lwNextWeightLocation;
//			lwNextWeightLocation.probability = itCarryForwardLocation->probability;
//			for (size_t i = 0; i < dim__; ++i)
//				lwNextWeightLocation.displacement[i] = ((MAX[i] - MIN[i]) == 0.)
//				? 0. : (itCarryForwardLocation->displacement[i] + MAX_enlarged[i] - 2. * MIN[i]) / (MAX_enlarged[i] - MIN[i]);
//			// slwCarryForward carries the remapped points
//			slwCarryForward.insert(lwNextWeightLocation);
//		}
//
//#if _VERBOSE_TEST_
//		std::cout << "\n itSortedMorton" << std::endl;
//		for (auto itSortedMorton(slwCarryForward.begin()); itSortedMorton != slwCarryForward.end(); ++itSortedMorton)
//			std::cout << "[" << itSortedMorton->displacement[0] << ", " << itSortedMorton->displacement[1] << "], ";
//		std::cout << std::endl;
//#endif
//
//		if (ipIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//#if _VERBOSE_
//		std::cout << "TreeDepth: " << DCount;
//		std::cout << ", uiBitMask " << uiBitMask << ", patch diam "  << dAccuracySpread;  //<< spread << ", "
//#endif
//
//		//unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int)vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
//		//unsigned int stCubatureDegree = 5; // not sure if this works well? 
//		slwCarryForward.ReducePointSet(uiBitMask, stCubatureDegree); // , ipIntParams);
//
//		slwLocationsWeightsLeft.clear();
//		for (CWeightLocation2D::CCloudWeightLocation::const_iterator itCarryForwardLocation(slwCarryForward.begin());
//			itCarryForwardLocation != slwCarryForward.end();
//			++itCarryForwardLocation)
//		{
//			// rescale back to original cube
//			CWeightLocation2D::CWeightLocation lwNextWeightLocation;
//			lwNextWeightLocation.probability = itCarryForwardLocation->probability;
//			for (size_t i = 0; i < dim__; ++i)
//				lwNextWeightLocation.displacement[i] = ((MAX[i] - MIN[i]) == 0.) ? MIN[i]
//				: (MAX_enlarged[i] - MIN[i]) * itCarryForwardLocation->displacement[i] + 2. * MIN[i] - MAX_enlarged[i];
//			slwLocationsWeightsLeft.push_back(lwNextWeightLocation);
//		}
//		slwCarryForward.clear();
//		vlwLocationsWeightsCarryForward.clear();
//#if _VERBOSE_
//		std::cout << ", points remaining " << slwLocationsWeightsLeft.size() << std::endl;
//#endif
//		DCount++;
//#if _VERBOSE_TEST_
//		std::cout << "\n itScaledBack" << std::endl;
//		for (CWeightLocation2D::CVectorWeightLocation::const_iterator itScaledBack(slwLocationsWeightsLeft.begin());
//			itScaledBack != slwLocationsWeightsLeft.end();
//			++itScaledBack)
//			std::cout << "[" << itScaledBack->displacement[0] << ", " << itScaledBack->displacement[1] << "], ";
//		std::cout << std::endl;
//#endif
//	}
//	return dAccumulatedSolution*std::exp(-r*t);
//	//return dAccumulatedSolution;
//}
//
//
