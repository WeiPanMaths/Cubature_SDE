/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#ifndef CubaturePDESolver_h__
#define CubaturePDESolver_h__
//#include "SharedFunctionsAndTypes.h" // CVectorWeightLocation
#include "MainProgramParameters.h"
#include "CloudWeightLocation.h"
#include "Utils.h"

// solves u_t + .5 u_xx = 0
double CubaturePDESolver(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
						 CMainProgramParameters& ipIntParams, double e, unsigned int & DCount);

template <typename F, typename Datatype>
double CubaturePDESolverWithAdaptiveApprox(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
	CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int& DCount, F& funcwrapper, Datatype & data)
{
	e = ipIntParams.tolerance = (e == 0. && ipIntParams.tolerance == 0.)
		? ipIntParams.dErrorAcceptanceFactor * abs(CompareCubatureAndTrueSolution(ipIntParams.
			dProportionOfRemainingTimeConsumedThisStep, t, x, vlwNormalApproxMean0Var1))
		: (e != 0)
		? e
		: ipIntParams.tolerance;

	//LogArgLessExp myfunc;
	//FunctionData data;
	//ApproximationMorton<LogArgLessExp, FunctionData > funcwrapper(myfunc, 6, std::pow(10., -12));	// wrapper function algorithm

	CSetWeightLocation slwLocationsWeightsLeft;
	CWeightLocation arg;
	arg.displacement = x;
	arg.probability = 1;
	slwLocationsWeightsLeft.insert(arg);
	double dTimeLeft(t);
	double dAccumulatedSolution(0);

	while (slwLocationsWeightsLeft.size() != 0)
	{
		double sd = sqrt(dTimeLeft);
		double sd1 = sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
		double sd2 = sqrt(dTimeLeft *= (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));

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
				point pt_ = { itExistingLocation->displacement + sd * (*itOuterNormal).displacement };
				double temp = funcwrapper(pt_, data);
				ans += (*itOuterNormal).probability * temp;

				for (
					CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
					itInnerNormal != vlwNormalApproxMean0Var1.end();
					++itInnerNormal
					)
				{
					ipIntParams.iCountBoundaryFunctionEvaluations++;
					point pt_outer = { itExistingLocation->displacement
							+ sd1 * (*itOuterNormal).displacement
							+ sd2 * (*itInnerNormal).displacement };
					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
						funcwrapper(pt_outer, data); // path concatenation
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
					lwNextWeightLocation.displacement = itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement;
					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
					slwCarryForward.insert(lwNextWeightLocation);
				}
			}
		}
		double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
			rbegin()->displacement);

		if (ipIntParams.bVerbose)
			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;

		unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int)vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
		double spread = dAccuracySpread * ipIntParams.dAdjustClusterDiameter;
		slwCarryForward.ReducePointSet
		(spread
			, iNoCubatureToDimension
			, ipIntParams);

		slwLocationsWeightsLeft.swap(slwCarryForward);
		slwCarryForward.clear();
		DCount++;
	}

	return dAccumulatedSolution;
}


// solves  u_t + r u_x + .5 sigma^2  u_xx = 0
double CubaturePDESolver(double x, double t, double r, double sigma, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
						 CMainProgramParameters& ipIntParams, double e, unsigned int & DCount);

// solves black scholes pde call
double CubaturePDESolver(double x, double strike, double t, double r, double sigma,
						 const CVectorWeightLocation & vlwNormalApproxMean0Var1,
						 CMainProgramParameters& ipIntParams, double e);

// solves black scholes pde with reduced parameter set
double CubaturePDESolver(CMainProgramParameters& ipIntParams, const CVectorWeightLocation& vlwNormalApproxMean0Var1, double e);

// solve black scholes pde at expiry
//double CubaturePDESolver(double x, double strike, double barrier, double t, double r, double sigma, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters& ipIntParams, double e);
double CubaturePDESolver(double x, double strike, double barrier, double t, double r, double sigma, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters& ipIntParams, double e, unsigned int & DCount);
double CubaturePDESolverThreeStepOneStep(double x, double t, double r, double sigma, const CVectorWeightLocation& vlwNormalApproxMean0Var1, CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount );
double CubaturePDESolverThreeStepOneStep(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1, CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int & DCount );

#if dim_ == 1
template <typename F, typename Datatype>
double CubaturePDESolverThreeStepOneStepAA(double x, double t, const CVectorWeightLocation& vlwNormalApproxMean0Var1,
	CMainProgramParameters& ipIntParams, double e/*= 0.*/, unsigned int& DCount, F& funcwrapper, Datatype& data)
{
	//CDebugTools& dbgLog = ipIntParams.dbgLogFile;

#ifdef TJLVERBOSE
	//	dbgLog.m_out << "TJLVERBOSE defined in libs, Producing Debug Data" << std::endl;
#endif

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
		double stepOneSD = sqrt(dTimeLeft * ipIntParams.dProportionOfRemainingTimeConsumedThisStep);
		double stepTwoSD = sqrt(dAfterStepOneTimeLeft);
		double stepTwoSD1 = sqrt(ipIntParams.dProportionOfRemainingTimeConsumedThisStep * dAfterStepOneTimeLeft);
		double stepTwoSD2 = sqrt(dAfterStepOneTimeLeft * (1. - ipIntParams.dProportionOfRemainingTimeConsumedThisStep));
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
				lwNextStepOneWeightLocation.displacement = itExistingLocation->displacement + stepOneSD * (*itStepOneNormal).displacement;
				lwNextStepOneWeightLocation.probability = itExistingLocation->probability * (*itStepOneNormal).probability;
				slwStepOneLocationsWeights.insert(lwNextStepOneWeightLocation);
			}

			bool stepFwdFlag(0);
			double val(0);

			// for each intermediate one step...do one-step-two-step comparison
			for (CSetWeightLocation::iterator itExistingLocation = slwStepOneLocationsWeights.begin();
				itExistingLocation != slwStepOneLocationsWeights.end() && stepFwdFlag == 0;
				++itExistingLocation)
			{
				double ans(0), ans1(0);
				// compute both approximations
				for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
					itOuterNormal != vlwNormalApproxMean0Var1.end();
					++itOuterNormal)
				{
					point pt_ = { itExistingLocation->displacement + stepTwoSD * (*itOuterNormal).displacement };
					double temp = funcwrapper(pt_, data);
					ans += (*itOuterNormal).probability * temp;

					for (CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
						itInnerNormal != vlwNormalApproxMean0Var1.end();
						++itInnerNormal)
					{
						ipIntParams.iCountBoundaryFunctionEvaluations++;
						point pt_outer = { itExistingLocation->displacement + stepTwoSD1 * (*itOuterNormal).
								displacement + stepTwoSD2 * (*itInnerNormal).displacement };
						ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability *
							funcwrapper(pt_outer, data); // path concatenation
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
			}
			else
				dAccumulatedSolution += val;	// don't need to further expand so add to solution

			slwStepOneLocationsWeights.clear();

		} // end for each to be expanded

		if (ipIntParams.bVerbose)
			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;

		unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int)vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();

		//std::cout<< "number of particles is " << slwLocationsWeightsLeft.size() << std::endl;
#if 0
		// patch determination
		double patchEpsilon(1.0e-10);
		int n = 1 + iNoCubatureToDimension;
		double spreadConst = std::sqrt(2.) * std::pow(patchEpsilon * std::sqrt(boost::math::constants::pi<double>()) * (double)(n)*boost::math::factorial<double>((n - 1) / 2), 1. / (double)(n));
		double spread = spreadConst * std::pow(dTimeLeft, 0.5);

		//spread = std::pow(dTimeLeft,0.5); // test

		slwCarryForward.ReducePointSet
		(spread
			, iNoCubatureToDimension
			, ipIntParams);

#else
		double dAccuracySpread = stepOneSD * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.rbegin()->displacement);
		double spread = dAccuracySpread * ipIntParams.dAdjustClusterDiameter;
		slwCarryForward.ReducePointSet
		(spread
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
#endif

#endif // CubaturePDESolver_h__
