#pragma once
#ifndef TemplateCubaturePDESolver_h__
#define TemplateCubaturePDESolver_h__

#include "MainProgramParameters.h"
#include "CloudWeightLocation2D.h"
#include "CloudWeightLocation.h"
#include "wpUtilities.h"	// points
#include "morton.h"

#include <iostream>

#define _VERBOSE_ 1

namespace MyPDESolver
{
	struct CubaturePDESolverParams
	{
		unsigned int uiRecombinationDegree;
		double dProportionOfRemainingTimeConsumedThisStep;
		double dAdjustClusterDiameter;
		unsigned int uiMaxEvaluationTreeDepth;
		unsigned int uiCountBoundaryFunctionEvaluations;
		double dTolerance;
	};

	class BrownianMotionWithDrift2D
	{
		typedef wpUtilities::Point_POD<dim__> Point_POD;
		typedef morton::block<Point_POD>::BLOCK BLOCK;
		typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
		typedef morton::block<Point_POD>::compare COMPARE;
		typedef wpUtilities::Point<dim__> Point;
		double r;
		double sigma1;
		double sigma2;
		double rho;
	public:
		BrownianMotionWithDrift2D(double _r, double _sigma1, double _sigma2, double _rho)
			: r(_r), sigma1(_sigma1), sigma2(_sigma2), rho(_rho)
		{}

		Point operator()(const Point & x0, const Point & control_path, const double t)
		{
			double drift1 = (r - 0.5*sigma1*sigma1)*t;
			double drift2 = (r - 0.5*sigma2*sigma2)*t;
			double sd1 = sigma1 * std::sqrt(t);
			double sd2 = sigma2 * std::sqrt(t);

			double _temp1 = x0[0] + drift1 + sd1 * control_path[0];
			double _temp2 = x0[1] + drift2 + sd2 * (rho * control_path[0]
				+ std::sqrt(1. - rho * rho) * control_path[1]);
			Point ans = { _temp1, _temp2 };
			return ans;
		}

		Point computeSpreadInEachComponent(const double t, const CWeightLocation2D::CVectorWeightLocation & oneStepCubatureFormula)
		{
			std::vector<Point> _solutions;

			for (auto it = oneStepCubatureFormula.begin(); it != oneStepCubatureFormula.end(); ++it)
				_solutions.push_back(this->operator()(Point({ 0, 0 }), (*it).displacement, t));

			Point max = *(_solutions.begin());
			Point min = *(_solutions.begin());
			for (auto it = _solutions.begin(); it != _solutions.end(); ++it)
			{
				auto x = (*it)[0];
				auto y = (*it)[1];

				max[0] = std::max(x, max[0]);
				max[1] = std::max(y, max[1]);
				min[0] = std::min(x, min[0]);
				min[1] = std::min(y, min[1]);
			}

			return Point({ max[0] - min[0], max[1] - min[1] });
		}
	};

	template <typename ODESolver, typename BoundaryFunction, typename BoundaryFunctionParams>
	class CubaturePDESolver
	{
		typedef wpUtilities::Point_POD<dim__> Point_POD;
		typedef morton::block<Point_POD>::BLOCK BLOCK;
		typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
		typedef morton::block<Point_POD>::compare COMPARE;
		typedef wpUtilities::Point<dim__> Point;
	public:
		CubaturePDESolver(CubaturePDESolverParams & solverparams, CWeightLocation2D::CVectorWeightLocation oneStepCubatureFormula, ODESolver & odeSolver, BoundaryFunction& boundaryFunc, BoundaryFunctionParams& boundaryFuncParams)
			: mBoundaryFunction(boundaryFunc)
			, mBoundaryFunctionParams(boundaryFuncParams)
			, mOneStepCubatureFormula(oneStepCubatureFormula)
			, mODESolver(odeSolver)
			, mSolverParams(solverparams)
		{}

		~CubaturePDESolver() {}

		double operator() (const double, const Point &);

	private:
		BoundaryFunction& mBoundaryFunction;
		const CWeightLocation2D::CVectorWeightLocation	mOneStepCubatureFormula;
		BoundaryFunctionParams& mBoundaryFunctionParams;
		ODESolver& mODESolver;
		CubaturePDESolverParams& mSolverParams;
	};
}


namespace MyPDESolver
{
	template <typename ODESolver, typename BoundaryFunction, typename BoundaryFunctionParams>
	double CubaturePDESolver<ODESolver, BoundaryFunction, BoundaryFunctionParams>::operator()(const double t, const Point & x0)
	{
		unsigned int DCount(0);
		CWeightLocation2D::CVectorWeightLocation slwLocationsWeightsLeft;
		CWeightLocation2D::CVectorWeightLocation vlwLocationsWeightsCarryForward;
		CWeightLocation2D::CWeightLocation arg;
		arg.displacement = x0;
		arg.probability = 1;
		slwLocationsWeightsLeft.push_back(arg);
		double dTimeLeft(t);
		double dAccumulatedSolution(0);

		while (slwLocationsWeightsLeft.size() != 0)
		{
			for (CWeightLocation2D::CVectorWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
				itExistingLocation != slwLocationsWeightsLeft.end();
				++itExistingLocation)
			{
				double ans(0), ans1(0);
				// compute both approximations
				for (CWeightLocation2D::CVectorWeightLocation::const_iterator itOuterNormal(mOneStepCubatureFormula.begin());
					itOuterNormal != mOneStepCubatureFormula.end();
					++itOuterNormal)
				{
					mSolverParams.uiCountBoundaryFunctionEvaluations++;
					// move forward in one step
					Point _outerpt(this->mODESolver(itExistingLocation->displacement, itOuterNormal->displacement, dTimeLeft));
					double _temp = this->mBoundaryFunction(_outerpt, this->mBoundaryFunctionParams);
					ans += (*itOuterNormal).probability * _temp;

					for (CWeightLocation2D::CVectorWeightLocation::const_iterator itInnerNormal(mOneStepCubatureFormula.begin());
						itInnerNormal != mOneStepCubatureFormula.end();
						++itInnerNormal)
					{
						mSolverParams.uiCountBoundaryFunctionEvaluations++;
						// move forward in two steps
						Point __inner_first_step(this->mODESolver(itExistingLocation->displacement, itOuterNormal->displacement, 
											dTimeLeft* mSolverParams.dProportionOfRemainingTimeConsumedThisStep));
						Point __innerpt(this->mODESolver(__inner_first_step, itInnerNormal->displacement,
											dTimeLeft *(1. - mSolverParams.dProportionOfRemainingTimeConsumedThisStep)));
						double __temp = this->mBoundaryFunction(__innerpt, this->mBoundaryFunctionParams);
						ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * __temp;
					}
				}
				double dError(0);
				if (((dError = abs(ans - ans1)) <= mSolverParams.dTolerance) || (DCount >= mSolverParams.uiMaxEvaluationTreeDepth))
					dAccumulatedSolution += itExistingLocation->probability * ans1;
				else
				{
					for (CWeightLocation2D::CVectorWeightLocation::const_iterator itOuterNormal(mOneStepCubatureFormula.begin());
						itOuterNormal != mOneStepCubatureFormula.end();
						++itOuterNormal)
					{
						CWeightLocation2D::CWeightLocation lwNextWeightLocation;
						lwNextWeightLocation.displacement = this->mODESolver(itExistingLocation->displacement, itOuterNormal->displacement,
																dTimeLeft * mSolverParams.dProportionOfRemainingTimeConsumedThisStep);
						lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
						vlwLocationsWeightsCarryForward.push_back(lwNextWeightLocation);
					}
				}
			}

			// rescale all carry forward values to [1,2)^2
			// get the max and min, which are needed for the rescaling, also compute the maximum diameter of the points
			Point MAX = { 0,0 };
			if (vlwLocationsWeightsCarryForward.size() > 0)
				MAX = vlwLocationsWeightsCarryForward.begin()->displacement;
			Point MIN(MAX);

			for (size_t i = 0; i < dim__; ++i)
			{
				for (CWeightLocation2D::CVectorWeightLocation::const_iterator itCarryForwardLocation(vlwLocationsWeightsCarryForward.begin());
					itCarryForwardLocation != vlwLocationsWeightsCarryForward.end();
					++itCarryForwardLocation)
				{
					MAX[i] = std::max(itCarryForwardLocation->displacement[i], MAX[i]);
					MIN[i] = std::min(itCarryForwardLocation->displacement[i], MIN[i]);
				}
			}

			// make max slightly bigger than 2 -- for fixing exponent in bit representation
			// TODO: does this work in negative space?
			Point MAX_enlarged(MAX);
			std::for_each(MAX_enlarged.begin(), MAX_enlarged.end(), [](double & element) { element += 1e-15; });
			double dNonzeroComponents = 0;
			Point ptSpreads = this->mODESolver.computeSpreadInEachComponent(mSolverParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft, mOneStepCubatureFormula);
			// dAccuracySpread - recombination patch size, set to rescale as sqrt(time remaining)
			// compute dAccuracySpread consistent with rescaling.
			// double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.rbegin()->displacement);     // original 1D dAccuracySpread
			double dAccuracySpread(0.);
			for (size_t i = 0; i < dim__; ++i)
			{
				dAccuracySpread += ((MAX[i] - MIN[i]) == 0.)
					? 0. : std::pow(ptSpreads[i] / (MAX_enlarged[i] - MIN[i]), 2);
				dNonzeroComponents += ((MAX[i] - MIN[i]) == 0.) ? 0. : 1;
			}
			dAccuracySpread = std::sqrt(dAccuracySpread);

			// get the spread in each component, since we scale all points to [1, 2)^2, spread must be < 1
			double spread(0);
			try
			{
				spread = (dNonzeroComponents > 0) ? dAccuracySpread * mSolverParams.dAdjustClusterDiameter / std::sqrt(dNonzeroComponents) : 0;
				if (spread > 1)
					throw(2);
			}
			catch (int i)
			{
				std::cout << "spread >= 1";
				throw i;
			}
			// Obtain the bit mask value that gives the spread:  argmax_n {2^(-n) =< spread}
			size_t uiBitMask = (spread > 0.) ? 52 + std::floor(std::log2(spread)) : 0;

			CWeightLocation2D::CCloudWeightLocation slwCarryForward;

			// we make the choice to refactor the points to the dyadic set [1, 2)^d and add the points to slwCarryForward
			// this is so that we can fix the exponent bits in order to use morton ordering to do partitioning in the ReducePointSet member function
			for (CWeightLocation2D::CVectorWeightLocation::const_iterator itCarryForwardLocation(vlwLocationsWeightsCarryForward.begin());
				itCarryForwardLocation != vlwLocationsWeightsCarryForward.end();
				++itCarryForwardLocation)
			{
				CWeightLocation2D::CWeightLocation lwNextWeightLocation;
				lwNextWeightLocation.probability = itCarryForwardLocation->probability;
				for (size_t i = 0; i < dim__; ++i)
					lwNextWeightLocation.displacement[i] = ((MAX[i] - MIN[i]) == 0.)
					? 0. : (itCarryForwardLocation->displacement[i] + MAX_enlarged[i] - 2. * MIN[i]) / (MAX_enlarged[i] - MIN[i]);
				// slwCarryForward carries the remapped points
				slwCarryForward.insert(lwNextWeightLocation);
			}
#if _VERBOSE_
			std::cout << "TreeDepth: " << DCount;
			std::cout << ", uiBitMask " << uiBitMask << ", patch diam " << dAccuracySpread;  //<< spread << ", "
#endif

		//unsigned int iNoCubatureToDimension = (ipIntParams.INumerator() * (int)vlwNormalApproxMean0Var1.size()) / ipIntParams.IDenominator();
		//unsigned int stCubatureDegree = 5; // not sure if this works well? 
			slwCarryForward.ReducePointSet(uiBitMask, this->mSolverParams.uiRecombinationDegree); // , ipIntParams);

			slwLocationsWeightsLeft.clear();
			for (CWeightLocation2D::CCloudWeightLocation::const_iterator itCarryForwardLocation(slwCarryForward.begin());
				itCarryForwardLocation != slwCarryForward.end();
				++itCarryForwardLocation)
			{
				// rescale back to original cube
				CWeightLocation2D::CWeightLocation lwNextWeightLocation;
				lwNextWeightLocation.probability = itCarryForwardLocation->probability;
				for (size_t i = 0; i < dim__; ++i)
					lwNextWeightLocation.displacement[i] = ((MAX[i] - MIN[i]) == 0.) ? MIN[i]
					: (MAX_enlarged[i] - MIN[i]) * itCarryForwardLocation->displacement[i] + 2. * MIN[i] - MAX_enlarged[i];
				slwLocationsWeightsLeft.push_back(lwNextWeightLocation);
			}
			slwCarryForward.clear();
			vlwLocationsWeightsCarryForward.clear();
#if _VERBOSE_
			std::cout << ", points remaining " << slwLocationsWeightsLeft.size() << std::endl;
#endif
			dTimeLeft *= (1. - mSolverParams.dProportionOfRemainingTimeConsumedThisStep); // next
			DCount++;
		}
		return dAccumulatedSolution;
	}
}

//
///// boundary function is assumed to be a unary function object
//template <typename BoundaryFunc>
//class PDESolver
//{
//public:
//	PDESolver(BoundaryFunc& terminalFunction, CVectorWeightLocation vlwNormalApproxMean0Var1, CMainProgramParameters ipIntParams, double e)
//		: mTerminalFunction(terminalFunction)
//		, mIpIntParams(ipIntParams)
//		, mVlwNormalApproxMean0Var1(vlwNormalApproxMean0Var1)
//		, mError(e)
//	{}
//
//	~PDESolver() {}
//
//
//	double operator() (double spot)
//	{
//		return CubaturePDESolver(spot);
//	}
//
//
//	unsigned int getCountBoundaryFunctionEvaluations() const
//	{
//		return mIpIntParams.iCountBoundaryFunctionEvaluations;
//	}
//
//
//private:
//
//	BoundaryFunc& mTerminalFunction;
//	const CVectorWeightLocation		mVlwNormalApproxMean0Var1;
//	CMainProgramParameters			mIpIntParams;
//	double							mError;
//
//
//	double CubaturePDESolver(double spot);
//
//};
//
//
//template <typename BoundaryFunc>
//class PDESolverDriftBM
//{
//public:
//	PDESolverDriftBM(BoundaryFunc& terminalFunction, CVectorWeightLocation vlwNormalApproxMean0Var1, CMainProgramParameters ipIntParams, double e)
//		: mTerminalFunction(terminalFunction)
//		, mIpIntParams(ipIntParams)
//		, mVlwNormalApproxMean0Var1(vlwNormalApproxMean0Var1)
//		, mError(e)
//	{}
//
//	~PDESolverDriftBM() {}
//
//	//template <typename OptData>
//	double operator() (double log_spot)
//	{
//		return CubaturePDESolver(log_spot);
//	}
//	/*
//		double operator() ( OptData & opt)
//		{
//			return CubaturePDESolver(std::log(opt.spot), opt);
//		}*/
//
//
//	unsigned int getCountBoundaryFunctionEvaluations() const
//	{
//		return mIpIntParams.iCountBoundaryFunctionEvaluations;
//	}
//
//private:
//
//	BoundaryFunc& mTerminalFunction;
//	const CVectorWeightLocation		mVlwNormalApproxMean0Var1;
//	CMainProgramParameters			mIpIntParams;
//	double							mError;
//
//	// unary boundary condition
//	double CubaturePDESolver(double log_spot);
//
//};
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////	implementation
//
//
///// brownian motion with unary terminal condition
//template <typename BoundaryFunc>
//double PDESolver<BoundaryFunc>::CubaturePDESolver(double spot)
//{
//	unsigned int DCount(0);
//
//	CSetWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation arg;
//	arg.displacement = spot;
//	arg.probability = 1;
//	slwLocationsWeightsLeft.insert(arg);
//	double dTimeLeft(mIpIntParams.dTime);
//	double dAccumulatedSolution(0);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		double sd = mIpIntParams.dVol * sqrt(dTimeLeft);
//		double sd1 = mIpIntParams.dVol * sqrt(mIpIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd2 = mIpIntParams.dVol * sqrt(dTimeLeft *= (1. - mIpIntParams.dProportionOfRemainingTimeConsumedThisStep));
//
//		CCloudWeightLocation slwCarryForward;
//
//		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			itExistingLocation != slwLocationsWeightsLeft.end();
//			++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CVectorWeightLocation::const_iterator itOuterNormal(mVlwNormalApproxMean0Var1.begin());
//				itOuterNormal != mVlwNormalApproxMean0Var1.end();
//				++itOuterNormal)
//			{
//				double outerPathValue = itExistingLocation->displacement + sd * (*itOuterNormal).displacement;
//				ans += (*itOuterNormal).probability * mTerminalFunction(outerPathValue);
//
//				for (CVectorWeightLocation::const_iterator itInnerNormal(mVlwNormalApproxMean0Var1.begin());
//					itInnerNormal != mVlwNormalApproxMean0Var1.end();
//					++itInnerNormal)
//				{
//					mIpIntParams.iCountBoundaryFunctionEvaluations++;
//					double innerPathValue = itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement
//						+ sd2 * (*itInnerNormal).displacement;
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * mTerminalFunction(innerPathValue);
//				}
//			}
//
//			double dError;
//
//			if (((dError = abs(ans - ans1)) <= mError) || (DCount >= mIpIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CVectorWeightLocation::const_iterator itOuterNormal(mVlwNormalApproxMean0Var1.begin());
//					itOuterNormal != mVlwNormalApproxMean0Var1.end();
//					++itOuterNormal)
//				{
//					CWeightLocation lwNextWeightLocation;
//					lwNextWeightLocation.displacement = itExistingLocation->displacement + sd1 * (*itOuterNormal).displacement;
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					slwCarryForward.insert(lwNextWeightLocation);
//				}
//			}
//		}
//
//		if (mIpIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//		double dAccuracySpread = sd1 * (mVlwNormalApproxMean0Var1.begin()->displacement - mVlwNormalApproxMean0Var1.
//			rbegin()->displacement);
//
//		slwCarryForward.ReducePointSet
//		(dAccuracySpread * mIpIntParams.dAdjustClusterDiameter
//			, (unsigned int)(mIpIntParams.INumerator() * mVlwNormalApproxMean0Var1
//				.size()) / mIpIntParams.IDenominator()
//			, mIpIntParams);
//
//		slwLocationsWeightsLeft.swap(slwCarryForward);
//		slwCarryForward.clear();
//		//	std::cout << slwLocationsWeightsLeft.size() << std::endl;
//		DCount++;
//	}	// end while
//
//	//std::cout << "tree depth reached is " << DCount << std::endl;
//
//	return dAccumulatedSolution;
//
//} /// end bm solver
//
//
///// brownian motion with drift with unary terminal condition
//template <typename BoundaryFunc>
//double PDESolverDriftBM<BoundaryFunc>::CubaturePDESolver(double log_spot)
//{
//	unsigned int DCount(0);
//
//	CSetWeightLocation slwLocationsWeightsLeft;
//	CWeightLocation arg;
//	arg.displacement = log_spot;
//	arg.probability = 1;
//	slwLocationsWeightsLeft.insert(arg);
//	double dTimeLeft(mIpIntParams.dTime);
//	double dAccumulatedSolution(0);
//	double df = std::exp(-mIpIntParams.dRate * dTimeLeft);
//
//	while (slwLocationsWeightsLeft.size() != 0)
//	{
//		// drift
//		double drift = (mIpIntParams.dRate - .5 * mIpIntParams.dVol * mIpIntParams.dVol) * dTimeLeft;
//		double drift1 = drift * mIpIntParams.dProportionOfRemainingTimeConsumedThisStep;
//		double drift2 = drift * (1. - mIpIntParams.dProportionOfRemainingTimeConsumedThisStep);
//
//		// std dev
//		double sd = mIpIntParams.dVol * sqrt(dTimeLeft);
//		double sd1 = mIpIntParams.dVol * sqrt(mIpIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
//		double sd2 = mIpIntParams.dVol * sqrt(dTimeLeft *= (1. - mIpIntParams.dProportionOfRemainingTimeConsumedThisStep));
//
//		CCloudWeightLocation slwCarryForward;
//
//		for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
//			itExistingLocation != slwLocationsWeightsLeft.end();
//			++itExistingLocation)
//		{
//			double ans(0), ans1(0);
//			// compute both approximations
//			for (CVectorWeightLocation::const_iterator itOuterNormal(mVlwNormalApproxMean0Var1.begin());
//				itOuterNormal != mVlwNormalApproxMean0Var1.end();
//				++itOuterNormal)
//			{
//				ans += (*itOuterNormal).probability * mTerminalFunction(itExistingLocation->displacement + drift + sd * (*itOuterNormal).displacement);
//
//				for (CVectorWeightLocation::const_iterator itInnerNormal(mVlwNormalApproxMean0Var1.begin());
//					itInnerNormal != mVlwNormalApproxMean0Var1.end();
//					++itInnerNormal)
//				{
//					mIpIntParams.iCountBoundaryFunctionEvaluations++;
//					ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * mTerminalFunction(itExistingLocation->displacement + drift1 + sd1 * (*itOuterNormal).displacement
//						+ drift2 + sd2 * (*itInnerNormal).displacement); // path concatenation
//				}
//			}
//			if ((abs(ans - ans1) <= mError) || (DCount >= mIpIntParams.iMaxEvaluationTreeDepth))
//			{
//				dAccumulatedSolution += itExistingLocation->probability * ans1;
//			}
//			else
//			{
//				for (CVectorWeightLocation::const_iterator itOuterNormal(mVlwNormalApproxMean0Var1.begin());
//					itOuterNormal != mVlwNormalApproxMean0Var1.end();
//					++itOuterNormal)
//				{
//					CWeightLocation lwNextWeightLocation;
//					lwNextWeightLocation.displacement = itExistingLocation->displacement + drift1 + sd1 * (*itOuterNormal).displacement;
//					lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
//					slwCarryForward.insert(lwNextWeightLocation);
//				}
//			}
//		}
//
//		if (mIpIntParams.bVerbose)
//			std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;
//
//		double dAccuracySpread = sd1 * (mVlwNormalApproxMean0Var1.begin()->displacement - mVlwNormalApproxMean0Var1.
//			rbegin()->displacement);
//
//		slwCarryForward.ReducePointSet
//		(dAccuracySpread * mIpIntParams.dAdjustClusterDiameter
//			, (mIpIntParams.INumerator() * (unsigned int)mVlwNormalApproxMean0Var1.size()) / mIpIntParams.IDenominator()
//			, mIpIntParams);
//
//		slwLocationsWeightsLeft.swap(slwCarryForward);
//		slwCarryForward.clear();
//		DCount++;
//	}
//
//	return dAccumulatedSolution * df;
//}	/// end bm drift solver

#endif