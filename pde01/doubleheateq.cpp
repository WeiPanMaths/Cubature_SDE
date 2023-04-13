// binary_recursion.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "CloudWeightLocation.h" // CVectorWeightLocation
#include "SharedFunctionsAndTypes.h" // operator >>
#include "CubaturePDESolver.h" // CubaturePDESolver
#include "CubaturePDESolver2D.h" // CubaturePDESolver2D
#include "BasketSpreadQuantLib.h"
#include "ExactPDESolution.h" // ExactPDESolution
#include "MainProgramParameters.h" //CMainProgramParameters
#include "DoublePrecisionGaussianQuadrature.h" //Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "trapfpe.h"
#include <omp.h>
#include "AdaptiveApproximation.h"
#include "TemplateCubaturePDESolver.h"

#define LENGTH(a) (sizeof(a)/sizeof(a[0]))
#define sanity 0

// spread put on asset spreads
struct BenchmarkBoundaryData
{
	wpUtilities::Point<dim__> _strike;
	BenchmarkBoundaryData(const wpUtilities::Point<dim__>& dStrike) : _strike(dStrike) {}

	double BoundaryData(const double arg1, const double arg2) const
	{

		QL_REQUIRE(_strike[1] > _strike[0], "strike values inconsistent with assumption");
		double underlying = std::exp(arg1) - std::exp(arg2);
		//double long_call_payoff = ((underlying - dStrike[0]) < 0) ? 0 : underlying - dStrike[0];
		//double short_call_payoff = ((underlying - dStrike[1]) < 0) ? 0 : underlying - dStrike[1];
		//return long_call_payoff - short_call_payoff;
		double short_put_payoff = ((_strike[0] - underlying) < 0) ? 0 : _strike[0] - underlying;
		double long_put_payoff = ((_strike[1] - underlying) < 0) ? 0 : _strike[1] - underlying;
		return long_put_payoff - short_put_payoff;
	}

	double operator() (const wpUtilities::Point<dim__>& pt, FunctionData& data) const
	{
		//std::cout << "here" << std::endl;
		return BoundaryData(pt[0], pt[1]);
	}
};

struct SpreadBoundaryData
{
	wpUtilities::Point<dim__> _strike;
	double _r;
	double _t;
	SpreadBoundaryData(const wpUtilities::Point<dim__>& dStrike, const double dr, const double t) : _strike(dStrike), _r(dr), _t(t) {}

	double BoundaryData(const double arg1, const double arg2) const
	{

		QL_REQUIRE(_strike[1] > _strike[0], "strike values inconsistent with assumption");
		double underlying = std::exp(arg1) - std::exp(arg2);
		//double long_call_payoff = ((underlying - dStrike[0]) < 0) ? 0 : underlying - dStrike[0];
		//double short_call_payoff = ((underlying - dStrike[1]) < 0) ? 0 : underlying - dStrike[1];
		//return long_call_payoff - short_call_payoff;
		double short_put_payoff = ((_strike[0] - underlying) < 0) ? 0 : _strike[0] - underlying;
		double long_put_payoff = ((_strike[1] - underlying) < 0) ? 0 : _strike[1] - underlying;
		return (long_put_payoff - short_put_payoff)*std::exp(-_r*_t);
	}

	double operator() (const wpUtilities::Point<dim__>& pt, FunctionData& data) const
	{
		//std::cout << "here" << std::endl;
		return BoundaryData(pt[0], pt[1]);
	}
};


void TestTemplateCubatureSolver()
{
	MyBasket::BasketOptionTwoData values[] = {
		{MyBasket::SpreadBasket, Option::Put, 1.0,  133.0, 130.0, 0.0, 0.0, 0.0,  0.1, 0.25, 0.20,  0.5, 3.7970, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Put, 4.0,  133.0, 130.0, 0.0, 0.0, 0.0,  0.1, 0.25, 0.20,  0.5, 3.7970, 1.0e-3}
	};

	for (Size i = 0; i < LENGTH(values); i += 2)
	{
		MyBasket::MyBasketSolver solver(values[i]);
		MyBasket::MyBasketSolver solver2(values[i + 1]);
		double ans_quantlib = -solver.CalculationKirk();
		ans_quantlib += solver2.CalculationKirk();

		std::cout << "\ncubature " << std::endl;

		double x1 = std::log(values[i].s1);
		double x2 = std::log(values[i].s2);
		double t = values[i].t;
		double r = values[i].r;
		double sigma1 = values[i].v1;
		double sigma2 = values[i].v2;
		double rho = values[i].rho;

		wpUtilities::Point<dim__> strike = { values[i].strike, values[i + 1].strike };

		FunctionData data;
		SpreadBoundaryData boundaryData(strike, r, t);
		MyPDESolver::BrownianMotionWithDrift2D bm(r, sigma1, sigma2, rho);
		MyPDESolver::CubaturePDESolverParams solverparams;
		CWeightLocation2D::CVectorWeightLocation onestepcubature;
		QuadratureDoublePrecisionDeg5(onestepcubature);

		solverparams.uiRecombinationDegree = 5;
		solverparams.dProportionOfRemainingTimeConsumedThisStep = 0.25;
		solverparams.dAdjustClusterDiameter = 0.75;
		solverparams.uiMaxEvaluationTreeDepth = 15;
		solverparams.uiCountBoundaryFunctionEvaluations = 0;
		solverparams.dTolerance = 1e-5;

		double dSecondsElapsed(0), ans(0);
		MyPDESolver::CubaturePDESolver<MyPDESolver::BrownianMotionWithDrift2D, SpreadBoundaryData, FunctionData> cubsolver(solverparams, onestepcubature, bm, boundaryData, data);
		Point x0 = { x1, x2 };
		{
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			ans = cubsolver(t, x0);
		}

		//double ans_exact = values[i].result; 
#if sanity
		double ans_exact = ExactSolution(x1, x2, r, sigma1, sigma2, rho, t);
#endif
		std::cout << "\nElapsed time in seconds " << dSecondsElapsed << "\n";
		std::cout << "deg5 cubature approx: " << ans << std::endl;
#if sanity
		std::cout << "exact solution      : " << ans_exact << std::endl;
		std::cout << "rel error           : " << std::abs(ans - ans_exact) / ans_exact << std::endl;
#else
		std::cout << "QuantLib approx     : " << ans_quantlib << std::endl;
		//std::cout << "QuantLib finite diff: " << fd_approx << std::endl;
		std::cout << "relerr (vs quantlib): " << std::abs(ans - ans_quantlib) / ans_quantlib << std::endl;
		//std::cout << "| cub - fd |/fd   : " << std::abs(ans - fd_approx) / fd_approx << std::endl;
		//std::cout << "| fd - kirk|/kirk : " << std::abs(ans_kirk - fd_approx) / ans_kirk << std::endl;
#endif
	}
}

void Test(CMainProgramParameters & mppProgramParameters)
{
	/*
		Data from:
		Excel spreadsheet www.maths.ox.ac.uk/~firth/computing/excel.shtml
		and
		"Option pricing formulas", E.G. Haug, McGraw-Hill 1998 pag 56-58
		European two asset max basket options
	*/
	MyBasket::BasketOptionTwoData values[] = {
		//      basketType,   optionType, strike,    s1,    s2,   q1,   q2,    r,    t,   v1,   v2,  rho, result, tol
		// data from http://www.maths.ox.ac.uk/~firth/computing/excel.shtml

		//      basketType,   optionType, strike,    s1,    s2,   q1,   q2,    r,    t,   v1,   v2,  rho,  result, tol
		// data from "Option pricing formulas" VB code + spreadsheet

		/* "Option pricing formulas", E.G. Haug, McGraw-Hill 1998 pag 59-60
			Kirk approx. for a european spread option on two futures*/
#if 0
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20, -0.5, 4.7530, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20,  0.0, 3.7970, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20,  0.5, 2.5537, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20, -0.5, 5.4275, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20,  0.0, 4.3712, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20,  0.5, 3.0086, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25, -0.5, 5.4061, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25,  0.0, 4.3451, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25,  0.5, 2.9723, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20, -0.5,10.7517, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20,  0.0, 8.7020, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20,  0.5, 6.0257, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20, -0.5,12.1941, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20,  0.0, 9.9340, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20,  0.5, 7.0067, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25, -0.5,12.1483, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25,  0.0, 9.8780, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25,  0.5, 6.9284, 1.0e-3}
#endif
		{MyBasket::SpreadBasket, Option::Put, 1.0,  133.0, 130.0, 0.0, 0.0, 0.0,  0.1, 0.25, 0.20,  0.5, 3.7970, 1.0e-3},
		{MyBasket::SpreadBasket, Option::Put, 4.0,  133.0, 130.0, 0.0, 0.0, 0.0,  0.1, 0.25, 0.20,  0.5, 3.7970, 1.0e-3}
	};

	for (Size i = 0; i < LENGTH(values); i += 2)
	{
		MyBasket::MyBasketSolver solver(values[i]);
		MyBasket::MyBasketSolver solver2(values[i+1]);
		double ans_quantlib = -solver.CalculationKirk();
		ans_quantlib += solver2.CalculationKirk();
		//double fd_approx = solver.CalculationFD();
		//solver.CalculationMC();

		std::cout << "\ncubature " << std::endl;
#if !sanity
		double x1 = std::log(values[i].s1);
		double x2 = std::log(values[i].s2);
		double t = values[i].t;
		double r = values[i].r;
		double sigma1 = values[i].v1;
		double sigma2 = values[i].v2;
		double rho = values[i].rho;
#else
		double x1(0), x2(0), t(.5), r(0.2), sigma1(.2), sigma2(.1), rho(0.5);
#endif
		wpUtilities::Point<dim__> strike = { values[i].strike, values[i + 1].strike };

		FunctionData data;
		BenchmarkBoundaryData boundaryData(strike);
		double err = std::pow(10., -8);
		unsigned int degree(6);
		ApproximationMorton<BenchmarkBoundaryData, FunctionData > approxBoundaryData(boundaryData, degree, err);	// wrapper function algorithm

		unsigned int DCount(0);
		double dSecondsElapsed(0), ans(0);
		unsigned int  stCubatureDegree(5);
		double e(1e-5);
		mppProgramParameters.dProportionOfRemainingTimeConsumedThisStep = 0.25;
		mppProgramParameters.dAdjustClusterDiameter = 0.75;
		mppProgramParameters.iMaxEvaluationTreeDepth = 15;
		{
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			//ans = CubaturePDESolver2D(x1, x2, strike, t, r, sigma1, sigma2, rho, mppProgramParameters, e, DCount, stCubatureDegree);
			ans = CubaturePDESolver2D(x1, x2, strike, t, r, sigma1, sigma2, rho, mppProgramParameters, e, DCount, stCubatureDegree, approxBoundaryData, data);
		}

		//double ans_exact = values[i].result; 
#if sanity
		double ans_exact = ExactSolution(x1, x2, r, sigma1, sigma2, rho, t);
#endif
		std::cout << "\nElapsed time in seconds " << dSecondsElapsed << "\n";
		std::cout << "deg5 cubature approx: " << ans << std::endl;
#if sanity
		std::cout << "exact solution      : " << ans_exact << std::endl;
		std::cout << "rel error           : " << std::abs(ans - ans_exact) / ans_exact << std::endl;
#else
		std::cout << "QuantLib approx     : " << ans_quantlib << std::endl;
		//std::cout << "QuantLib finite diff: " << fd_approx << std::endl;
		std::cout << "relerr (vs quantlib): " << std::abs(ans - ans_quantlib) / ans_quantlib << std::endl;
		//std::cout << "| cub - fd |/fd   : " << std::abs(ans - fd_approx) / fd_approx << std::endl;
		//std::cout << "| fd - kirk|/kirk : " << std::abs(ans_kirk - fd_approx) / ans_kirk << std::endl;
#endif
		std::cout << "\nTree depth: " << DCount << std::endl;
		std::cout << "No. of boundary function evaluations :" << mppProgramParameters.iCountBoundaryFunctionEvaluations << std::endl;
		std::cout << "no of extra exact function calls: " << approxBoundaryData.getNoOfExactFunctionCalls() << std::endl;
		std::cout << "no of extra approx function calls: " << approxBoundaryData.getNoOfApproxFunctionCalls() << std::endl;
		std::cout << "number of patches " << approxBoundaryData.getNoOfInterpolations() << std::endl;
		std::cout << "size of stored data " << approxBoundaryData.getSizeOfPointCloud() << std::endl;
		std::cout << "no of max depth reached " << approxBoundaryData.getmyHitBottom() << std::endl;
		std::cout << "\n\n";
	}
}

// use two step cubature for the pde
struct BenchmarkCubature
{
	LogArgLessExp mBoundary;
	CVectorWeightLocation mVlw;
	CMainProgramParameters mProgramParameters;

	BenchmarkCubature(CVectorWeightLocation vlw, CMainProgramParameters pp): mVlw(vlw), mProgramParameters(pp) 
	{
		std::cout << "constructor" << std::endl;
		mProgramParameters.iMaxEvaluationTreeDepth = 30;
	}

	double return_option_value(double x)
	{
		unsigned int DCount(0);
		double value = CubaturePDESolverWithAdaptiveApprox(x, 1., mVlw, mProgramParameters, 1.e-12, DCount, mBoundary, FunctionData());
		//return std::max(ExactPDESolution(x, 1) - 0.05, 0.);
		return std::max(value - 0.05, 0.);
	}

	double operator() (const point& pt) 
	{
		//return std::max(ExactPDESolution(pt[0], 1) - 0.05, 0.);
		return return_option_value(pt[0]);
	}

	double operator() (double pt) 
	{
		//return std::max(ExactPDESolution(pt, 1) - 0.05, 0.);
		return return_option_value(pt);
	}

	double operator() (const point& pt, FunctionData& data) 
	{
		//return std::max(ExactPDESolution(pt[0], 1) - 0.05, 0.);
		return return_option_value(pt[0]);
	}
};

#if dim_ == 1
void doubleheat(CMainProgramParameters& mppProgramParameters, double cubature_err, size_t degree, double error, bool useAA)
{
	CVectorWeightLocation vlwNormalApproxMean0Var1;
	Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());
	BenchmarkExact myfunc_exact;
	BenchmarkCubature myfunc_cubature(vlwNormalApproxMean0Var1, mppProgramParameters);
	FunctionData data;
	ApproximationMorton<BenchmarkCubature, FunctionData > funcwrapper(myfunc_cubature, degree, error);	// wrapper function algorithm
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);

		double position(mppProgramParameters.dLocation), time(mppProgramParameters.dTime), vol(mppProgramParameters.dVol), rate(mppProgramParameters.dRate), strike(mppProgramParameters.dStrike), barrier(mppProgramParameters.dBarrier);
		double dApproxValue(0), dExactValue(0);
		double dSecondsElapsed;

		unsigned int DCount(0);		// added this for tracking the depth of tree reached
		// useAA = false;
		if (!useAA) {
			std::cout << "standard cubature\n";
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			// two step cubature with closed form boundary CubaturePDESolverThreeStepOneStepAA
			dApproxValue = CubaturePDESolverWithAdaptiveApprox(position, time, vlwNormalApproxMean0Var1, mppProgramParameters, cubature_err, DCount, myfunc_exact, data);
		}
		else
		{
			std::cout << "cubature with AA\n";
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			// three step cubature with adaptive approximated boundary  CubaturePDESolverWithAdaptiveApprox
			dApproxValue = CubaturePDESolverThreeStepOneStepAA(position, time, vlwNormalApproxMean0Var1, mppProgramParameters, cubature_err, DCount, funcwrapper, data);
		}
		//bm with step function boundary
		// need to work out the exact PDE solution
		//dExactValue = ExactPDESolution(position, time);
		dExactValue = 0;

		std::cout << "At (" << position << ", " << time << ") " << std::endl;
		//std::cout << "dExactValue " << dExactValue << "\n";
		std::cout << "dApproxValu " << dApproxValue << "\n";
		//std::cout << "the value of the difference  is " << dExactValue - dApproxValue << "\n";

		std::cout << "the CubaturePDESolver evaluated the boundary function " << mppProgramParameters.iCountBoundaryFunctionEvaluations
			<< " times\n";
		std::cout << "the CubaturePDESolver created "
			<< (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.iNoGaussianCubaturePoints - 1)) << " nodes\n";
		std::cout << "Elapsed time in seconds " << dSecondsElapsed << "\nElapsed time/node in seconds "
			<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.iNoGaussianCubaturePoints - 1))
			<< std::endl;
		std::cout << "Depth of tree reached " << DCount << std::endl;

		if (useAA)
		{
			std::cout << "No. of approx func calls " << funcwrapper.getNoOfApproxFunctionCalls() << std::endl;
			std::cout << "No. of exact func calls " << funcwrapper.getNoOfExactFunctionCalls() << std::endl;
			std::cout << "No. of interpolated patches " << funcwrapper.getNoOfInterpolations() << std::endl;
			std::cout << "size of stored data " << funcwrapper.getSizeOfPointCloud() << std::endl;
			std::cout << "no of max depth reached " << funcwrapper.getmyHitBottom() << std::endl;
			//std::cout << "\nInterpolated patches\n";
			//funcwrapper.getInterpolatedPatches();
		}
	}
}
#endif

int _tmain(int argc, const  _TCHAR** argv)
{
	CMainProgramParameters mppProgramParameters;
	if (mppProgramParameters.MapCommandLineToProgramParameters(argc, argv) != 0)
		return 1;
#if dim_ == 1
	do {
		

		int cub_err(0);
		std::cout << "enter cubature error order\n";
		std::cin >> cub_err;

		std::cout << "enter dLocation\n";
		std::cin >> mppProgramParameters.dLocation;

		bool useAA;
		std::cout << "useAA? \n";
		std::cin >> useAA;

		size_t degAA(8);
		int errOrderAA(12);
		if (useAA)
		{
			std::cout << "AA degree\n";
			std::cin >> degAA;
			std::cout << "AA err order\n";
			std::cin >> errOrderAA;
		}

		doubleheat(mppProgramParameters, pow(10,-cub_err) , degAA, pow(10, -errOrderAA), useAA);
	} while (1);
#else
	std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
	std::cout.precision(17);
	TestTemplateCubatureSolver();
	//Test(mppProgramParameters);
#endif
	return 0;
}

