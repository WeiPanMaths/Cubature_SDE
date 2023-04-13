// binary_recursion.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "CloudWeightLocation.h" // CVectorWeightLocation
#include "SharedFunctionsAndTypes.h" // operator >>
#include "CubaturePDESolver.h" // CubaturePDESolver
#include "ExactPDESolution.h" // ExactPDESolution
#include "MainProgramParameters.h" //CMainProgramParameters
#include "DoublePrecisionGaussianQuadrature.h" //Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "BlackScholesGreeks.h"
#include "trapfpe.h"
//#include "TARNTwoStepPDESolver.h"	// test two period TARN solver

void TestTimeSteppingRatioForDifferentQuadratures(CMainProgramParameters& mppProgramParameters, CVectorWeightLocation& vlwNormalApproxMean0Var1);

int _tmain(int argc, const  _TCHAR** argv)
{
	//trapfpe();
	CMainProgramParameters mppProgramParameters;
	if (mppProgramParameters.MapCommandLineToProgramParameters(argc, argv) != 0)
		return 1;
	CVectorWeightLocation vlwNormalApproxMean0Var1;
	Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());

	// test the pde solver
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);

// computes the best error which will be used for the error specification
#if 0
		/*
		std::cout << "Error With Factor " << mppProgramParameters.dProportionOfRemainingTimeConsumedThisStep
			<< " and "
			<< mppProgramParameters.NoGaussianCubaturePoints()
			<< " Cubature points of approx integral to the true solution = "
			<< CompareCubatureAndTrueSolution
			( mppProgramParameters.dProportionOfRemainingTimeConsumedThisStep
			, mppProgramParameters.dTime
			, mppProgramParameters.dLocation
			, mppProgramParameters.dStrike
			, mppProgramParameters.dBarrier
			, mppProgramParameters.dRate
			, mppProgramParameters.dVol			
		    , vlwNormalApproxMean0Var1)
			<< std::endl;*/ //terry's original

		// binary check
		std::cout << "Error With Factor " << mppProgramParameters.dProportionOfRemainingTimeConsumedThisStep
			<< " and "
			<< mppProgramParameters.NoGaussianCubaturePoints()
			<< " Cubature points of approx integral to the true solution = "
			<<  CompareCubatureAndTrueSolution(mppProgramParameters.dProportionOfRemainingTimeConsumedThisStep, 
			mppProgramParameters.dTime, mppProgramParameters.dLocation, mppProgramParameters.dRate
			, mppProgramParameters.dVol, vlwNormalApproxMean0Var1) << std::endl;
#endif

		double position(mppProgramParameters.dLocation), time(mppProgramParameters.dTime), vol(mppProgramParameters.dVol), rate(mppProgramParameters.dRate), strike(mppProgramParameters.dStrike), barrier(mppProgramParameters.dBarrier);
		double dApproxValue(0), dExactValue(0);
		double dSecondsElapsed;

		unsigned int DCount(0);		// added this for tracking the depth of tree reached
		;
		double fullThreeStepApproxValue(0);

		{
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			//dApproxValue = CubaturePDESolver(position, strike,barrier, time, rate, vol, vlwNormalApproxMean0Var1, mppProgramParameters, 0. ); // AtExpCall
			//dApproxValue = CubaturePDESolver(position, strike,barrier, time, rate, vol, vlwNormalApproxMean0Var1, mppProgramParameters, 0., DCount ); // AtExpCall with depth tracking
			//dApproxValue = CubaturePDESolver( 0, 0, 1, 1, 0, 1, vlwNormalApproxMean0Var1, mppProgramParameters, 0., DCount ); // AtExpCall test with hard coded values

			//bm  with step function boundary
			dApproxValue = CubaturePDESolver(position,time,vlwNormalApproxMean0Var1,mppProgramParameters,1.e-10,DCount);
			//DCount = 0;
			//fullThreeStepApproxValue = CubaturePDESolverThreeStepOneStep(position,time,vlwNormalApproxMean0Var1,mppProgramParameters,0.,DCount);
		}
		
		//dExactValue = ExactPDESolution(position,strike,barrier, time, rate, vol ); // AtExpCall
		//dExactValue = ExactPDESolution( 0, 0, 1, 1, 0, 1 );

		//bm with step function boundary
		dExactValue = ExactPDESolution(position,time);

		std::cout << "At (" << position << ", " << time << ") " << std::endl;
		std::cout << "the value of CubaturePDESolver is " << dApproxValue << "\n";
		//std::cout << "the value of CubaturePDEsolverThreeStepOneStep is " << fullThreeStepApproxValue << "\n";
		std::cout << "the value of ExactPDESolution  is " << dExactValue << "\n";
		std::cout << "the value of the difference    is " << dExactValue - dApproxValue << "\n";
		std::cout << "the value of the relative difference is " << abs(dExactValue - dApproxValue) / abs(1.+dExactValue) << "\n";
		//std::cout << "the value of the difference (three step)    is " << dExactValue - fullThreeStepApproxValue << "\n";
		std::cout << "the CubaturePDESolver evaluated the boundary function " << mppProgramParameters.
			iCountBoundaryFunctionEvaluations << " times\n";
		std::cout << "the CubaturePDESolver created "
			<< (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
																		  iNoGaussianCubaturePoints -
																		  1)) << " nodes\n";
		std::cout << "Elapsed time in seconds " << dSecondsElapsed << " Elapsed time/node in seconds "
			<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
																							iNoGaussianCubaturePoints -	1))
			<< std::endl;
		std::cout << "Depth of tree reached " << DCount << std::endl;

	}
#if 0
// test moving forward

	TestTimeSteppingRatioForDifferentQuadratures(mppProgramParameters, vlwNormalApproxMean0Var1);
#endif
	return 0;
}


void TestTimeSteppingRatioForDifferentQuadratures(CMainProgramParameters& mppProgramParameters, CVectorWeightLocation& vlwNormalApproxMean0Var1)
{
	// test the normal time stepping - how quickly can one go to the boundary and preserve the accuracy -
	{
		double factor = .4;
		for (int j = 0; j != 10; j++)
		{
			double t=mppProgramParameters.dTime;
			for (int k = 1; k != 15; ++k)
			{
				//double x=0.1, ans;
				//double ans = CompareCubatureAndTrueSolution(factor, t, mppProgramParameters.dLocation, vlwNormalApproxMean0Var1);

#if 0
				// at expiry call
				double ans = CompareCubatureAndTrueSolution( factor, t, mppProgramParameters.dLocation, mppProgramParameters.dStrike, mppProgramParameters.dBarrier	, mppProgramParameters.dRate, mppProgramParameters.dVol	, vlwNormalApproxMean0Var1);
#else
				// binary
				double ans = CompareCubatureAndTrueSolution( factor, t, mppProgramParameters.dLocation, mppProgramParameters.dRate, mppProgramParameters.dVol, vlwNormalApproxMean0Var1);

#endif
				std::cout << "with "
					<< (int)vlwNormalApproxMean0Var1.size()
					<< " points the difference between Cubature and answer is "
					<< ans << std::endl;
				t *= 1. - factor;	// this will shrink the time to maturity to test when we are very close to the boundary what is the accuracy that we can get
			}
			std::cout << "Timer factor was " << factor << std::endl;
			factor = .8 * factor + .2;
		}
	}
}



// Derivative test
#if 0
double h = 0.001;
double deltaExact = BSDeltaExactSolution(mppProgramParameters);

//fwd diff
double deltaApproxFwdDiff = BSDeltaFwdDifferenceCubatureApprox(mppProgramParameters,vlwNormalApproxMean0Var1,h);
std::cout << "the Cubature approx to delta using fwd diff is " << deltaApproxFwdDiff << std::endl;
std::cout << "the exact value of delta is " << deltaExact << std::endl;
std::cout << "accuracy of Cubature fwd diff delta approx is " << deltaApproxFwdDiff - deltaExact << std::endl;

//central diff
double deltaApprox = BSDeltaCubatureApprox(mppProgramParameters,vlwNormalApproxMean0Var1,h);
std::cout << "the Cubature approx to delta using central diff is " << deltaApprox << std::endl;
std::cout << "the exact value of delta is " << deltaExact << std::endl;
std::cout << "accuracy of Cubature central diff delta approx is " << deltaApprox - deltaExact << std::endl;

// five point stencil
double deltaApproxFivePoint = BSDeltaFivePointCubatureApprox(mppProgramParameters,vlwNormalApproxMean0Var1,h);
std::cout << "the Cubature approx to delta using 5point rule is " << deltaApproxFivePoint << std::endl;
std::cout << "the exact value of delta is " << deltaExact << std::endl;
std::cout << "accuracy of Cubature 5point delta approx is " << deltaApproxFivePoint - deltaExact << std::endl;

std::cout << "Elapsed time in seconds " << dSecondsElapsed << " Elapsed time/node in seconds "
<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
					  iNoGaussianCubaturePoints -
					  1))
					  << std::endl;
#endif


// Call: dApproxValue = CubaturePDESolver(position, strike, time, rate, vol, vlwNormalApproxMean0Var1, mppProgramParameters, 0. );
// Call: dExactValue = ExactPDESolution(position,strike, time, rate, vol );
