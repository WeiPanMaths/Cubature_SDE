// binary_recursion.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "CloudWeightLocation.h" // CVectorWeightLocation
#include "SharedFunctionsAndTypes.h" // operator >>
#include "CubaturePDESolver.h" // CubaturePDESolver
#include "ExactPDESolution.h" // ExactPDESolution
#include "MainProgramParameters.h" //CMainProgramParameters
#include "DoublePrecisionGaussianQuadrature.h" //Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "trapfpe.h"


void benchmark(CMainProgramParameters& mppProgramParameters, size_t degree, double error, bool useAA)
{
	CVectorWeightLocation vlwNormalApproxMean0Var1;
	Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());
	LogArgLessExp myfunc;
	FunctionData data;
	ApproximationMorton<LogArgLessExp, FunctionData > funcwrapper(myfunc, degree, error);	// wrapper function algorithm
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);

		double position(mppProgramParameters.dLocation), time(mppProgramParameters.dTime), vol(mppProgramParameters.dVol), rate(mppProgramParameters.dRate), strike(mppProgramParameters.dStrike), barrier(mppProgramParameters.dBarrier);
		double dApproxValue(0), dExactValue(0);
		double dSecondsElapsed;

		unsigned int DCount(0);		// added this for tracking the depth of tree reached

		if (!useAA) {
			std::cout << "standard cubature\n";
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			//bm  with step function boundary
			// CubaturePDESolverWithAdaptiveApprox
			dApproxValue = CubaturePDESolverThreeStepOneStepAA(position, time, vlwNormalApproxMean0Var1, mppProgramParameters, 1.e-10, DCount, myfunc, data);
		}
		else
		{
			//CubaturePDESolverThreeStepOneStepAA
			std::cout << "cubature with AA\n";
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			dApproxValue = CubaturePDESolverThreeStepOneStepAA(position, time, vlwNormalApproxMean0Var1,mppProgramParameters, 1.e-10, DCount, funcwrapper, data);
		}
		//bm with step function boundary
		dExactValue = ExactPDESolution(position, time);

		std::cout << "At (" << position << ", " << time << ") " << std::endl;
		std::cout << "dExactValue " << dExactValue << "\n";
		std::cout << "dApproxValu " << dApproxValue << "\n";
		std::cout << "the value of the difference  is " << dExactValue - dApproxValue << "\n";

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
			std::cout << "\nInterpolated patches\n";
			funcwrapper.getInterpolatedPatches();
		}
	}
}

int _tmain(int argc, const  _TCHAR** argv)
{
	do {
		CMainProgramParameters mppProgramParameters;
		if (mppProgramParameters.MapCommandLineToProgramParameters(argc, argv) != 0)
			return 1;

		std::cout << "enter dLocation\n";
		std::cin >> mppProgramParameters.dLocation;

		bool useAA;
		std::cout << "useAA? \n";
		std::cin >> useAA;

		size_t degAA(8);
		int errOrderAA(10);
		if (useAA)
		{
			std::cout << "AA degree\n";
			std::cin >> degAA;
			std::cout << "AA err order\n";
			std::cin >> errOrderAA;
		}

		benchmark(mppProgramParameters, degAA, pow(10,-errOrderAA), useAA);
	} while (1);
	
	return 0;
}

