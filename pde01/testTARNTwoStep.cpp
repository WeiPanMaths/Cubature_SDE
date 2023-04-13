// binary_recursion.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "ExactPDESolution.h" // ExactPDESolution
#include "CloudWeightLocation.h" // CVectorWeightLocation
//#include "SharedFunctionsAndTypes.h" // operator >>
#include "MainProgramParameters.h" //CMainProgramParameters
#include "DoublePrecisionGaussianQuadrature.h" //Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "TARNTwoStepPDESolver.h"	// test two period TARN solver
#include "TemplateCubaturePDESolver.h"
#include <fstream>
#include <sstream>

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

		/*std::ofstream myfile;
		std::stringstream filename;
		filename << "experiment_d" << d << ".txt";
		myfile.open (filename.str());*/
		
	
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);

		double spot(mppProgramParameters.dLocation), time(mppProgramParameters.dTime), vol(mppProgramParameters.dVol), rate(mppProgramParameters.dRate), strike(mppProgramParameters.dStrike), barrier(mppProgramParameters.dBarrier);
		double dApproxValue(0), dExact(0);
		double dSecondsElapsed;
		
		spot = 4;
		unsigned int d = 16;
		int noSplineCountBFuncEvaluations(0), splineCountBFuncEvaluations(0);
		{
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);			
			int cub_err_exp(10);
			int int_err_exp(10);
			double cub_err = std::pow(10.0, -cub_err_exp);
			double int_err = std::pow(10.0, -int_err_exp);
			
			dApproxValue = TARNOneStepPDESolver(spot,vlwNormalApproxMean0Var1, mppProgramParameters, cub_err );
			
			//TestApproxAtExpiry(std::log(spot), d, int_err, mppProgramParameters);
			//dApproxValue = OneStepSplinePDESolver(position,vlwNormalApproxMean0Var1, mppProgramParameters, d, cub_err, int_err );
			//dApproxValue = TwoStepSplinePDESolver(position,vlwNormalApproxMean0Var1, mppProgramParameters, d, cub_err, int_err );
		}
		
		dExact =  0;//ExactPDESolution(position,1.,5.,0.5,0.05,.1);

		std::cout << "position " << spot << " time " << time << std::endl;

		std::cout << "At (" << spot << ", " << time << ") " << std::endl;
		std::cout << "the value of PDESolver is " << dApproxValue << "\n";
		std::cout << "the value of exact is " << dExact << "\n";
		std::cout << "difference between exact and approx is " << std::abs(dApproxValue - dExact) << "\n";
		
		std::cout << "the pde solver spline evaluated bound function " << mppProgramParameters.iCountBoundaryFunctionEvaluations << " times\n";
		
		std::cout << "the PDESolver created "
			<< (mppProgramParameters.iCountBoundaryFunctionEvaluations / mppProgramParameters.iNoGaussianCubaturePoints) << " nodes\n";
		std::cout << "Elapsed time in seconds " << dSecondsElapsed 
			<< " Elapsed time/node in seconds "
			<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / mppProgramParameters.iNoGaussianCubaturePoints)
			<< std::endl;	


		//myfile << "cub_err " << cub_err << "\n";
		//myfile << "int_err " << int_err << "\n";	
		///*myfile << "the PDESolver with spline evaluated the boundary function " << splineCountBFuncEvaluations << " times\n";
		//myfile << "the PDESolver w/o  spline evaluated the boundary function " << noSplineCountBFuncEvaluations << " times\n";*/
		//myfile << "difference between number of count b fun evaluations is " << splineCountBFuncEvaluations - noSplineCountBFuncEvaluations << " times\n";
		//myfile << "absolute difference between the two solutions " << std::abs( dApproxValue1Step - dApproxValue1StepSpline ) << "\n\n";
		//myfile.close();

	}

	return 0;
}
