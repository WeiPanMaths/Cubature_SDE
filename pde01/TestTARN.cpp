// binary_recursion.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "CloudWeightLocation.h"					// CVectorWeightLocation
#include "MainProgramParameters.h"					//CMainProgramParameters
#include "DoublePrecisionGaussianQuadrature.h"		//Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "TARN.h"
#include <string>

void TestTimeSteppingRatioForDifferentQuadratures(CMainProgramParameters& mppProgramParameters, CVectorWeightLocation& vlwNormalApproxMean0Var1);

int _tmain1(int argc, const  _TCHAR** argv)
{
	//trapfpe();
	CMainProgramParameters mppProgramParameters;
	if (mppProgramParameters.MapCommandLineToProgramParameters(argc, argv) != 0)
		return 1;

	// load cubature points	//////////////////////////
	CVectorWeightLocation vlwNormalApproxMean0Var1;
	Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());

	// set parameters	////////////////////////////
	const size_t no_of_pay = 1;
	const double k = 1.;
	double log_spot = 0;
	double m = 3;
	std::vector<double> dexpiries(no_of_pay,.5);
	std::vector<double> strikes(no_of_pay,k);
	std::vector<double> cubature_err(no_of_pay, std::pow(10.,-10));
	if (no_of_pay >1) cubature_err[0] *= 1.;			// for two step tarn - modify the error tolerance for the one
	else cubature_err[0] *= 1.;
	double rate = 0;
	unsigned int interpolation_degree = 8;
	double interpolation_error = std::pow(10., -10);

	// prep solver	//////////////////////////////////
	TARN::TARNModelData modelData(vlwNormalApproxMean0Var1, mppProgramParameters, dexpiries[0], cubature_err, interpolation_degree, interpolation_error);
	TARN::TARNOptionData optionData(rate, dexpiries, strikes);
	TARN::Data data(optionData,modelData);
	

	// execute	/////////////////////////////////////
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);
		double dSecondsElapsed(0);
		double dApproxValue(0);

		TARN::TARN<no_of_pay> myTARNPricing(data);
		
		{
			CScopedWindowsPerformanceTimer timer(dSecondsElapsed);			
			dApproxValue = myTARNPricing.valuation(log_spot,m);
		}
		

		std::cout << "the value of PDESolver is "   << dApproxValue << "\n";
		std::cout << "Elapsed time in seconds "		<< dSecondsElapsed << std::endl;	

		//std::cout << "exact solution " << ThesisOneStepSolution(m,k) << std::endl;
/*		
		std::cout << "the PDESolver created "
			<< (iCount / mppProgramParameters.iNoGaussianCubaturePoints) << " nodes\n";
		std::cout << "Elapsed time in seconds " << dSecondsElapsed 
			<< " Elapsed time/node in seconds "
			<< dSecondsElapsed / (iCount / mppProgramParameters.iNoGaussianCubaturePoints)
			<< std::endl;	
*/
	}

	return 0;
}
