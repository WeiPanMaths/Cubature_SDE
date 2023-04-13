#include "TemplateCubaturePDESolver.h"	//pde solver template
#include "TerminalCondition.h"
#include "TARNTwoStepPDESolver.h"
#include "CloudWeightLocation.h"
#include <iostream>
#include "OptionData.h"


double TARNOneStepPDESolver(double spot, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters & ipIntParams, double e )
{
	OptionData tarnOptionData;
	tarnOptionData.strike = 1.;
	tarnOptionData.barrier = 5.;
	tarnOptionData.time = 0.5;
#if 1	
	e = 1.e-11;

	TARNTimeOneTC TarnTimeOneBoundary( tarnOptionData, vlwNormalApproxMean0Var1, ipIntParams,e) ;

	//e *= 100;
	PDESolverDriftBM< TARNTimeOneTC > TarnTwoStepOption( TarnTimeOneBoundary, vlwNormalApproxMean0Var1,ipIntParams, e );
	
	double ans = TarnTwoStepOption(std::log(spot));
	
	TarnTimeOneBoundary.NumFuncCall();

	ipIntParams.iCountBoundaryFunctionEvaluations = TarnTwoStepOption.getCountBoundaryFunctionEvaluations();
#else

	TARNTimeTwoExact TARNExact(tarnOptionData, ipIntParams);
	double ans = TARNExact(std::log(spot));

#endif
	return ans; 
}

# if 0
/// test interpolation of the expiry option function
double TestApproxAtExpiry(double x, unsigned int d, double int_err, CMainProgramParameters & ipIntParams )
{
	OptionData tarnOptionData;
	tarnOptionData.strike = 1.;
	tarnOptionData.barrier = 5.;
	tarnOptionData.time = 0.5;

	TARNTimeTwoExact myTarnTimeTwoExact( tarnOptionData, ipIntParams);

	Approximation<TARNTimeTwoExact> splineApprox(myTarnTimeTwoExact, d, int_err);

	double y = 0.5;
	double increment = 2./100.;
	for( int i = 1; i <= 100; i++ )
	{
		y += increment;
		splineApprox(y);
	
	}
	
	splineApprox.showNumFuncCall();
	////////////// errors /////////////////////////////////////////
	
	double maxerr = 0;
	int index = 0;
	double problem_original = 0;
	double problem_approx = 0;
	double problem_y = 0;
	y = 0.5;
	increment = 2./2000000.;
	for( int i = 1; i <= 2000000; i++ )
	{
		y += increment;
		double original = myTarnTimeTwoExact(y);
		double approx = splineApprox(y);
		double err = std::abs( (original - approx )  );

		if( err > int_err && index == 0 )
		//if( err > maxerr )
		{
			index = i;
			problem_y = y;
			problem_original = original;
			problem_approx = approx;
		}

		if (err > maxerr) 
			maxerr = err ;	
	}

	splineApprox.showSetInterpolationsSize();
	splineApprox.showNumFuncCall();
	//splineApprox.showSize();
	//return (splineApprox(x) - hockeystick(x))/ hockeystick(x);
	std::cout << "maxerr is " << maxerr << std::endl;
	std::cout << "problem index " << index << std::endl;
	std::cout << "problem y " << problem_y << std::endl;
	std::cout << "problem original " << problem_original << std::endl;
	std::cout << "problem approx " << problem_approx << std::endl;
	std::cout << "problem abs err " << std::abs( problem_original - problem_approx ) << std::endl;
	std::cout << "problem rel err " << std::abs( (problem_original - problem_approx)/problem_original ) << std::endl;
	return splineApprox(x);
}


double TestApproximationClass(double x, unsigned int d, double int_err)
{
	HockeyStickFunctionTC hockeystick;
	Approximation<HockeyStickFunctionTC> splineApprox(hockeystick, d, int_err);

	double y = -13;
	double increment = (26.)/256.;
	for( int i = 1; i <= 256; i++ )
	{
		y += increment;
		splineApprox(y);
	
	}
	
	splineApprox.showNumFuncCall();
	////////////// errors /////////////////////////////////////////
	
	double maxerr = 0;
	int index = 0;
	double problem_original = 0;
	double problem_approx = 0;
	double problem_y = 0;
	y = -16;
	increment = (16.+1.3)/2000000.;
	for( int i = 1; i <= 2000000; i++ )
	{
		y += increment;
		double original = hockeystick(y);
		double approx = splineApprox(y);
		double err = std::abs( (original - approx )  );

		if( err > int_err && index == 0 )
		//if( err > maxerr )
		{
			index = i;
			problem_y = y;
			problem_original = original;
			problem_approx = approx;
		}

		if (err > maxerr) maxerr = err ;	
	}

	splineApprox.showSetInterpolationsSize();
	splineApprox.showNumFuncCall();
	//splineApprox.showSize();
	//return (splineApprox(x) - hockeystick(x))/ hockeystick(x);
	std::cout << "maxerr is " << maxerr << std::endl;
	std::cout << "problem index " << index << std::endl;
	std::cout << "problem y " << problem_y << std::endl;
	std::cout << "problem original " << problem_original << std::endl;
	std::cout << "problem approx " << problem_approx << std::endl;
	std::cout << "problem abs err " << std::abs( problem_original - problem_approx ) << std::endl;
	std::cout << "problem rel err " << std::abs( (problem_original - problem_approx)/problem_original ) << std::endl;
	return splineApprox(0);
}
#endif
//
//double OneStepPDESolver(double x, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters & ipIntParams, double e )
//{
//	CallOptionTC calloption;
//	ipIntParams.iCountBoundaryFunctionEvaluations = 0;		// reset count value
//	PDESolverDriftBM< CallOptionTC > t0Solution( calloption, vlwNormalApproxMean0Var1, ipIntParams, e);
//	
//	double ans = t0Solution(x);
//	
//	ipIntParams.iCountBoundaryFunctionEvaluations = t0Solution.getCountBoundaryFunctionEvaluations();
//
//	return ans;
//}
//
//double TwoStepSplinePDESolver(double x, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters & ipIntParams, unsigned int d, double e, double int_err)
//{
//	// half time set to 0.5, total time is 1, so no need to set the time here for this particular experiment
//	HockeyStickFunctionTC hockeystick;
//	PDESolver< HockeyStickFunctionTC > halftimeSoln(hockeystick, vlwNormalApproxMean0Var1, ipIntParams, 1.e-11);
//
//	Approximation< PDESolver<HockeyStickFunctionTC> >  halftimeApprox(halftimeSoln, d, 1.e-10);
//	PDESolver< Approximation<PDESolver<HockeyStickFunctionTC>> >   
//		fulltimeSoln( halftimeApprox, vlwNormalApproxMean0Var1, ipIntParams, 1.e-9 );
//	
//	double ans = fulltimeSoln(x);
//	//double ans = halftimeSoln(x);
//
//	halftimeApprox.showSetInterpolationsSize();
//	halftimeApprox.showNumFuncCall();		// number of times halftimeSoln is called;
//	halftimeApprox.showSize();
//	// t0 - t1 count,  ie. count of number of spline evaluations.
//	ipIntParams.iCountBoundaryFunctionEvaluations = fulltimeSoln.getCountBoundaryFunctionEvaluations();
//
//	return ans;
//}
//
//double OneStepSplinePDESolver(double x, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters & ipIntParams, unsigned int d, double e, double int_err)
//{
//	HockeyStickFunctionTC hockeystick;
//	std::cout << "degree " << d << std::endl;
//	Approximation<HockeyStickFunctionTC> splineApprox(hockeystick, d, int_err);
//	ipIntParams.iCountBoundaryFunctionEvaluations = 0;							// reset count value
//	PDESolver< Approximation<HockeyStickFunctionTC> > t0Solution( splineApprox, vlwNormalApproxMean0Var1, ipIntParams, e);
//	
//	double ans = t0Solution(x);
//
//	splineApprox.showSetInterpolationsSize();		// show interpolated patches
//	splineApprox.showNumFuncCall();				// show number of times hockeystick was called
//	splineApprox.showSize();
//	ipIntParams.iCountBoundaryFunctionEvaluations = t0Solution.getCountBoundaryFunctionEvaluations();
//
//#if 1
//	double max_int_err(0);
//	double y = -8.01;
//	double increment = 15.6/ 2000000.;
//	double problem_y;
//	double problem_original;
//	double problem_approx;
//	double problem_index;
//	for( int i = 1; i <= 2000000; i++ )
//	{
//		y += increment;
//		double value = hockeystick(y);
//		double approx_value = splineApprox(y);
//		double temperr = std::abs( value - approx_value );
//		if( temperr > max_int_err )
//		{
//			problem_index = i;
//			problem_y = y;
//			problem_original = value;
//			problem_approx = approx_value;
//			max_int_err = temperr;
//		}
//	}
//
//	std::cout << "\n TEST INTERPOLATION " << std::endl;
//	std::cout << "problem y " << problem_y << std::endl;
//	std::cout << "problem index " << problem_index << std::endl;
//	std::cout << "Problem original " << problem_original << std::endl;
//	std::cout << "problem approx   " << problem_approx << std::endl;
//	std::cout << "max err          " <<	max_int_err << std::endl << std::endl; ;
//#endif
//	return ans;
//}