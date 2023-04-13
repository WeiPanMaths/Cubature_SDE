// Adaptive Approximation of BS Call as a function of spot and strike

#include "stdafx.h"
//#include "CloudWeightLocation.h" // CVectorWeightLocation
//#include "SharedFunctionsAndTypes.h" // operator >>
//#include "CubaturePDESolver.h" // CubaturePDESolver
//#include "ExactPDESolution.h" // ExactPDESolution
//#include "MainProgramParameters.h" //CMainProgramParameters
//#include "DoublePrecisionGaussianQuadrature.h" //Normal01QuadratureDoublePrecision
#include "ScopedWindowsPerformanceTimer.h"
#include "AdaptiveApproximation.h"
#include <random>
#include "csvfile.h"


int LearnBSCall(int argc, const  _TCHAR** argv);

std::string create_filename(const std::string function, const size_t degree, const int tolpow)
{
#if RMSE
	auto rmse = "rmse";
#else
	auto rmse = "relerr";
#endif
	return function + "_deg_" + std::to_string(degree) + "_tolpow_" + std::to_string(tolpow) + "_" + rmse;
}

struct Payoff
{
	double operator() (const point& pt) const
	{
		//return 0.5 * exp(-0.5*((pt[0]) * (pt[0])));
		//return -exp(pt[0]);
		return (pt[0] < double(0)) ? double(1) - exp(pt[0]) : double(0);
		//return (pt[0] < double(0)) ? double(0) : chebyshev_interval_shift(sigmoid(pt[0]), 1, 52);
		//return pt[0]; // (pt[0] < double(0)) ? 0 : pt[0];
		//return chebyshev_interval_shift(sigmoid(pt[0]), 1, 52);
	}
};


struct FinalFunction
{
	double operator() (const point_2d& pt) const
	{
		/*double log_s = pt[0];
		double m = pt[1];
		double ans = (m > 0)? m : 0.;
		return ans;*/
		//auto ans = black_scholes_boost(pt[0], pt[1], 1.0);
		
		//gaussian kernel
		//double ans = 0.5 * exp(-20. * ((pt[0] - .75) * (pt[0] - .75) + (pt[1] - .75) * (pt[1] - .75)));  // 0.5 / 0.001
		
		//double ans(0.);
		double ans(0);
		// basic piecewise polynomial ///////////////////////////////////////////////////
		//if (!(pt[0] > 0.75) && !(pt[1] > 0.75))   // [0.5, 0.75)**2
		//	ans = pt[1] * pt[1] * pt[1] + 5.* pt[0] * pt[1] * pt[1] + 3.* pt[0] * pt[0];

		//else if (!(pt[0] > 0.75) && pt[1] > 0.75) //[0.5,0.75)×[0.75,1) 
		//	ans = -3. * pt[0] * pt[0] * std::pow(pt[1] - 0.75, 3) + pt[0] * pt[1] - pt[0];

		//else if (pt[0] > 0.75 && !(pt[1] > 0.75)) // [0.75,1)×[0.5,0.75)
		//	ans = -5. * std::pow(pt[0] - .75, 4) * pt[1] + 10. * pt[0] * (pt[1]-0.75);

		//else if (pt[0] > 0.75 && pt[1] > 0.75)    // [0.75, 1)**2
		//	ans = -5. * std::pow(pt[0] - .75, 4) * pt[1] + 10. * pt[0] * (pt[1] - 0.75);

		// same piecewise polynomial (rotated anticlockwise 45 degrees) //////////////////////////////
		double x = pt[0] - 0.75;  // centre of the box is the 0.75 0.75
		double y = pt[1] - 0.75;
		if (x >= y && -x > y)   // [0.5, 0.75)**2 --> x >= y, -x > y
			ans = pt[1] * pt[1] * pt[1] + 5.* pt[0] * pt[1] * pt[1] + 3.* pt[0] * pt[0];

		else if (-x >= y && x < y) //[0.5,0.75)×[0.75,1) -->  -x >= y, x < y
			ans = -3. * pt[0] * pt[0] * std::pow(pt[1] - 0.75, 3) + pt[0] * pt[1] - pt[0];

		else if (x > y && x >= -y) // [0.75,1)×[0.5,0.75)  --> x > y,  x >= -y 
			ans = -5. * std::pow(pt[0] - .75, 4) * pt[1] + 10. * pt[0] * (pt[1]-0.75);

		else if (x <= y && x > -y )    // [0.75, 1)**2  --> x <= y,  x > -y
			ans = -5. * std::pow(pt[0] - .75, 4) * pt[1] + 10. * pt[0] * (pt[1] - 0.75);

		return ans;
	}
};


//struct CubatureBSCall {
//
//	CMainProgramParameters mppProgramParameters;
//	CVectorWeightLocation vlwNormalApproxMean0Var1;
//
//	CubatureBSCall(int argc, const  _TCHAR** argv) {
//		try{
//			mppProgramParameters.MapCommandLineToProgramParameters(argc, argv);
//		}
//		catch (...){
//			std::cout << "exception caught" << std::endl;
//		}
//		
//		Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());
//	}
//
//	double operator() (const point_2d& pt)
//	{
//		auto LinearTransform = [](double currentpoint)->double
//		{
//			return 12.*(currentpoint - 0.75);
//		};
//		auto TransformPoint = [](double currentpoint)->double
//		{
//			return exp(12.*(currentpoint - 0.75));
//		};
//
//		// computes the best error which will be used for the error specification  mppProgramParameters.dTime)
//		double position(pt[0]), time(1.), vol(0.2), rate(0.);
//		double dApproxValue(0), dExactValue(0);
//		//double dSecondsElapsed;
//		//std::cout << "params " << position << " " << strike << " " << vol << std::endl;
//		{
//			//CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
//			//bm  call
//			//dApproxValue = CubaturePDESolver(position, time, vlwNormalApproxMean0Var1, mppProgramParameters, 1.e-10, DCount);
//			// position here is log spot
//			//dApproxValue = CubaturePDESolver(LinearTransform(position), TransformPoint(pt[1]), time, rate, vol, vlwNormalApproxMean0Var1, mppProgramParameters, 1.e-10);
//		}
//		
//		//bm with step function boundary 
//		dExactValue = BlackScholesCall(TransformPoint(pt[0]), TransformPoint(pt[1]), rate, vol, time);
//
//		//std::cout << "At (" << pt[0] << ", " << pt[1] << ") " << std::endl;
//		//std::cout << "the value of CubaturePDESolver is " << dApproxValue << "\n";
//		////std::cout << "the value of CubaturePDEsolverThreeStepOneStep is " << fullThreeStepApproxValue << "\n";
//		//std::cout << "the value of ExactPDESolution  is " << dExactValue << "\n";
//		//std::cout << "the value of the difference    is " << dExactValue - dApproxValue << "\n";
//		//std::cout << "the value of the rel difference is " << std::abs(dExactValue - dApproxValue) / std::abs(dExactValue + 0.001) << "\n";
//		////std::cout << "the value of the difference (three step)    is " << dExactValue - fullThreeStepApproxValue << "\n";
//		//std::cout << "the CubaturePDESolver evaluated the boundary function " << mppProgramParameters.
//		//	iCountBoundaryFunctionEvaluations << " times\n";
//		//std::cout << "the CubaturePDESolver created "
//		//	<< (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
//		//	iNoGaussianCubaturePoints -
//		//	1)) << " nodes\n";
//		//std::cout << "Elapsed time in seconds " << dSecondsElapsed << " Elapsed time/node in seconds "
//		//	<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
//		//	iNoGaussianCubaturePoints - 1))
//		//	<< std::endl;
//
//		return dExactValue; //dApproxValue;  // 
//	}
//
//};

#define SAVE 0
#if dim_==2
template<typename F, typename Func, typename D>
int TEST_APPROX(F& test, Func& myfunc, D& data, const size_t no_of_pts, const size_t degree, const int tol)
/*
test -- wrapper function
myfunc -- actual function
*/
{
	std::mt19937 Engine; // Mersenne twister MT19937
	std::uniform_real_distribution<double> UniformDistribution(0.5, 1.);

	auto LogSpotGenerator = std::bind(UniformDistribution, Engine);
	auto CreateNewPoint = [&LogSpotGenerator](point_2d& currentpoint)->void
	{
		currentpoint[0] = LogSpotGenerator();
		currentpoint[1] = LogSpotGenerator();
	};

	std::vector<point_2d> interpolation_pts(no_of_pts);
	for (auto it = std::begin(interpolation_pts); it != std::end(interpolation_pts); ++it)
		CreateNewPoint(*it);

	std::cout << "... learning ... " << std::endl;

	double dSeconds(0);
	{
		CScopedWindowsPerformanceTimer timer(dSeconds);
		for (size_t i = 0; i < no_of_pts; ++i)
			test(interpolation_pts[i], data);
	}
	std::cout << "learning time " << dSeconds << std::endl;

	std::cout << "no of exact function calls: " << test.getNoOfExactFunctionCalls() << std::endl;
	std::cout << "no of approx function calls: " << test.getNoOfApproxFunctionCalls() << std::endl << std::endl;
	std::cout << "number of patches " << test.getNoOfInterpolations() << std::endl;

#if SAVE
	try
	{
		csvfile csv(create_filename("learning", degree, tol) + ".csv", ","); // throws exceptions!
		csv << test.getNoOfExactFunctionCalls() << test.getNoOfApproxFunctionCalls() << test.getNoOfInterpolations() 
			<< test.getSizeOfPointCloud() << test.getmyHitBottom() << endrow;
	}
	catch (const std::exception &ex)
	{
		std::cout << "Exception was thrown: " << ex.what() << std::endl;
	}
#endif 
	// evaluate errors
	test.resetNoOfExactFunctionCalls();
	test.resetNoOfApproxFunctionCalls();

	if (1) 
	{
		std::cout << "..... evaluating new test set ... " << std::endl;
		point_2d maxerrpt;
		double max_exact;
		double max_approx;
		double rmse(0);
		try
		{
			double max_rel_err = 0;
			double max_abs_err = 0;

			auto save_data = [&maxerrpt, &max_exact, &max_approx](point_2d newpt, double exactvalue, double approxvalue)->void
			{
				maxerrpt[0] = newpt[0];
				maxerrpt[1] = newpt[1];
				max_exact = exactvalue;
				max_approx = approxvalue;
			};
#if 0
			csvfile csv(create_filename("results", degree, tol) + ".csv", ","); // throws exceptions!
#endif
			for (size_t i = 0; i < no_of_pts; ++i)
			{
				point_2d newpt;
				CreateNewPoint(newpt);
				double exactvalues = myfunc(newpt);
				double approxvalues = test(newpt, data);
				double abstemp = std::abs( exactvalues - approxvalues );		// works if the function is nonzero
				double reltemp = abstemp / std::abs(exactvalues + 1.e-15);
				//std::cout << newpt[0] << ", " << newpt[1] << ", "<< abstemp << std::endl;
				rmse += abstemp*abstemp;

				bool absflag = false;
				if ((abstemp > max_abs_err) && absflag)
				{
					max_abs_err = abstemp;
					save_data(newpt, exactvalues, approxvalues);
				}
				if ((reltemp > max_rel_err) && !absflag)
				{
					max_rel_err = reltemp;
					save_data(newpt, exactvalues, approxvalues);
				}
#if 0
				if (i % 30 == 0)
					csv << newpt[0] << newpt[1] << exactvalues << approxvalues << endrow;
#endif
			}
		}
		catch (const std::exception &ex)
		{
			std::cout << "Exception was thrown: " << ex.what() << std::endl;
		}

		rmse = std::sqrt(rmse / (double)no_of_pts);
		std::cout << "rmse: " << rmse << std::endl;
		//std::cout << "max abs error: " << std::abs(max_exact - max_approx) << std::endl;
		std::cout << "max rel error: " << std::abs(max_exact - max_approx) / std::abs(max_exact + 1.e-15)  << std::endl;
		std::cout << "exact value used in computing max abs err = " << max_exact << std::endl;
		std::cout << "approx value used in computing max abs err = " << max_approx << std::endl << std::endl;
		std::cout << "no of extra exact function calls: " << test.getNoOfExactFunctionCalls() << std::endl;
		std::cout << "no of extra approx function calls: " << test.getNoOfApproxFunctionCalls() << std::endl;
		std::cout << "number of patches " << test.getNoOfInterpolations() << std::endl;
		std::cout << "size of stored data " << test.getSizeOfPointCloud() << std::endl;
		std::cout << "no of max depth reached " << test.getmyHitBottom() << std::endl;

		std::vector<double> patches;
		test.getInterpolatedPatchesCoeffs(patches);

#if SAVE
		try
		{
			csvfile csv(create_filename("evaluating", degree, tol) + ".csv", ","); // throws exceptions!
			csv << test.getNoOfExactFunctionCalls() << test.getNoOfApproxFunctionCalls() << test.getNoOfInterpolations()
				<< std::abs(max_exact - max_approx) << std::abs(max_exact - max_approx) / std::abs(max_exact + 1.e-15) << rmse << endrow;

			csvfile csv_patch(create_filename("patches", degree, tol) + ".csv", ",");
			int i(0);
			int num_coeffs((2 + degree)*(1 + degree) / 2);
			for (auto itr = patches.begin(); itr != patches.end(); ++itr, ++i) {
				csv_patch << *itr;
				if ((i + 1) % (4+num_coeffs) == 0)
					csv_patch << endrow;
			}
		}
		catch (const std::exception &ex)
		{
			std::cout << "Exception was thrown: " << ex.what() << std::endl;
		}
#endif
	}

	return 0;
}
#endif


#if dim_ == 1
template<typename F, typename Func, typename D>
int TEST_1D_APPROX(F& test, Func& myfunc, D& data, const size_t no_of_pts, const size_t degree, const int tol)
/*
test -- wrapper function
myfunc -- actual functionxit
*/
{
	double lo_range = 6.;
	std::mt19937 Engine; // Mersenne twister MT19937
	std::uniform_real_distribution<double> UniformDistribution(-lo_range, lo_range);

	auto LogSpotGenerator = std::bind(UniformDistribution, Engine);
	auto CreateNewPoint = [&LogSpotGenerator](point& currentpoint)->void
	{
		currentpoint[0] = LogSpotGenerator();
	};

	std::vector<point> interpolation_pts(no_of_pts);
	for (auto it = std::begin(interpolation_pts); it != std::end(interpolation_pts); ++it)
		CreateNewPoint(*it);

	std::cout << "... learning ... " << std::endl;

	double dSeconds(0);
	{
		CScopedWindowsPerformanceTimer timer(dSeconds);
		for (size_t i = 0; i < no_of_pts; ++i)
			test(interpolation_pts[i], data);
	}
	std::cout << "learning time " << dSeconds << std::endl;

	std::cout << "no of exact function calls: " << test.getNoOfExactFunctionCalls() << std::endl;
	std::cout << "no of approx function calls: " << test.getNoOfApproxFunctionCalls() << std::endl << std::endl;
	std::cout << "number of patches " << test.getNoOfInterpolations() << std::endl;


	// evaluate errors
	test.resetNoOfExactFunctionCalls();
	test.resetNoOfApproxFunctionCalls();

	if (1)
	{
		std::cout << "..... evaluating new test set ... " << std::endl;
		point maxerrpt;
		double max_exact;
		double max_approx;
		double rmse(0);
		try
		{
			double max_rel_err = 0;
			double max_abs_err = 0;

			auto save_data = [&maxerrpt, &max_exact, &max_approx](point newpt, double exactvalue, double approxvalue)->void
			{
				maxerrpt[0] = newpt[0];
				max_exact = exactvalue;
				max_approx = approxvalue;
			};

			for (size_t i = 0; i < no_of_pts; ++i)
			{
				point newpt;
				CreateNewPoint(newpt);
				double exactvalues = myfunc(newpt);
				double approxvalues = test(newpt, data);
				double abstemp = std::abs(exactvalues - approxvalues);		// works if the function is nonzero
				double reltemp = abstemp / std::abs(exactvalues + 1.e-15);
				//std::cout << newpt[0] << ", " << newpt[1] << ", "<< abstemp << std::endl;
				rmse += abstemp * abstemp;

				bool absflag = false;
				if ((abstemp > max_abs_err) && absflag)
				{
					max_abs_err = abstemp;
					save_data(newpt, exactvalues, approxvalues);
				}
				if ((reltemp > max_rel_err) && !absflag)
				{
					max_rel_err = reltemp;
					save_data(newpt, exactvalues, approxvalues);
				}
			}
		}
		catch (const std::exception & ex)
		{
			std::cout << "Exception was thrown: " << ex.what() << std::endl;
		}

		rmse = std::sqrt(rmse / (double)no_of_pts);
		std::cout << "rmse: " << rmse << std::endl;
		std::cout << "max abs error: " << std::abs(max_exact - max_approx) << std::endl;
		std::cout << "max rel error: " << std::abs(max_exact - max_approx) / std::abs(max_exact + 1.e-15) << std::endl;
		std::cout << "exact value used in computing max abs err = " << max_exact << std::endl;
		std::cout << "approx value used in computing max abs err = " << max_approx << std::endl << std::endl;
		std::cout << "no of extra exact function calls: " << test.getNoOfExactFunctionCalls() << std::endl;
		std::cout << "no of extra approx function calls: " << test.getNoOfApproxFunctionCalls() << std::endl;
		std::cout << "number of patches " << test.getNoOfInterpolations() << std::endl;
		std::cout << "size of stored data " << test.getSizeOfPointCloud() << std::endl;
		std::cout << "no of max depth reached " << test.getmyHitBottom() << std::endl;
	}
	std::cout << "\nInterpolated patches\n";
	test.getInterpolatedPatches();
	//std::cout << "\nAll patches\n";
	//test.getAllPatches();

	return 0;
}
#endif 

struct FunctionData{};

#if dim_ ==2
//void testBSCallApprox(const unsigned int degree, const unsigned int no_of_points, int err_pow, int argc, const  _TCHAR** argv) {
//	FunctionWrapper<CubatureBSCall, point_2d> myfunc(&CubatureBSCall(argc, argv));
//	FunctionData data;
//	double err = std::pow(10., -err_pow);
//
//	ApproximationMorton<FunctionBase, FunctionData > test2(myfunc, degree, err);
//	std::cout << "\n ---- deg " << degree << ", err tol pow " << err_pow << " ---- \n";
//	double dSecondsElapsed(0);
//	{
//		CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
//		TEST_APPROX(test2, myfunc, data, no_of_points, degree, err_pow);
//	}
//	//test2.TEST_print_subdivide();
//	std::cout << "-------total elapsed time in seconds -------------\n" << dSecondsElapsed << std::endl;
//}

void testGaussianKernelApprox(const unsigned int degree, const unsigned int no_of_points, int err_pow) {
	FunctionWrapper<FinalFunction, point_2d> myfunc(&FinalFunction());
	FunctionData data;
	double err = std::pow(10., -err_pow);

	ApproximationMorton<FunctionBase, FunctionData > test2(myfunc, degree, err);	// wrapper function algorithm
	std::cout << "---- degree " << degree << ", err tol pow " << err_pow << " ---- \n";
	double dSecondsElapsed = 0;
	{
		CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
		TEST_APPROX(test2, myfunc, data, no_of_points, degree, err_pow);
	}
	//test2.TEST_print_subdivide();
	std::cout << "-------total elapsed time in seconds -------------\n" << dSecondsElapsed << std::endl;
}
#endif

#if dim_ == 1
void test1Dhockystick(const unsigned int degree, const unsigned int no_of_points, int err_pow)
{
	FunctionWrapper<Payoff, point> myfunc(&Payoff());
	FunctionData data;
	double err = std::pow(10., -err_pow);

	ApproximationMorton<FunctionBase, FunctionData > test2(myfunc, degree, err);	// wrapper function algorithm
	std::cout << "---- degree " << degree << ", err tol pow " << err_pow << " ---- \n";
	double dSecondsElapsed = 0;
	{
		CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
		TEST_1D_APPROX(test2, myfunc, data, no_of_points, degree, err_pow);
	}
	std::cout << "-------total elapsed time in seconds -------------\n" << dSecondsElapsed << std::endl;
}
#endif

int _tmain_(int argc, const  _TCHAR** argv)
{
	std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
	std::cout.precision(15); 

	//point_2d pt = { 0., 1. };
	//LearnBSCall(argc, argv);
	//CubatureBSCall myCall(argc, argv);
	//myCall(pt);
	/*unsigned int degrees[] = { std::strtol(argv[1], NULL, 0) };
	int tols[] = { std::strtol(argv[2], NULL, 0) };*/
	int _temp_;
	std::cout << "press any key to start ";
	std::cin >> _temp_;
	unsigned int degrees[] = { _temp_ }; // { 3, 4, 5, 6 }; //, 5, 6, 7}; //{ 3 }; // { 5 }; 
	int tols[] = { 6 }; // , 7, 8, 9}; //{ 10 }; // 
	const unsigned int no_of_points = 3000000;

	for (auto itdeg = std::begin(degrees); itdeg != std::end(degrees); itdeg++)
		for (auto ittols = std::begin(tols); ittols != std::end(tols); ittols++)
				//test1Dhockystick(*itdeg, no_of_points, *ittols);
			testGaussianKernelApprox(*itdeg, no_of_points, *ittols);
			//testBSCallApprox(*itdeg, no_of_points, *ittols, argc, argv);
	
	//testCubature(argc, argv);
	return 0;
}

# if 0

void testCubature(int argc, const  _TCHAR** argv)
{
	CubatureBSCall solver(argc, argv);
	std::mt19937 Engine; // Mersenne twister MT19937
	std::uniform_real_distribution<double> UniformDistribution(0.5, 1.);
	size_t no_of_points = 5;

	auto LogSpotGenerator = std::bind(UniformDistribution, Engine);
	auto CreateNewPoint = [&LogSpotGenerator](point_2d& currentpoint)->void
	{
		currentpoint[0] = LogSpotGenerator();
		currentpoint[1] = LogSpotGenerator();
	};

	std::vector<point_2d> interpolation_pts(no_of_points);
	for (auto it = std::begin(interpolation_pts); it != std::end(interpolation_pts); ++it)
		CreateNewPoint(*it);

	std::cout << "\n\n";
	for (size_t i = 0; i < no_of_points; ++i)
		solver(interpolation_pts[i]);
}


int LearnBSCall(int argc, const  _TCHAR** argv)
{
	CMainProgramParameters mppProgramParameters;
	if (mppProgramParameters.MapCommandLineToProgramParameters(argc, argv) != 0)
		return 1;
	CVectorWeightLocation vlwNormalApproxMean0Var1;
	Normal01QuadratureDoublePrecision(vlwNormalApproxMean0Var1, mppProgramParameters.NoGaussianCubaturePoints());

	// test the pde solver
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(15);

		// computes the best error which will be used for the error specification  mppProgramParameters.dTime)
		double position(mppProgramParameters.dLocation), time(0.5), vol(0.2), //vol(mppProgramParameters.dVol), 
			rate(0.), //rate(mppProgramParameters.dRate), 
			strike(mppProgramParameters.dStrike),
			barrier(mppProgramParameters.dBarrier);
		double dApproxValue(0), dExactValue(0);
		double dSecondsElapsed;
		std::cout << "params " << position << " " << strike << " " << vol << std::endl;
		{
			//CScopedWindowsPerformanceTimer timer(dSecondsElapsed);
			//bm  call
			//dApproxValue = CubaturePDESolver(position, time, vlwNormalApproxMean0Var1, mppProgramParameters, 1.e-10, DCount);
			// position here is log spot
			dApproxValue = CubaturePDESolver(position, strike, time, rate, vol,
				vlwNormalApproxMean0Var1,
				mppProgramParameters, 1.e-10);
		}

		//bm with step function boundary
		dExactValue = BlackScholesCall(exp(position), strike, rate, vol, time);

		std::cout << "At (" << position << ", " << time << ") " << std::endl;
		std::cout << "the value of CubaturePDESolver is " << dApproxValue << "\n";
		//std::cout << "the value of CubaturePDEsolverThreeStepOneStep is " << fullThreeStepApproxValue << "\n";
		std::cout << "the value of ExactPDESolution  is " << dExactValue << "\n";
		std::cout << "the value of the difference    is " << dExactValue - dApproxValue << "\n";
		//std::cout << "the value of the difference (three step)    is " << dExactValue - fullThreeStepApproxValue << "\n";
		std::cout << "the CubaturePDESolver evaluated the boundary function " << mppProgramParameters.
			iCountBoundaryFunctionEvaluations << " times\n";
		std::cout << "the CubaturePDESolver created "
			<< (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
			iNoGaussianCubaturePoints -
			1)) << " nodes\n";
		std::cout << "Elapsed time in seconds " << dSecondsElapsed << " Elapsed time/node in seconds "
			<< dSecondsElapsed / (mppProgramParameters.iCountBoundaryFunctionEvaluations / (mppProgramParameters.
			iNoGaussianCubaturePoints - 1))
			<< std::endl;

	}

	return 0;
}
#endif