#pragma once

#include <vector>
#include "CloudWeightLocation.h"
#include <boost/math/special_functions/erf.hpp>
#include "TARNOptionData.h"
#include "TARNPricingFunctions.h"
#include <string>

using boost::math::erfc;
using boost::math::erf;

//#include "AdaptiveApproximation.h"
#include "MortonTreeApproximation.h"

namespace TARNPricing
{
	struct ConstFunction
	{
		double operator()(const point_2d& pt) const
		{
			return 1.0;
		}
	};

	template <size_t no_of_tarn_steps>
	class TARN
	{
		typedef		std::vector<double>		dvec;
		Data*		mpData;

	public:

		TARN(Data * tarnData) : mpData(tarnData) {}
		~TARN(void) {}
	
		double valuation( const double log_spot, const double cash_remaining)
		{	
			point_2d pt = {log_spot, cash_remaining};
			//Payoff<no_of_tarn_steps> firstPay(mpData);
			//FunctionWrapper< Payoff<no_of_tarn_steps> > payoffWrapper(&Payoff<no_of_tarn_steps>(mpData));
			FunctionWrapper< ConstFunction > payoffWrapper(&ConstFunction());
			//FunctionBase* pfirstPay = &payoffWrapper;
			MortonTreeApproximation firstPayA(&payoffWrapper, mpData->model_data.interpolation_degree, mpData->model_data.interpolation_error_tolerance);
			mpData->model_data.cubature_error_tolerance = mpData->model_data.error_tolerance[0];
			PDESolver solver(mpData);
			std::string filename;
			filename += "data.txt";
			solver.set_filename(filename);
			//auto ans = solver.evaluate(pt, firstPay, true, true, false);
			auto ans = solver.evaluate(pt, firstPayA);
			std::cout << "exact calls : " << firstPayA.getNoOfExactFunctionCalls();
			std::cout << "appro calls : " << firstPayA.getNoOfApproxFunctionCalls();
			std::cout << "no of interp: " << firstPayA.getNoOfInterpolations();
			return ans;
		}

		/// read in data from file and continue to compute
		double valuation( const double log_spot, const double cash_remaining, std::string & filename)
		{
			point_2d pt = {log_spot, cash_remaining};
			//Payoff<no_of_tarn_steps> firstPay(mpData);
			MortonTreeApproximation firstPay(Payoff<no_of_tarn_steps>(mpData));
			mpData->model_data.cubature_error_tolerance = mpData->model_data.error_tolerance[ 0 ];
			PDESolver new_solver(mpData);
			new_solver.restore_computation(filename.c_str());
			auto ans = new_solver.evaluate(pt, firstPay, true, false, true);
			return ans;
		}
	};

}	// end namespace
