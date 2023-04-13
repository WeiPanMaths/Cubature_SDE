#pragma once

#include <vector>
#include "MainProgramParameters.h"


//////
//		CMainProgramParameters should only provide vol 
//		rest of the parameters should be with option data
/////

namespace TARNPricing
{
	struct TARNOptionData
	{
		TARNOptionData( double r, std::vector<double> & dexp, std::vector<double> & k )
			: rate(r), dtimes(dexp)	, strikes(k) 
		{ 
			numberOfPayoffs = strikes.size();	
		}

		double					rate;
		std::vector<double>		dtimes;	
		std::vector<double>		strikes;
		size_t					numberOfPayoffs;
	};

	struct TARNModelData
	{
		TARNModelData( const CVectorWeightLocation vlw, CMainProgramParameters params, double initial_dt, 
			std::vector<double> err, unsigned int interpolation_deg, double interpolation_err )
		: vector_weight_location(vlw)
		, program_parameters(params)
		, dtime(initial_dt)
		//, cubature_error_tolerance(err)
		, error_tolerance(err)
		, interpolation_error_tolerance(interpolation_err)
		, interpolation_degree(interpolation_deg)
	{}
		
		const CVectorWeightLocation		vector_weight_location;
		CMainProgramParameters			program_parameters;
		std::vector<double>				error_tolerance;
		double							cubature_error_tolerance;	// temp value holder for each step, to be modified by calling function
		double							dtime;

		double			interpolation_error_tolerance;				// TODO: may need to turn this into vector to allow for different tolerances
		unsigned int	interpolation_degree;
	};

	// wrapper to combine both optiondata and modeldata
	struct Data
	{
		Data(TARNOptionData opt, TARNModelData model) : option_data(opt), model_data(model)	{}

		TARNOptionData	option_data;
		TARNModelData	model_data;
	};
}