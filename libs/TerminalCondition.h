#pragma once
#ifndef TerminalCondition_h__
#define TerminalCondition_h__

#include "MainProgramParameters.h"
#include "TemplateCubaturePDESolver.h"	// TwoStepTARNFuncHTerminalCondition requires this
#include <functional>
#include <boost/math/special_functions/erf.hpp>
using boost::math::erfc;
// #include "Approximation.h"
#include "SharedFunctionsAndTypes.h"
#include "OptionData.h"



struct TARNTimeTwoTC : public std::unary_function<double, double>
{
	TARNTimeTwoTC( OptionData optData ) : myOptionData(optData) {}
	result_type operator() ( argument_type log_asset )
	{
		double diff = std::exp(log_asset) - myOptionData.strike;
		//return std::exp( - optData.rate * optData.time) * ( ( std::max( diff, 0. ) < optData.barrier )?  diff : 0. );
		return ( std::max( diff, 0. ) < myOptionData.barrier )?  diff : 0. ;
	}

	OptionData myOptionData;
};


struct TARNTimeTwoExact
{
	TARNTimeTwoExact( OptionData optData, CMainProgramParameters modelData )
		: strike(optData.strike)
		, barrier(optData.barrier)
		, r(modelData.dRate)
		, vol(modelData.dVol)
		, t(optData.time)
	{}

	double operator() ( double log_spot )
	{
		double diff = barrier-std::exp(log_spot);
		return (diff > 0)? AtExpCall(std::exp(log_spot), strike, diff, r, vol, t) : 0;
		//return (barrier > 0)? AtExpCall(std::exp(log_spot), strike, barrier, r, vol, t) : 0;

		//return AtExpCall(std::exp(log_spot), strike, barrier, r, vol, t);
	}

	double strike;
	double barrier;
	double r;
	double vol;
	double t;
};


struct TARNTimeOneTC : public std::unary_function<double, double>
{
	TARNTimeOneTC(OptionData optData, CVectorWeightLocation vlwNormalApproxMean0Var1, CMainProgramParameters ipIntParams, double e)
		: myOptionData(optData)
		, myTimeTwoTC(optData)
		, myTimeTwoSoln( myTimeTwoTC, vlwNormalApproxMean0Var1, ipIntParams, e/1000.)
		//, myTimeTwoSolnApprox( myTimeTwoSoln, 8, e/10.)
		, myTimeTwoExact(myOptionData, ipIntParams)
		//, myTimeTwoExactApprox( myTimeTwoExact, 8, e/10. )
	{}

	result_type operator() ( argument_type log_spot )
	{
		double diff = std::exp(log_spot) - myOptionData.strike;
		
		///////////
		// experiment using exact soln
		return  (  std::max( diff, 0.0 ) < myOptionData.barrier )? diff + myTimeTwoExact(log_spot) : 0.;
		
		//return (  std::max( diff, 0.0 ) < myOptionData.barrier )? diff + myTimeTwoExactApprox(log_spot) : 0.;


		///////////
		// experiment using cubature soln
		//myTimeTwoTC.myOptionData.barrier = myOptionData.barrier - std::max(diff, 0.0);
		//return ( 0 < myTimeTwoTC.myOptionData.barrier)? diff + myTimeTwoSolnApprox(log_spot) : 0;
		//return ( 0 < myTimeTwoTC.myOptionData.barrier)? diff + myTimeTwoSoln(log_spot) : 0;
	}

	void NumFuncCall() const
	{
		//myTimeTwoExactApprox.showNumFuncCall();
		//myTimeTwoSolnApprox.showNumFuncCall();
	}

	OptionData myOptionData;
	TARNTimeTwoTC myTimeTwoTC;
	PDESolverDriftBM<TARNTimeTwoTC>  myTimeTwoSoln;
	//Approximation< PDESolverDriftBM<TARNTimeTwoTC> >  myTimeTwoSolnApprox;
	TARNTimeTwoExact myTimeTwoExact;
	//Approximation< TARNTimeTwoExact > myTimeTwoExactApprox;
};


//struct TARNTimeOneTC : public std::unary_function<double, double>
//{
//	//Approximation< PDESolver<TARNPeriodTwoTC> >  periodTwoSolnApprox ;
//
//	result_type operator() ( argument_type x )
//	{
//		double strike = 1.;
//		double barrier = 5.;
//		double assetValue = std::exp(x);
//		double diff = assetValue - strike;
//		double barrier_diff = barrier - std::max( diff, 0.0 );
//		return ( 0 < barrier_diff )? diff + AtExpCall(assetValue, strike, barrier_diff, 0.05, 0.1, 0.5) : 0. ;
//		//return ( 0 < barrier_diff )? diff : 0. ;
//	}
//};



//// call option payoff with strike set to 1. and rate 0. so no discount factor
struct CallOptionTC : public std::unary_function<double, double>
{
	result_type operator() ( argument_type x )
	{
		return std::max ( std::exp( x ) - 1. , 0. ) * std::exp(-0.05);
	}
};


//// for testing splines, first just do call option and check
struct HockeyStickFunctionTC : public std::unary_function<double, double> 
{
	result_type operator()( argument_type x )
	{
		//return (x > optData.dStrike)? (x - optData.dStrike) : 0;
		return (x < double(0)) ? double(1) - exp(x) : double(0);
		//return (erfc(x) - exp(.25 + x) * erfc(.5 + x)) / 2.;		// halftime exact solution
		//return std::sin(std::abs(x-4./3.));
		//return std::pow(x,9);
		//return std::sin(x-4./3.);
		//return x*x;
	}
};

#endif




//
//struct TwoStepTARNFuncHTerminalCondition {
//public:
//	TwoStepTARNFuncHTerminalCondition( PDESolver<AtExpTerminalCondition> _terminalFunction ) : terminalFunction(_terminalFunction) {}
//
//	double operator() ( double, CMainProgramParameters & );
//
//private:
//	PDESolver<AtExpTerminalCondition> terminalFunction;
//};
