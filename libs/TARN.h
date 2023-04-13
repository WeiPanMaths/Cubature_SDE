#pragma once

#include <vector>
#include <iostream>
#include "MainProgramParameters.h"
#include "CloudWeightLocation.h"
#include <boost/math/special_functions/erf.hpp>

using boost::math::erfc;
using boost::math::erf;

#include "AdaptiveApproximation.h"


//////
//		CMainProgramParameters should only provide vol 
//		rest of the parameters should be with option data
/////
namespace TARN
{
	struct TARNOptionData
	{
		TARNOptionData( double r, std::vector<double> & dexp, std::vector<double> & k )
			: rate(r)
			, dtimes(dexp)
			, strikes(k) 
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
		TARNModelData( const CVectorWeightLocation vlw, CMainProgramParameters params, double initial_dt, std::vector<double> err, unsigned int interpolation_deg, double interpolation_err )
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
		Data(TARNOptionData opt, TARNModelData model)
		: option_data(opt)
		, model_data(model)
		//, funcCallCounts(opt.numberOfPayoffs,0)
		{}

		TARNOptionData	option_data;
		TARNModelData	model_data;

		//std::vector<size_t>		funcCallCounts;		// TEST: for testing purposes
	};

	template <size_t no_of_tarn_steps> class TARN
	{
		typedef		std::vector<double>		Doublevector;
		Data		myTarnData;

		//static size_t evaluate_inner_count;					// this is for two step tarn test
		//static size_t evaluate_outer_count;
	public:

		TARN(Data & tarnData) : myTarnData(tarnData) {}
		~TARN(void) {}
	
		/// m is cash remaining
		/// t must be less than first expiry date
		double valuation( const double log_spot, const double m )
		{	
			double ans(0);
			point_2d pt = {log_spot, m};
			Payoff<no_of_tarn_steps> firstPay;
			myTarnData.model_data.cubature_error_tolerance = myTarnData.model_data.error_tolerance[ 0 ];

#if 1 			
			ApproximationMorton< Payoff<no_of_tarn_steps>, Data > 
				firstPayApprox(firstPay, myTarnData.model_data.interpolation_degree, myTarnData.model_data.interpolation_error_tolerance);
			
			if(no_of_tarn_steps > 1)
				ans = PDESolver::evaluate(pt, myTarnData, firstPayApprox);
			else
				ans = PDESolver::evaluate(pt, myTarnData, firstPay);
#endif
			
			//ans = PDESolver::evaluate(pt, myTarnData, firstPay);

			//std::cout << "total evaluation of second period solution is " << evaluate_inner_count << std::endl;
			return ans;
		}
		unsigned int get_boundary_fn_count(size_t i) const
		{
			//return myTarnData.funcCallCounts[i];
			return 0;
		}
	
	private:
		struct PDESolver;
		
		// coupon payment
		static inline double CouponPayment ( const point_2d & pt, const double k )
		{
			double log_s = pt[0];
			double m = pt[1];			
			double ans = (m > 0)? std::min( std::exp(log_s)-k, m ) : 0.;
			return ans;
		}

		// computes new cash remaining given coupon amount
		static inline double CalculateAccrual ( const double & cash_remaining, const double & coupon )
		{
			double value = cash_remaining - std::max( coupon, 0. );				
			return cash_remaining - coupon;
		}
		
		// TEST: pay whatever is remaining
		static inline double FinalPaymentFunction ( const point_2d & pt, const double k )
		{
			double log_s = pt[0];
			double m = pt[1];
			double ans = (m > 0)? m : 0.;
			return ans;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <size_t i> struct Payoff
		{
			Payoff () 
				: nextPayoff()
				, nextApprox(nextPayoff, myTarnData.model_data.interpolation_degree, myTarnData.model_data.interpolation_error_tolerance)
			{}

			double operator() ( const point_2d & pt, Data & data ) 
			{
				size_t index = data.option_data.numberOfPayoffs - i;				// locate the right strike and time interval length
				double K = data.option_data.strikes[ index ];

				double couponPayment = CouponPayment( pt, K );
			
				// setup next step
				point_2d new_pt = { pt[0], CalculateAccrual(pt[1] , couponPayment) };		// accrue cash remaining on positive side of coupon payment										
				data.model_data.dtime = data.option_data.dtimes[ index ];					// set the time parameter for next step solver

				//return (new_pt[1] > 0)? ( couponPayment + PDESolver::evaluate(new_pt, data, nextPayoff )) : couponPayment;
				return (new_pt[1] > 0)? ( couponPayment + PDESolver::evaluate(new_pt, data, nextApprox)) : couponPayment;
			}

			Payoff<i-1> nextPayoff;
			ApproximationMorton< Payoff<i-1>, Data > nextApprox;
		};

#if 1
		template <>	struct Payoff<2>
		{
			double operator() ( const point_2d & pt, Data & data ) 
			{
				size_t index = data.option_data.numberOfPayoffs - 2;				// locate the right strike and time interval length
				double K = data.option_data.strikes[ index ];

				double couponPayment = CouponPayment( pt, K );
			
				// setup next step
				point_2d new_pt = { pt[0], CalculateAccrual(pt[1] , couponPayment) };	// accrue cash remaining on positive side of coupon payment										
				data.model_data.dtime = data.option_data.dtimes[ index ];				// set the time parameter for next step solver
				
				double err = data.model_data.error_tolerance[ index +1 ];				// TODO: NEED TO MODIFY THIS
				data.model_data.cubature_error_tolerance = err;

				OneStepPayoff finalpay;
				return (new_pt[1] > 0)? ( couponPayment + PDESolver::evaluate(new_pt, data, finalpay )) : couponPayment;
				//return (new_pt[1] > 0)? ( couponPayment + new_pt[1]) : couponPayment;
			}
		};

#else
		// TEST: 
		// this funciton is the exact solution of  E(   {min(exp(z) - K , m)}^+ ) which is the exact solution
		// of second step of TARN. Use this for debugging purposes
		static inline double OneStepExactPaymentFunction ( const point_2d & pt, const double k )
		{
			double log_s = pt[0];
			double m = pt[1];	
			double T_sqr = 0.5;
			
			double integrand1 = ( log_s - std::log(k+m) ) / std::sqrt( T_sqr * 2. );
			double integrand2 = ( T_sqr + log_s - std::log(k+m)) / std::sqrt( T_sqr * 2. );
			double ans = (m > 0)?
				.5 * ( m + m * erf( integrand1 ) - k * erfc( integrand1 )  + std::exp(.5 * T_sqr + log_s) * erfc(integrand2))
					: 0.;
			return ans;
		}

		//// TEST: using exact solution for the second time period for testing reasons
		template <>
		struct Payoff<2>
		{
			double operator() ( const point_2d & pt, Data & data ) 
			{
				//++evaluate_outer_count;
				//if(evaluate_inner_count%10000000 == 0) std::cout << "outer count " << evaluate_inner_count << std::endl;

				size_t index = data.option_data.numberOfPayoffs - 2;				// locate the right strike and time interval length
				double K = data.option_data.strikes[ index ];

				double couponPayment = CouponPayment( pt, K );
			
				// setup next step
				point_2d new_pt = { pt[0], CalculateAccrual(pt[1] , couponPayment) };	// accrue cash remaining on positive side of coupon payment										
				data.model_data.dtime = data.option_data.dtimes[ index ];				// set the time parameter for next step solver
				
				double err = data.model_data.error_tolerance[ index +1 ];			// TODO: NEED TO MODIFY THIS
				data.model_data.cubature_error_tolerance = err;

				SecondStepExactPayoff finalpay;
				return (new_pt[1] > 0)? ( couponPayment + finalpay(new_pt, data) ) : couponPayment;
				//return (new_pt[1] > 0)? ( couponPayment + new_pt[1] ) : couponPayment;
			}
		};
#endif
		template <>	struct Payoff<1>
		{
			double operator() ( const point_2d & pt, Data & data )
			{
				size_t index = data.option_data.numberOfPayoffs - 1;
				double K = data.option_data.strikes[ index ];
				return CouponPayment(pt, K );
			}
		};
		struct OneStepPayoff
		{
			double operator() ( const point_2d & pt, Data & data )
			{
#if 0
				++evaluate_inner_count;
				if(evaluate_inner_count%10000000 == 0) std::cout << "inner count " << evaluate_inner_count << std::endl;
#endif
				size_t index = data.option_data.numberOfPayoffs - 1;
				double K = data.option_data.strikes[ index ];
				return FinalPaymentFunction(pt,K);				// pay remaining coupon amount - this makes the two step easy to compute
				//return CouponPayment(pt, K );				// orginal payment function
			}

		};

		// this funciton is the exact solution of  E(  {min(exp(z) - K , m)}^+ ) which is the exact solution
		// of second step of TARN. Use this for debugging purposes
		struct SecondStepExactPayoff 
		{
			double operator() ( const point_2d & pt, Data & data )
			{
				//++evaluate_inner_count;
				size_t index = data.option_data.numberOfPayoffs - 1;
				double K = data.option_data.strikes[ index ];
				return OneStepExactPaymentFunction(pt, K );
			}

		};
		
		struct PDESolver
		{
			// TEST: added tarnstepindex for testing boundary function count
			template <typename Function>
			//static double evaluate ( size_t tarnStepIndex, const point_2d & pt, Data & data, Function & boundary )
			static double evaluate ( const point_2d & pt, Data & data, Function & boundary )
			{
				// defining references to make variable names clearer and easier to identify
				const CVectorWeightLocation & vlwNormalApproxMean0Var1 = data.model_data.vector_weight_location;
				CMainProgramParameters & mIpIntParams = data.model_data.program_parameters;
				//double & mError = data.model_data.error_tolerance[tarnStepIndex];			//TEST: different tolerance for different sections
				double mError = data.model_data.cubature_error_tolerance;

				unsigned int DCount(0);
	
				CSetWeightLocation slwLocationsWeightsLeft;
				CWeightLocation arg;
				arg.displacement = pt[0];
				arg.probability = 1;
				slwLocationsWeightsLeft.insert(arg);
				double dTimeLeft( data.model_data.dtime );
				double dAccumulatedSolution(0);
				double rate = data.option_data.rate;
				double df = std::exp( - rate * dTimeLeft);
				double vol = mIpIntParams.dVol;

				while (slwLocationsWeightsLeft.size() != 0)
				{
					// drift
					double drift  =0;// (rate - .5 * vol* vol) * dTimeLeft;
					double drift1 =0;// drift * mIpIntParams.dProportionOfRemainingTimeConsumedThisStep;
					double drift2 =0;// drift * (1. - mIpIntParams.dProportionOfRemainingTimeConsumedThisStep); 
			
					// std dev
					double sd	= vol * sqrt(dTimeLeft);
					double sd1	= vol * sqrt(mIpIntParams.dProportionOfRemainingTimeConsumedThisStep * dTimeLeft);
					double sd2	= vol * sqrt(dTimeLeft *= (1. - mIpIntParams.dProportionOfRemainingTimeConsumedThisStep));
			
					CCloudWeightLocation slwCarryForward;

					for (CSetWeightLocation::iterator itExistingLocation = slwLocationsWeightsLeft.begin();
						 itExistingLocation != slwLocationsWeightsLeft.end();
						 ++itExistingLocation)
					{
						double ans(0), ans1(0);
				
						for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
							 itOuterNormal != vlwNormalApproxMean0Var1.end();
							 ++itOuterNormal)
						{
							point_2d outer_pt = {itExistingLocation->displacement + drift + sd*(*itOuterNormal).displacement, pt[1]};

							ans += (*itOuterNormal).probability * boundary( outer_pt, data );

							for (CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
								itInnerNormal != vlwNormalApproxMean0Var1.end();
								++itInnerNormal)
							{
								///mIpIntParams.iCountBoundaryFunctionEvaluations ++;				//////////
								//++data.funcCallCounts[tarnStepIndex];								/// TEST:  
						
								point_2d inner_pt = { itExistingLocation->displacement + drift1 + sd1 * (*itOuterNormal).displacement + drift2 + sd2 * (*itInnerNormal).displacement, pt[1]};

								ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * boundary(  inner_pt, data ); // path concatenation
							}
						}
						if ((abs(ans - ans1) <= mError) || (DCount >= mIpIntParams.iMaxEvaluationTreeDepth))
						{
							dAccumulatedSolution += itExistingLocation->probability * ans1;
						}
						else
						{
							for (CVectorWeightLocation::const_iterator itOuterNormal(vlwNormalApproxMean0Var1.begin());
								 itOuterNormal != vlwNormalApproxMean0Var1.end();
								 ++itOuterNormal)
							{
								CWeightLocation lwNextWeightLocation;
								lwNextWeightLocation.displacement = itExistingLocation->displacement + drift1 + sd1 * (*itOuterNormal).displacement;
								lwNextWeightLocation.probability = itExistingLocation->probability * (*itOuterNormal).probability;
								slwCarryForward.insert(lwNextWeightLocation);
							}
						}
					}

					if (mIpIntParams.bVerbose)	
						std::cout << std::endl << dTimeLeft << ", Time to maturity " << std::endl;

					double dAccuracySpread = sd1 * (vlwNormalApproxMean0Var1.begin()->displacement - vlwNormalApproxMean0Var1.
										rbegin()->displacement);
					//slwCarryForward.ReducePointSet
					slwCarryForward.ReducePointSetMultithreaded
						( dAccuracySpread * mIpIntParams.dAdjustClusterDiameter
						, (mIpIntParams.INumerator() * vlwNormalApproxMean0Var1.size()) / mIpIntParams.IDenominator()
						, mIpIntParams);

					slwLocationsWeightsLeft.swap(slwCarryForward);
					slwCarryForward.clear();
					
					DCount++;		
				
				}	// end while loop

				return dAccumulatedSolution * df;
			}	
		};	// end struct PDESolver

	};	// end class TARN

	/*template<size_t no_of_tarn_steps>
	size_t TARN<no_of_tarn_steps>::evaluate_inner_count = 0;
	
	template<size_t no_of_tarn_steps>
	size_t TARN<no_of_tarn_steps>::evaluate_outer_count = 0;*/

}	// end namespace



//// exact solution for 1 step TARN init location 0, vol 1, t = 1, m > 0
//double ExactSolution(double m, double k)
//{
//	double root_e = std::sqrt( std::exp(1.) );
//	double log_over_root_2 = std::log(k+m) / std::sqrt(2.);
//
//
//	double ans = .5 * (root_e - k  + root_e * erf( log_over_root_2 - 1./std::sqrt(2.) ) - k * erf( log_over_root_2 )
//				+ m * erfc( log_over_root_2 ) );
//
//	return ans;
//
//}

double ExactSolution(double s, double m, double k)
{
	if (m <0) return 0.;
	else {
		double log_k_m = std::log(k+m);
		double one_over_root_2 = 1./std::sqrt(2.);

		double ans = std::exp(s + .5) * (1. + erf( (log_k_m - s - 1.)*one_over_root_2 ) );
		ans	-= k * (1. + erf( (log_k_m - s )*one_over_root_2 ) );
		ans += m * erfc( (log_k_m - s) * one_over_root_2 );
		ans *= .5;

		return ans;
	}
	
}

// for s = 0
double ThesisOneStepSolution(double m, double k)
{
	double t = .5;
	double d1 = std::log(k+m) / std::sqrt(t);
	double ans = std::exp( t / 2.) * .5 * (1. + erf( (d1 - std::sqrt(t))/std::sqrt(2.) ) );
	ans += m - (k+m)*.5 * (1. + erf( d1 /std::sqrt(2.) ) );
	return ans;
}