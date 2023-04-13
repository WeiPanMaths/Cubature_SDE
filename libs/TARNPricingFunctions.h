#pragma once

#include "TARNOptionData.h"
#include "Utils.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/set.hpp>
#include <fstream>
#include <string>
#include "MortonTreeApproximation.h"

namespace TARNPricing
{
	// coupon payment
	double CouponPayment ( const point_2d & pt, const double k )
	{
		double log_s = pt[0];
		double m = pt[1];			
		double ans = (m > 0)? std::min( std::exp(log_s)-k, m ) : 0.;
		return ans;
	}

	// computes new cash remaining given coupon amount
	double CalculateAccrual ( const double & cash_remaining, const double & coupon )
	{
		double value = cash_remaining - std::max( coupon, 0. );				
		return cash_remaining - coupon;
	}
		
	// TEST: pay whatever is remaining
	double FinalPaymentFunction ( const point_2d & pt, const double k )
	{
		double log_s = pt[0];
		double m = pt[1];
		double ans = (m > 0)? m : 0.;
		return ans;
	}

	template <size_t i> struct Payoff
	{
		Payoff (Data* data) : mpData(data), nextPayoff(Payoff<i-1>(data), data->model_data.interpolation_error_tolerance) {}

		double operator() ( const point_2d & pt ) 
		{
			size_t index = mpData->option_data.numberOfPayoffs - i;		// locate the right strike and time interval length
			double K = mpData->option_data.strikes[ index ];
			double couponPayment = CouponPayment( pt, K );
			
			// setup next step
			point_2d new_pt = { pt[0], CalculateAccrual(pt[1] , couponPayment) };		// accrue cash remaining on positive side of coupon payment
			mpData->model_data.dtime = mpData->option_data.dtimes[ index ];				// set the time parameter for next step solver

			double err = mpData->model_data.error_tolerance[ index +1 ];				// TODO: NEED TO MODIFY THIS
			mpData->model_data.cubature_error_tolerance = err;
			PDESolver solver(mpData);
			return (new_pt[1] > 0)? ( couponPayment + solver.evaluate(new_pt, nextPayoff, false)) : couponPayment;
		}

		Data*		mpData;
		MortonTreeApproximation nextPayoff;
		//Payoff<i-1> nextPayoff;
	};

	template <>	struct Payoff<1>
	{
		Payoff(Data* data) : mpData(data){}

		double operator() ( const point_2d & pt)
		{
			auto & data = *mpData;
			size_t index = data.option_data.numberOfPayoffs - 1;
			double K = data.option_data.strikes[ index ];
			return FinalPaymentFunction(pt, K );							// pay whatever is remaining	// this makes testing easier
			//return CouponPayment(pt, K);
		}
		Data*	mpData;
	};

	struct PDESolver
	{
		/// serialization items ////////////////////////////////////
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & slwLocationsWeightsLeft;
			ar & dAccumulatedSolution;
			ar & DCount;
			ar & dTimeLeft;
			ar & filename;
		}

		PDESolver() : DCount(0), slwLocationsWeightsLeft(), mpData(0), dTimeLeft(0), dAccumulatedSolution(0) {}
		PDESolver(Data* pdata) : mpData(pdata), DCount(0), slwLocationsWeightsLeft(), dAccumulatedSolution(0) { dTimeLeft = mpData->model_data.dtime; }
		~PDESolver(){}
		
		Data* mpData;
		unsigned int DCount;
		CSetWeightLocation slwLocationsWeightsLeft;
		double dTimeLeft;
		double dAccumulatedSolution;
		std::string filename;
		
		void	set_filename(std::string & f)
		{
			filename = f;
		};
		void	save_computation(const PDESolver &s, const char * filename, bool save_flag =false )
		{
			if(save_flag) {
				std::ofstream ofs(filename);
				boost::archive::text_oarchive oa(ofs);
				oa << s;
			}
		}
		void	restore_computation(const char * filename)
		{
			std::ifstream ifs(filename);
			boost::archive::text_iarchive ia(ifs);
			ia >> *this;
		}
		template <typename Function>
		double	evaluate ( const point_2d & pt, Function & boundary, bool save_comp = false, bool kill_comp = false, bool from_file = false )
		{
			auto & data = *mpData;
			// defining references to make variable names clearer and easier to identify
			const CVectorWeightLocation & vlwNormalApproxMean0Var1 = data.model_data.vector_weight_location;
			CMainProgramParameters & mIpIntParams = data.model_data.program_parameters;
			double mError = data.model_data.cubature_error_tolerance;

			CWeightLocation arg;
			arg.displacement = pt[0];
			arg.probability = 1;
			if(!from_file) slwLocationsWeightsLeft.insert(arg);
			//double dTimeLeft( data.model_data.dtime );		// this needs to be taken out as well
			//double dAccumulatedSolution(0);
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

						ans += (*itOuterNormal).probability * boundary( outer_pt);

						for (CVectorWeightLocation::const_iterator itInnerNormal(vlwNormalApproxMean0Var1.begin());
							itInnerNormal != vlwNormalApproxMean0Var1.end();
							++itInnerNormal)
						{
							point_2d inner_pt = { itExistingLocation->displacement + drift1 + sd1 * (*itOuterNormal).displacement + drift2 + sd2 * (*itInnerNormal).displacement, pt[1]};
							ans1 += (*itOuterNormal).probability * (*itInnerNormal).probability * boundary(  inner_pt); // path concatenation
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
				slwCarryForward.ReducePointSetMultithreaded
					( dAccuracySpread * mIpIntParams.dAdjustClusterDiameter
					, (mIpIntParams.INumerator() * vlwNormalApproxMean0Var1.size()) / mIpIntParams.IDenominator()
					, mIpIntParams);
				slwLocationsWeightsLeft.swap(slwCarryForward);
				slwCarryForward.clear();
				DCount++;

				/////////////////////
				// save information here
				this->save_computation(*this, filename.c_str(), save_comp);
				
				/// test save and load
				if(save_comp)
				{
					std::cout << this->DCount << " " << this->slwLocationsWeightsLeft.size() << " " << this->dAccumulatedSolution << std::endl;
				}
				if(DCount == 1 && kill_comp) 
					return dAccumulatedSolution;
			}	// end while loop
			return dAccumulatedSolution * df;
		}
	};	// end struct PDESolver


}