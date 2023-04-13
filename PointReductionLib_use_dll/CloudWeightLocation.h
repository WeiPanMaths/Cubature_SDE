/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once
#ifndef CloudWeightLocation_h__
#define CloudWeightLocation_h__

#include "TreeBufferHelper.h"// CTreeBufferHelper
#include "MainProgramParameters.h"// CMainProgramParameters
//#include "SharedFunctionsAndTypes.h"// CWeightLocation CWeightLocationCompare
#include <valarray> // std::valarray<double>
#include <set> // std::set< CWeightLocation, CWeightLocationCompare > BASE
#include <map>
#include <vector>// std::vector< BASE::iterator >
//#include "Cmove1.h"
//#include <stdlib.h>
#include "recombine.h"
#include <boost/thread/thread.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


struct CWeightLocation
{
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & probability;
		ar & displacement;
	}

	mutable double probability;
	double displacement;
};
template <size_t depth> // depth >=1
struct CWeightMoments : public std::valarray<double>
{
	CWeightMoments(const CWeightLocation& rhs =  CWeightLocation())
		: std::valarray<double> (depth + 1)
	{
		(*this)[0] = rhs.probability;
		(*this)[1] = rhs.displacement;
		for (int i = 2; i < depth + 1; )
			(*this)[i] = pow(rhs.displacement, i);
	}
};


struct CWeightLocationPlus : public std::binary_function<CWeightLocation, CWeightLocation, CWeightLocation>
{	// functor for operator+
	CWeightLocation operator()(const CWeightLocation& _Left, const CWeightLocation& _Right) const
	{	// apply operator+ to operands
		CWeightLocation answer;
		answer.displacement = _Left.displacement + _Right.displacement;
		answer.probability = _Left.probability + _Right.probability;
		return answer;
	}
};
struct CWeightLocationCompare : public std::binary_function<bool, CWeightLocation, CWeightLocation>
{	// functor for operator+
	bool operator()(const CWeightLocation& _Left, const CWeightLocation& _Right) const
	{	// apply operator< to operands
		return
			(_Left.displacement < _Right.displacement)
			? true
			: (_Left.displacement == _Right.displacement && _Left.probability < _Right.probability)
			? true
			: false;
	}
};
typedef std::vector<CWeightLocation> CVectorWeightLocation;
typedef std::set<CWeightLocation, CWeightLocationCompare> CSetWeightLocation;

class CCloudWeightLocation :
public std::set<CWeightLocation, CWeightLocationCompare>
{
	typedef std::set<CWeightLocation, CWeightLocationCompare> BASE;
public:
	CCloudWeightLocation(void);
	~CCloudWeightLocation(void);

	// reduces the pointset based on reducing the number of points in each interval of given width to the accuracy of the cubature
	CCloudWeightLocation& ReducePointSet(const double& dAccuracySpread, const unsigned int& iNoDimensionsToCubature,
										 CMainProgramParameters& ipIntParams);

	CCloudWeightLocation& ReducePointSetMultithreaded(const double& dAccuracySpread, const unsigned int& iNoDimensionsToCubature,
										 CMainProgramParameters& ipIntParams);

	// breaks the cloud up into components // depends on the model
	void PartitionCloud(const double& dAccuracySpread, std::vector<std::vector<BASE::iterator> >&
						vvplwPartitionedCloud);

	void DoAbstractPruning(std::vector<BASE::iterator>& vplwToBeProcessed, const unsigned int& iNoDimensionsToCubature,
						   CMainProgramParameters& ipIntParams);

	void DoAbstractPruningMultithreaded(std::vector<std::set<CWeightLocation, CWeightLocationCompare>::iterator>& vplwToBeProcessed,
							const unsigned int& iNoDimensionsToCubature, CMainProgramParameters& ipIntParams, 
							std::vector<double>& NewWeights, std::map<size_t, size_t> & miIn2Out);

	void DoAbstractPruningUpdate(std::vector<BASE::iterator>& vplwToBeProcessed,
							std::vector<double>& NewWeights, std::map<size_t, size_t> & miIn2Out,
							size_t & noActiveWeightsLocations);

	static unsigned int IdentifyRebalancing(CTreeBufferHelper& bhBufInfo, std::vector<double>& vdWeightsBuffer,
											std::vector<BASE::iterator>& vplwToBeProcessed,
											std::vector<std::valarray<double> >& vdPointsBuffer,
											unsigned int& ICountCalls, std::map<size_t, size_t>& miTreePosition,
											std::vector<double>& weights);

	// remove locations and change weights in the cloud
	void UpdateCloud(CTreeBufferHelper& bhBufInfo, std::map<size_t, size_t>& miTreePosition,
					 std::vector<BASE::iterator>& vplwToBeProcessed, std::vector<double>& weights);
	//static unsigned int IdentifyLocationsRemainingAndTheirNewWeights( 
	//	CTreeBufferHelper &bhBufInfo, 
	//	std::map<size_t , size_t > &miTreePosition, 
	//	std::vector<double> &vdWeightsBuffer, 
	//	std::vector< std::valarray<double> > &vdPointsBuffer,
	//	std::vector < double >& weights, 
	//	unsigned int &ICountCalls 
	//	){exit(1);return 0;}//dummy definition
	// using the tree data in bhBufInfo and using the first bhBufInfo.iInitialNoLeaves
	// in the buffer as leaf values populates nodes in teh tree with the 
	// sums recursively to produce a tree with the specific number of roots
	static void ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
															   std::vector<double>& vdWeightsBuffer,
															   std::vector<std::valarray<double> >& vdPointsBuffer);

	static void ExpandRawLocationWeightsToConditionedVectorsOfDoubles(const CTreeBufferHelper& bhBufInfo,
																	  std::vector<double>& vdWeightsBuffer,
																	  const std::vector<BASE::iterator>&
																	  vplwToBeProcessed,
																	  std::vector<std::valarray<double> >&
																	  vdPointsBuffer);

	static void COGVar(CWeightLocation& wlCOG, const CTreeBufferHelper& bhBufInfo,
					   std::vector<double>& vdWeightsBuffer, const std::vector<BASE::iterator>& vplwToBeProcessed,
					   double& dApproxStandardDeviation);

	static void LocationToVector(const size_t iDepthOfVector, std::valarray<double>& dPointsBuffer,
								 const std::set<CWeightLocation, CWeightLocationCompare>::iterator& plwToBeProcessed,
								 CWeightLocation& wlCOG, double& dApproxStandardDeviation);

/*	class ToPVoid
	{
	public:
		// The function call
		void* operator ()(BASE::iterator& plw) const
		{
			return &(plw->displacement);
		}
	};
*/
	class ToWeight
	{
	public:
		// The function call
		double operator ()(BASE::iterator& plw) const
		{
			return (plw->probability);
		}
	};
};


#endif // CloudWeightLocation_h__
