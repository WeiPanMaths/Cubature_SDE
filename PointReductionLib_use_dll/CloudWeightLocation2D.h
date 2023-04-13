#pragma once
#ifndef CloudWeightLocation2D_h__
#define CloudWeightLocation2D_h__

#include "TreeBufferHelper.h"// CTreeBufferHelper
#include "MainProgramParameters.h"// CMainProgramParameters
//#include "SharedFunctionsAndTypes.h"// CWeightLocation CWeightLocationCompare
//#include <valarray> // std::valarray<double>
#include <set> // std::set< CWeightLocation, CWeightLocationCompare > BASE
//#include <map>
#include <vector>// std::vector< BASE::iterator >
//#include "Cmove1.h"
//#include <stdlib.h>
#include "recombine.h"
#include "morton.h"
#include "wpUtilities.h"
//#include "Utils.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace CWeightLocation2D
{
#define dim__ 2
	typedef wpUtilities::Point_POD<dim__> Point_POD;
	typedef morton::block<Point_POD>::BLOCK BLOCK;
	typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
	typedef morton::block<Point_POD>::compare COMPARE;
	typedef wpUtilities::Point<dim__> Point;

	struct CWeightLocation
	{
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version)
		{
			ar& probability;
			ar& displacement;
		}

		mutable double probability;
		//double displacement;
		//wpUtilities::Point<dim__> displacement;
		Point displacement;
		//wpUtilities2D::Point displacement;
	};

	struct CWeightLocationCompare : public std::binary_function<bool, CWeightLocation, CWeightLocation>
	{	// functor for operator+
		bool operator()(const CWeightLocation& _Left, const CWeightLocation& _Right) const
		{
			/*wpUtilities::Point_POD<dim__> _Left_POD;
			wpUtilities::Point_POD<dim__> _Right_POD;*/
			Point_POD _Left_POD;
			Point_POD _Right_POD;
			//wpUtilities2D::Point_POD _Left_POD;
			//wpUtilities2D::Point_POD _Right_POD;

			std::copy(std::begin(_Left.displacement), std::end(_Left.displacement), std::begin(_Left_POD));
			std::copy(std::begin(_Right.displacement), std::end(_Right.displacement), std::begin(_Right_POD));

			//wpUtilities::BLOCK<dim__> _Left_BLOCK(_Left_POD);
			//wpUtilities::BLOCK<dim__> _Right_BLOCK(_Right_POD);
			//wpUtilities2D::BLOCK _Left_BLOCK(_Left_POD);
			//wpUtilities2D::BLOCK _Right_BLOCK(_Right_POD);
			BLOCK _Left_BLOCK(_Left_POD);
			BLOCK _Right_BLOCK(_Right_POD);

			return (_Left_BLOCK < _Right_BLOCK)
				? true
				: (_Left_BLOCK == _Right_BLOCK && _Left.probability < _Right.probability)
				? true
				: false;
		}
	};

	struct CWeightLocationMortonCompare : public std::binary_function<bool, CWeightLocation, CWeightLocation>
	{	
		CWeightLocationMortonCompare(size_t _uiBitMask) : uiBitMask(_uiBitMask) {}

		// Default 52 bit mask treats whole dyadic cube as one.
		// Here dyadic cube is fixed by the sign and exponent bits. We want the comparison to only consider fraction bits.
		CWeightLocationMortonCompare() : uiBitMask(52) {}  

		size_t uiBitMask;

		bool operator()(const CWeightLocation& _Left, const CWeightLocation& _Right) const
		{
			wpUtilities::Point_POD<dim__> _Left_POD;
			wpUtilities::Point_POD<dim__> _Right_POD;
			//wpUtilities2D::Point_POD _Left_POD;
			//wpUtilities2D::Point_POD _Right_POD;

			std::copy(std::begin(_Left.displacement), std::end(_Left.displacement), std::begin(_Left_POD));
			std::copy(std::begin(_Right.displacement), std::end(_Right.displacement), std::begin(_Right_POD));

			//wpUtilities::BLOCK<dim__> _Left_BLOCK(_Left_POD);
			//wpUtilities::BLOCK<dim__> _Right_BLOCK(_Right_POD);
			//wpUtilities2D::BLOCK _Left_BLOCK(_Left_POD);
			//wpUtilities2D::BLOCK _Right_BLOCK(_Right_POD);
			BLOCK _Left_BLOCK(_Left_POD);
			BLOCK _Right_BLOCK(_Right_POD);

			//wpUtilities::COMPARE<dim__> less(uiBitMask);
			COMPARE less(uiBitMask);
			//wpUtilities2D::COMPARE less(uiBitMask);

			return (less(_Left_BLOCK, _Right_BLOCK))
				? true
				//: ((!(less(_Right_BLOCK, _Left_BLOCK))) && _Left.probability < _Right.probability)
				//? true
				: false;
		}
	};

	typedef std::vector<CWeightLocation> CVectorWeightLocation;
	typedef std::set<CWeightLocation, CWeightLocationCompare> CSetWeightLocation;
	
	class CCloudWeightLocation :
		public CSetWeightLocation //std::set<CWeightLocation, CWeightLocationCompare>
	{
		//typedef std::set<CWeightLocation, CWeightLocationCompare> BASE;
		typedef CSetWeightLocation BASE;
	public:
		CCloudWeightLocation(void);
		~CCloudWeightLocation(void);

		// reduces the pointset based on reducing the number of points in each interval of given width to the accuracy of the cubature
		//CCloudWeightLocation& ReducePointSet(const double& dAccuracySpread, const unsigned int& stCubatureDegree,
		//	CMainProgramParameters& ipIntParams);

		CCloudWeightLocation& ReducePointSet(const size_t& uiBitMask, const unsigned int& stCubatureDegree); 
		//,CMainProgramParameters& ipIntParams);

		// breaks the cloud up into components // depends on the model
		void PartitionCloud(const size_t& uiBitMask, std::vector<std::vector<BASE::iterator> >&
			vvplwPartitionedCloud);

		void DoAbstractPruning(std::vector<BASE::iterator>& vplwToBeProcessed, const unsigned int& stCubatureDegree);
			//CMainProgramParameters& ipIntParams);

	};

	/*template <typename _Compare>
	class CCloudWeightLocationMorton : public std::set<CWeightLocation, _Compare>
	{
		typedef std::set<CWeightLocation, _Compare> BASE;
	public:
		CCloudWeightLocationMorton(void);
		~CCloudWeightLocationMorton(void);

		CCloudWeightLocationMorton& ReducePointSet(const size_t& uiBitMask, const unsigned int& stCubatureDegree,
			CMainProgramParameters& ipIntParams);

		void PartitionCloud(const size_t& uiBitMask, std::vector<std::vector<BASE::iterator> >&
			vvplwPartitionedCloud);

		void DoAbstractPruning(std::vector<BASE::iterator>& vplwToBeProcessed, const unsigned int& stCubatureDegree,
			CMainProgramParameters& ipIntParams);

	};*/

}
#endif // CloudWeightLocation_h__
