#include "stdafx.h"
#include "cloudweightlocation2D.h"
//#include <iostream>
#include <map>
#include <crtdbg.h>
//#include "Cmove1.h"
#include <cmath>
#include <algorithm>

namespace CWeightLocation2D {

	CCloudWeightLocation::CCloudWeightLocation(void) {}

	CCloudWeightLocation::~CCloudWeightLocation(void) {}

	// reduces the pointset based on reducing the number of points in each interval of given width to the accuracy of the cubature
	/*CCloudWeightLocation& CCloudWeightLocation::ReducePointSet(const double& dAccuracySpread,
		const unsigned int& stCubatureDegree,
		CMainProgramParameters& ipIntParams)
	{
		size_t uiBitMask = 52;
		return ReducePointSet(uiBitMask, stCubatureDegree, ipIntParams);
	}*/

	CCloudWeightLocation& CCloudWeightLocation::ReducePointSet(const size_t& uiBitMask,
		const unsigned int& stCubatureDegree)
		//CMainProgramParameters& ipIntParams)
	{
		// CDebugTools&dbgLog = ipIntParams.dbgLogFile;
		//SHOW(dAccuracySpread);
		//size_t uiBitMask = ....;  // find the right mask value given the spread

		std::vector<std::vector<BASE::iterator> > vvplwToBeProcessed;
		PartitionCloud(uiBitMask, vvplwToBeProcessed);
		//SHOW(vvplwToBeProcessed.size());
		for (size_t uiPatchNo = 0; uiPatchNo < vvplwToBeProcessed.size(); ++uiPatchNo)
		{
			std::vector<BASE::iterator>& vplwToBeProcessed(vvplwToBeProcessed[uiPatchNo]);
			// vplwToBeProcessed now references the iterator_pointers to the points and weights 
			// that should be reduced in this cluster
			//SHOW(vplwToBeProcessed.size());
			DoAbstractPruning(vplwToBeProcessed, stCubatureDegree); // , ipIntParams);
			//SHOW(vplwToBeProcessed.size());
		}

		return *this;
	}
#if 0  // deprecated	
	void print_bits(double _value)
	{
		union {
			double  dValue;
			char array[sizeof(double)];
		};

		dValue = _value;
		std::cout << dValue << std::endl;

		for (int i = 0; i < sizeof(double) * CHAR_BIT; ++i) {
			int relativeToByte = i % CHAR_BIT;
			bool isBitSet = (array[sizeof(double) - 1 - i / CHAR_BIT] &
				(1 << (CHAR_BIT - relativeToByte - 1))) == (1 << (CHAR_BIT - relativeToByte - 1));
			std::cout << (isBitSet ? "1" : "0");
		}

		std::cout << std::endl;
	}

	// from tekpool.wordpress.com/category/bit-count/
	int BitCount(unsigned int u)
	{
		unsigned int uCount;
		uCount = u - ((u >> 1) & 033333333333) - ((u >> 2) & 011111111111);
		return ((uCount + (uCount >> 3)) & 030707070707) % 63;
	}

	int First0Bit(int i)
	{
		i = ~i; //bitwise complement
		return BitCount((i & (-i)) - 1);
	}

	// least significant bit has position index 0
	size_t FirstZeroBit(double pt)
	{
		union
		{
			double _val;
			unsigned int _ival[2];
		};
		_val = pt;
		return (First0Bit(_ival[0]) < 32)
			? First0Bit(_ival[0])
			: 32 + First0Bit(_ival[1]);
		/*try
		{
			auto ans = (First0Bit(_ival[0]) < 32)
				? First0Bit(_ival[0])
				: 32 + First0Bit(_ival[1]);
			if (ans == 64)
				throw(1);
			else
				return ans;
		}
		catch (int err)
		{
			std::cout << "-nan \n";
			throw err;
		}*/
	}

	/*  this currently only works for [1, 2)^2  */
	//void next_key(const size_t uiBitMask, wpUtilities2D::Point& pt)
	void next_key(const size_t uiBitMask, Point& pt)
	{
		union {
			double dValues[dim__];
			unsigned long long iValues[dim__];
		};
		std::copy_n(pt.begin(), dim__, dValues);

		// save exponent bits
		unsigned long long exponent_f = (iValues[0] >> 52) << 52;
		unsigned long long exponent_s = (iValues[1] >> 52) << 52;

		// clear the exponents bits of the original for "morton number" addition
		iValues[0] = exponent_f ^ iValues[0];
		iValues[1] = exponent_s ^ iValues[1];

		iValues[0] = (iValues[0] >> uiBitMask);
		iValues[1] = (iValues[1] >> uiBitMask);

		// print_bits(dValues[0]);
		// print_bits(dValues[1]);

		auto zero_pos_f = FirstZeroBit(dValues[0]);
		auto zero_pos_s = FirstZeroBit(dValues[1]);
		// iValues[0] = iValues[0] << uiBitMask;
		// iValues[1] = iValues[1] << uiBitMask;
		// std::cout << zero_pos_f << ", " << zero_pos_s << std::endl;
		if (zero_pos_f < zero_pos_s)
		{
			// first zero is in the first component
			iValues[0] += (unsigned long long) 1 << zero_pos_f;
			iValues[0] = (iValues[0] >> zero_pos_f) << zero_pos_f;
			iValues[1] = (iValues[1] >> (zero_pos_f + 1)) << (zero_pos_f + 1);
		}
		else
		{
			iValues[1] += (unsigned long long) 1 << zero_pos_s;
			iValues[1] = (iValues[1] >> zero_pos_s) << zero_pos_s;
			iValues[0] = (iValues[0] >> zero_pos_s) << zero_pos_s;
		}
		iValues[0] = (iValues[0] << uiBitMask) + exponent_f;
		iValues[1] = (iValues[1] << uiBitMask) + exponent_s;
		std::copy_n(dValues, dim__, pt.begin());
	}
#endif
	void CCloudWeightLocation::PartitionCloud(const size_t& uiBitMask, std::vector<std::vector<BASE::iterator> >&
		vvplwPartitionedCloud)
	{
		BASE::iterator itBelow(begin());
		vvplwPartitionedCloud.clear();
		CWeightLocationMortonCompare less(uiBitMask);

		// copy to vector to provide random access iterators
		std::vector<CWeightLocation> vwlcloud;
		vwlcloud.assign(this->begin(), this->end());
		std::vector<CWeightLocation>::iterator vitBelow(vwlcloud.begin());

		static_assert(sizeof(unsigned long long) == sizeof(double), "unsigned long long size assumptions invalid");
		static_assert(sizeof(unsigned int) == 4, "unsigned int size assumptions invalid");

		while (itBelow != end())
		{
			// given the patch diameter, we need the zip/block equivalent of
			// *itBelow.displacement += abs( patch size )
#if 0
			CWeightLocation	wlTopOfSpan(*itBelow);
			wlTopOfSpan.probability = 0;

			next_key(uiBitMask, wlTopOfSpan.displacement);
			BASE::iterator	itAbove = upper_bound(wlTopOfSpan);
#endif
			std::vector<CWeightLocation>::iterator vitAbove = std::upper_bound(vitBelow, vwlcloud.end(), *itBelow, less);
			BASE::iterator	itAbove;
			if (vitAbove == vwlcloud.end())
				itAbove = end();
			else
				itAbove = find(*vitAbove);

			// populate a vector of iterators to a cluster of nearby points in the base class
			vvplwPartitionedCloud.resize(vvplwPartitionedCloud.size() + 1);
			for (BASE::iterator itLoop = itBelow; itLoop != itAbove; ++itLoop)
				vvplwPartitionedCloud.rbegin()->push_back(itLoop);

			itBelow = itAbove;
			vitBelow = vitAbove;
		} // end while
	}

}