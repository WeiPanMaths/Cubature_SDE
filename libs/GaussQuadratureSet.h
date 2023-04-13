/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once

#include <crtdbg.h>
#include <vector>
#include "gqgenhermite.h"

template <unsigned int iPrecision>
class CHighPrecisionGaussianQuadratures
{
public:
	typedef amp::ampf<iPrecision> CHPNumber;
	typedef std::pair<CHPNumber, CHPNumber> CWeightLocation;

	static inline CHPNumber Pi()
	{
		return amp::pi<iPrecision>();
	}
	static inline CHPNumber Maximum(const CHPNumber& arg1, const CHPNumber& arg2)
	{
		return amp::maximum<iPrecision>(arg1, arg2);
	}
	static inline CHPNumber Sqrt(const CHPNumber& arg)
	{
		return amp::sqrt<iPrecision>(arg);
	}
	static inline CHPNumber Exp(const CHPNumber& arg)
	{
		return amp::exp<iPrecision>(arg);
	}

private:
	unsigned int mNoPoints;	
	std::vector<CWeightLocation> mNormalQuadrature;
	typedef ap::template_1d_array<CHPNumber> CapArrayHPNumbers;

public:	
	CHighPrecisionGaussianQuadratures(const unsigned int n = 4)
		: mNormalQuadrature (size_t(n)),
		  mNoPoints (n)
	{
		CapArrayHPNumbers mLocations, mWeights;
		gqgenhermite::buildgausshermitequadrature<iPrecision>(n, mLocations, mWeights);

		_ASSERT(mWeights.gethighbound() == mLocations.gethighbound() &&
			mWeights.getlowbound() == mLocations.getlowbound());
		_ASSERT(mWeights.gethighbound() + 1 - mWeights.getlowbound() == n);

		for (
			int j = 0, i = mWeights.getlowbound(); i != mWeights.gethighbound() + 1; i++
			)
		{
			mWeights(i) /= Sqrt(Pi());//probability measure
			mLocations(i) *= Sqrt(2); // variance 1
			mNormalQuadrature[j++] = CWeightLocation(mLocations(i), mWeights(i));
		}	
	}

	const std::vector<CWeightLocation>& operator ()(void) const
	{
		return mNormalQuadrature;
	}

	unsigned int Degree(void) const
	{
		return mNoPoints;
	}
};


class CGaussQuadratureSet
{
public:
	CGaussQuadratureSet(void);
	~CGaussQuadratureSet(void);
};
