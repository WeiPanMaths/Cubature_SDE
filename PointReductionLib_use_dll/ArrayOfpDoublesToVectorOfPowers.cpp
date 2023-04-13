#include "stdafx.h"
#include "topowersdata.h"
#include "arrayofpdoublestovectorofpowers.h"
#include <cmath>


void ArrayOfpDoublesToVectorOfPowers(void* pIn, double* pOut, void* vpCBufferHelper)
{
	CBufferHelper* pCBufferHelper = (CBufferHelper*)vpCBufferHelper;
	const size_t no_of_locations = pCBufferHelper->NoPointsToBeprocessed;
	const size_t depth_of_vector = pCBufferHelper->SmallestReducibleSetSize - 1;

	//pIn is a null pointer to an array of null pointers that really point to doubles
	void** pVoidIn = (void**)pIn; // a pointer to the first element of an array of null pointers
	for (size_t i = 0; i < no_of_locations; ++i, ++pVoidIn)
		for (size_t j = 0; j < depth_of_vector; ++j, ++pOut)
			*pOut = pow(*((double*)(*pVoidIn)), (int)j);
}

void ArrayOfpDoublesToVectorOfPowers1(void* pIn, double* pOut, void* vpCConditionedBufferHelper)
{
	CConditionedBufferHelper* pCConditionedBufferHelper = (CConditionedBufferHelper*)vpCConditionedBufferHelper;
	const size_t no_of_locations = pCConditionedBufferHelper->NoPointsToBeprocessed;
	const size_t depth_of_vector = pCConditionedBufferHelper->SmallestReducibleSetSize - 1;
	
	CConditioning* pConditioning = (CConditioning*)pCConditionedBufferHelper->pvCConditioning;
	const double dDisplacement = pConditioning->dMean;
	const double dStdDev = pConditioning->dStdDev;

	void** pVoidIn = (void**)pIn;
	for (size_t i = 0; i < no_of_locations; ++i, ++pVoidIn)
		for (size_t j = 0; j < depth_of_vector; ++j, ++pOut)
			if (j < size_t(2) || (dDisplacement == double(0) && dStdDev == double(1)))
			{
				*pOut = pow(*((double*)(*pVoidIn)), (int)j);
			}
			else
			{
				*pOut = pow((*((double*)(*pVoidIn)) - dDisplacement) / dStdDev, (int)j);
			}
}
