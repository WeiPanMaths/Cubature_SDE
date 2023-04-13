#include "stdafx.h"
#include "recombine.h"
#include "BufferConstructor.h"
#include "TreeBufferHelper.h"
#include "Cmove1.h" // CLinearAlgebraReductionTool moLinearAlgebraReductionTool
#include <map>
#include <valarray>
#include "..\test\show.h"


void ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
													std::vector<double>& vdWeightsBuffer,
													std::vector<std::valarray<double> >& vdPointsBuffer)
{
	//SHOW( bhBufInfo.iInitialNoLeaves );
	//SHOW( bhBufInfo.end());

	for (size_t iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++)
	{
		size_t uiLeftParent = bhBufInfo.left(iIndex);
		size_t uiRightParent = bhBufInfo.right(iIndex);
		double left = vdWeightsBuffer[uiLeftParent];
		double right = vdWeightsBuffer[uiRightParent];
		double sum = left + right;		
		vdWeightsBuffer[iIndex] = sum;
		std::valarray<double>& dPointsBuffer = vdPointsBuffer[iIndex];
		if (left <= right)
			dPointsBuffer = vdPointsBuffer[uiLeftParent] * (left / sum) + vdPointsBuffer[uiRightParent] * (1 - (left /
																												sum));
		else
			dPointsBuffer = vdPointsBuffer[uiLeftParent] * (1 - (right / sum)) +
				vdPointsBuffer[uiRightParent] * (right / sum);
		//SHOW(right);
		//SHOW(left);
		//SHOW(dPointsBuffer);
	}
}
unsigned int IdentifyLocationsRemainingAndTheirNewWeights(CTreeBufferHelper& bhBufInfo, std::map<size_t, size_t>&
														  miTreePosition, std::vector<double>& vdWeightsBuffer,
														  std::vector<std::valarray<double> >& vdPointsBuffer,
														  std::vector<double>& weights, unsigned int& ICountCalls)
{
	weights.clear();
	weights.resize(bhBufInfo.iNoTrees);
	std::vector<double> points (bhBufInfo.iNoTrees * (bhBufInfo.iNoTrees - 1));
	std::vector<double> MassCog (bhBufInfo.iNoTrees);
	std::vector<ptrdiff_t> currentroots (bhBufInfo.iNoTrees);
	std::vector<int> maxset;
	bool SomeLinearAlgebraToDo = (bhBufInfo.end() >= bhBufInfo.iNoTrees);
	_ASSERT(SomeLinearAlgebraToDo);
	for (size_t iTreeIndexInFixedBuffer = 0; iTreeIndexInFixedBuffer < bhBufInfo.iNoTrees; iTreeIndexInFixedBuffer ++)
	{
		ptrdiff_t currentroot = currentroots[iTreeIndexInFixedBuffer] = iTreeIndexInFixedBuffer + bhBufInfo.end() -
			bhBufInfo.iNoTrees;
		miTreePosition[(size_t)currentroot] = iTreeIndexInFixedBuffer;
		weights[iTreeIndexInFixedBuffer] = vdWeightsBuffer[currentroot];
		for (size_t iM = 0; iM < bhBufInfo.iNoTrees - 1; iM ++)
			points[iTreeIndexInFixedBuffer * (bhBufInfo.iNoTrees - 1) + iM ] = (vdPointsBuffer[currentroot])[iM];
	}

	CLinearAlgebraReductionTool moLinearAlgebraReductionTool;
	_ASSERT((size_t)((int)bhBufInfo.iNoTrees) == bhBufInfo.iNoTrees);
	moLinearAlgebraReductionTool.INoCoords((int)bhBufInfo.iNoTrees - 1);
	moLinearAlgebraReductionTool.INoPoints((int)bhBufInfo.iNoTrees);
	while (SomeLinearAlgebraToDo)
	{
		moLinearAlgebraReductionTool.MoveMass(weights,
			points,
			MassCog,
			maxset);
		while (maxset.size())
		{
			size_t togoposition(maxset.back());
			maxset.pop_back();
			miTreePosition.erase(currentroots[togoposition]);
			currentroots[togoposition] = -1;
			// if there is at least one non-trivial tree split the last (and so deepest) one to fill vacant slot
			size_t tosplit(miTreePosition.rbegin()->first);
			if (!bhBufInfo.isleaf(tosplit))
			{
				size_t tosplitposition = miTreePosition[tosplit];
				miTreePosition.erase(tosplit);
				currentroots[tosplitposition] = -1;

				weights[togoposition] = weights[tosplitposition] * vdWeightsBuffer[bhBufInfo.left(tosplit)] /
					vdWeightsBuffer[tosplit];
				weights[tosplitposition] *= vdWeightsBuffer[bhBufInfo.right(tosplit)] / vdWeightsBuffer[tosplit];
				for (size_t iM = 0; iM < bhBufInfo.iNoTrees - 1; iM ++)
				{
					points[togoposition * (bhBufInfo.iNoTrees - 1) + iM ] = (vdPointsBuffer[bhBufInfo.
																			 left(tosplit)])[iM];
					points[tosplitposition * (bhBufInfo.iNoTrees - 1) + iM ] = (vdPointsBuffer[bhBufInfo.
																				right(tosplit)])[iM];
				}
				currentroots[togoposition] = bhBufInfo.left(tosplit);
				miTreePosition[bhBufInfo.left(tosplit)] = togoposition;
				currentroots[tosplitposition] = bhBufInfo.right(tosplit);
				miTreePosition[bhBufInfo.right(tosplit)] = tosplitposition;
			}
			else
				SomeLinearAlgebraToDo = false;
		}
		ICountCalls = moLinearAlgebraReductionTool.INoCallsLinAlg();
	}
	return ICountCalls;
}

void RECOMBINE_API Recombine(void* pData)
{
	// unpack the void pointer
	sRecombineInterface& data = *(sRecombineInterface*)pData;

	sCloud& InCloud = *data.pInCloud;
	size_t& NoActiveWeightsLocations = InCloud.NoActiveWeightsLocations;
	double*& WeightBuf = InCloud.WeightBuf;
	void*& LocationBuf = InCloud.LocationBuf;
	void*& InEnd = InCloud.end; // can point to additional data needed by the expander function

	sRCloudInfo& OutCloudInfo = *data.pOutCloudInfo;
	size_t& No_KeptLocations = OutCloudInfo.No_KeptLocations; // number actually returned
	double*& NewWeightBuf = OutCloudInfo.NewWeightBuf;  // a buffer containing the weights of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1
	size_t*& KeptLocations = OutCloudInfo.KeptLocations; // a buffer containing the offsets of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1
	void*& OutEnd = OutCloudInfo.end;

	// iNoDimensionsToCubature is the number of dimensions in the vectors (for 1st moment only this would be 2) 
	size_t& iNoDimensionsToCubature = data.degree;
	size_t SmallestReducibleSetSize(iNoDimensionsToCubature + 1);
	if (SmallestReducibleSetSize > NoActiveWeightsLocations)
	{
		// reduction already achieved
		double* pwIn = WeightBuf;
		double* pwOut = NewWeightBuf;
		size_t* pl = KeptLocations;
		No_KeptLocations = NoActiveWeightsLocations;
		for (size_t iIndex = 0; iIndex < NoActiveWeightsLocations; iIndex++)
		{
			*(pwOut++) = *(pwIn++);
			*(pl++) = iIndex;
		}
	}
	else
	{
		CTreeBufferHelper bhBufInfo(SmallestReducibleSetSize, NoActiveWeightsLocations);
		expander ArrayOfpPointsToVectorDoubles = &(*data.fn);
		std::vector<double> vdWeightsBuffer(bhBufInfo.end());
		std::copy(WeightBuf, WeightBuf + bhBufInfo.iInitialNoLeaves, vdWeightsBuffer.begin());
		std::valarray<double> vdArrayPointsBuffer (bhBufInfo.iInitialNoLeaves * (bhBufInfo.iNoTrees - 1));
		CConditionedBufferHelper arg3;
		arg3.NoPointsToBeprocessed = bhBufInfo.iInitialNoLeaves;
		arg3.SmallestReducibleSetSize = bhBufInfo.iNoTrees;
		arg3.pvCConditioning = InEnd;
		//SHOW(LocationBuf);
		ArrayOfpPointsToVectorDoubles(LocationBuf, &vdArrayPointsBuffer[0], &arg3);
		//SHOW(LocationBuf);
		//SHOW(vdArrayPointsBuffer);
		size_t vecdim(bhBufInfo.iNoTrees - 1); //The vector of powers starts with the zero'th power two left at the end for first moments
		// double zero(0);
		std::valarray<double> zero_vec(vecdim);
		zero_vec = 0;
		// map buffer to array of val arrays for compatibility reasons
		std::vector<std::valarray<double> > vdPointsBuffer(bhBufInfo.end());
		// populate the leaves
		for (size_t i = 0; i < bhBufInfo.iInitialNoLeaves; i++)
		{
			std::slice vaSlice ((bhBufInfo.iNoTrees - 1) * i , (bhBufInfo.iNoTrees - 1) , 1);
			vdPointsBuffer[i] = vdArrayPointsBuffer[ vaSlice ];
		}
		//SHOW(vdPointsBuffer);

		// now fill the forest using the leafs, leaving exactly iNoTrees roots
		ForestOfWeightedVectorsFromWeightedLeafVectors(bhBufInfo,
			vdWeightsBuffer,
			vdPointsBuffer);	

		unsigned int ICountCalls;
		std::map<size_t, size_t> miTreePosition;
		std::vector<double> weights;

		ICountCalls = IdentifyLocationsRemainingAndTheirNewWeights(bhBufInfo,
			miTreePosition,
			vdWeightsBuffer,
			vdPointsBuffer,
			weights,
			ICountCalls);	

		//sRCloudInfo & OutCloudInfo = * data.pOutCloudInfo;
		//size_t & No_KeptLocations = OutCloudInfo.No_KeptLocations; // number actually returned
		//double* & NewWeightBuf = OutCloudInfo.NewWeightBuf;  // a buffer containing the weights of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1
		//size_t* & KeptLocations = OutCloudInfo.KeptLocations; // a buffer containing the offsets of the kept Locations // capacity must be at least iNoDimensionsToCubature + 1
		//void * & OutEnd = OutCloudInfo.end;

		double* pw = NewWeightBuf;
		size_t* pl = KeptLocations;
		No_KeptLocations = 0;//not weights.size();
		for (size_t iIndex = 0; bhBufInfo.isleaf(iIndex); iIndex++)
		{
			if (miTreePosition.find(iIndex) != miTreePosition.end())
			{
				++No_KeptLocations;
				*(pw++) = weights[miTreePosition[iIndex]];
				*(pl++) = iIndex;
			}
		}
	}
	return;
}

//void RECOMBINE_API Recombine(void* recombine_interface)
//{
//	sRecombineInterface* arg = (sRecombineInterface*) recombine_interface;
//	if(arg->end != 0){
//		//Invalid Structure passed to this version of Recombine
//		return;
//	}
//
//}
