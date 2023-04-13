/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#include "stdafx.h"
#include "cloudweightlocation.h"
//#include <iostream>
#include <map>
#include <crtdbg.h>
//#include "Cmove1.h"
#include <cmath>


CCloudWeightLocation::CCloudWeightLocation(void) {}

CCloudWeightLocation::~CCloudWeightLocation(void) {}

// reduces the pointset based on reducing the number of points in each interval of given width to the accuracy of the cubature
CCloudWeightLocation& CCloudWeightLocation::ReducePointSet(const double& dAccuracySpread,
														   const unsigned int& iNoDimensionsToCubature,
														   CMainProgramParameters& ipIntParams)
{
	// CDebugTools&dbgLog = ipIntParams.dbgLogFile;
	//SHOW(dAccuracySpread);
	std::vector<std::vector<BASE::iterator> > vvplwToBeProcessed;
	PartitionCloud(dAccuracySpread, vvplwToBeProcessed);
	//SHOW(vvplwToBeProcessed.size());
	for (size_t uiPatchNo = 0; uiPatchNo < vvplwToBeProcessed.size(); ++uiPatchNo)
	{
		std::vector<BASE::iterator>& vplwToBeProcessed(vvplwToBeProcessed[uiPatchNo]);
		// vplwToBeProcessed now references the iterator_pointers to the points and weights 
		// that should be reduced in this cluster
		//SHOW(vplwToBeProcessed.size());
		DoAbstractPruning(vplwToBeProcessed, iNoDimensionsToCubature, ipIntParams);
		//SHOW(vplwToBeProcessed.size());
	}

	return *this;
}

// maps raw point location data to the vector of vectors locations and moments and a vector of weights for use with the svd solver 
void PrepareDataForCMove(const CSetWeightLocation& swlData, const int iNoPoints,
						 std::vector<double>& vdCMoveWeightsBuffer, std::vector<double>& vdCMovePointsBuffer,
						 const int& iNoCoords, const CWeightLocation& wlCOG)
{
	_ASSERT((iNoPoints > 0) && !(iNoPoints > (int)swlData.size()));

	vdCMoveWeightsBuffer.erase(vdCMoveWeightsBuffer.begin(), vdCMoveWeightsBuffer.end());
	vdCMovePointsBuffer.erase(vdCMovePointsBuffer.begin(), vdCMovePointsBuffer.end());
	_ASSERT((int)vdCMovePointsBuffer.capacity() >= iNoCoords * iNoPoints);
	_ASSERT((int)vdCMoveWeightsBuffer.capacity() >= iNoPoints);

	CSetWeightLocation::const_iterator itRawPoints = swlData.begin();
	for (int i = 0; i < iNoPoints; i++)
	{
		vdCMoveWeightsBuffer.push_back(itRawPoints->probability);
		double point = itRawPoints->displacement;
		itRawPoints++;
		for (int j = 0; j < iNoCoords; j++)
			vdCMovePointsBuffer.push_back(pow(point - ((j > 1) ? wlCOG.displacement : 0.), j));
	}
}

//void CCloudWeightLocation::DoAbstractPruning
//( 
// std::vector < BASE::iterator > &vplwToBeProcessed, 
// const unsigned int &iNoDimensionsToCubature,
// CMainProgramParameters &ipIntParams
// )
//{
//	CDebugTools&dbgLog = ipIntParams.dbgLogFile;
//	if (ipIntParams.bVerbose) std::cout << ":" << ((int) vplwToBeProcessed.size());
//
//	unsigned int ICountCalls(0);
//
//	if (iNoDimensionsToCubature + 1 <= vplwToBeProcessed.size())
//	{
//		// set parameters for the data
//		CTreeBufferHelper 
//			bhBufInfo(iNoDimensionsToCubature + 1, vplwToBeProcessed.size());
//		//                 iNoTrees,                  iNoInitialLeaves
//
//		// Create Buffers
//		std::vector <double> vdWeightsBuffer
//			( 
//			bhBufInfo.end(), 
//			double(-1) 
//			);
//		std::vector < std::valarray <double> > vdPointsBuffer
//			( 
//			bhBufInfo.end(), 
//			std::valarray <double>( double(0) , size_t( bhBufInfo.iNoTrees - 1 ) ) 
//			);
//		
//		// for the answers
//		std::vector < double > weights;
//		std::map <size_t , size_t> miTreePosition;
//
//		ICountCalls = IdentifyRebalancing(bhBufInfo, vdWeightsBuffer, vplwToBeProcessed, vdPointsBuffer, ICountCalls, miTreePosition, weights);
//
//
//		UpdateCloud(bhBufInfo, miTreePosition, vplwToBeProcessed, weights);
//
//	}
//
//	if (ipIntParams.bVerbose) 
//		std::cout<< "/" << ICountCalls;
//}





void CCloudWeightLocation::PartitionCloud(const double& dAccuracySpread, std::vector<std::vector<BASE::iterator> >&
										  vvplwPartitionedCloud)
{
	BASE::iterator itBelow(begin());// first weight and location
	//	BASE::iterator itAbove(); // upper bound weight and location

	vvplwPartitionedCloud.clear();

	while (itBelow != end())
	{
		// collect the points into clusters
		// with diameter dAccuracySpread
		CWeightLocation
			wlTopOfSpan(*itBelow);
		wlTopOfSpan.displacement += abs(dAccuracySpread);
		wlTopOfSpan.probability = 0;
		BASE::iterator
			itAbove = lower_bound(wlTopOfSpan);

		// populate a vector of iterators to a cluster of nearby points in the base class
		vvplwPartitionedCloud.resize(vvplwPartitionedCloud.size() + 1);
		for (BASE::iterator itLoop = itBelow; itLoop != itAbove; ++itLoop)
			vvplwPartitionedCloud.rbegin()->push_back(itLoop);

		itBelow = itAbove;
	}
}


//************************************
// Method:    ForestOfWeightedVectorsFromWeightedLeafVectors
// FullName:  CCloudWeightLocation::ForestOfWeightedVectorsFromWeightedLeafVectors
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: const CTreeBufferHelper &bhBufInfo
// Parameter: std::vector<double> &vdWeightsBuffer
// Parameter: std::vector< std::valarray<double> > &vdPointsBuffer
//************************************
void CCloudWeightLocation::ForestOfWeightedVectorsFromWeightedLeafVectors(const CTreeBufferHelper& bhBufInfo,
																		  std::vector<double>& vdWeightsBuffer,
																		  std::vector<std::valarray<double> >&
																		  vdPointsBuffer)
{
	for (size_t iIndex = bhBufInfo.iInitialNoLeaves; iIndex < bhBufInfo.end(); iIndex++)
	{
		size_t uiLeftParent = bhBufInfo.left(iIndex);
		size_t uiRightParent = bhBufInfo.right(iIndex);
		double left = vdWeightsBuffer[uiLeftParent];
		double right = vdWeightsBuffer[uiRightParent];
		double sum = left + right;
		double ratio = left / right;
		double rightp = right / sum;
		vdWeightsBuffer[iIndex] = sum;
		std::valarray<double>& dPointsBuffer = vdPointsBuffer[iIndex];
		dPointsBuffer = vdPointsBuffer[uiLeftParent];
		dPointsBuffer *= ratio;
		dPointsBuffer += vdPointsBuffer[uiRightParent];
		dPointsBuffer *= rightp;
	}
}

void CCloudWeightLocation::ExpandRawLocationWeightsToConditionedVectorsOfDoubles(const CTreeBufferHelper& bhBufInfo,
																				 std::vector<double>& vdWeightsBuffer,
																				 const std::vector<BASE::iterator>&
																				 vplwToBeProcessed,
																				 std::vector<std::valarray<double> >&
																				 vdPointsBuffer)
{
	// set conditioning data
	CWeightLocation wlCOG;
	double dApproxStandardDeviation;
	COGVar(wlCOG, bhBufInfo, vdWeightsBuffer, vplwToBeProcessed, dApproxStandardDeviation);

	// then locations expanded into vectors of moments (Points) normalised by the cog
	for (size_t iIndex = 0; bhBufInfo.isleaf(iIndex); iIndex ++)
	{
		vdWeightsBuffer[iIndex] = vplwToBeProcessed[iIndex]->probability;

		const size_t iDepthOfVector(bhBufInfo.iNoTrees - 1);
		const BASE::iterator& plwToBeProcessed(vplwToBeProcessed[iIndex]);
		std::valarray<double>& dPointsBuffer(vdPointsBuffer[iIndex]);
		LocationToVector(iDepthOfVector,
			dPointsBuffer,
			plwToBeProcessed,
			wlCOG,
			dApproxStandardDeviation);
	}
}

void CCloudWeightLocation::COGVar(CWeightLocation& wlCOG, const CTreeBufferHelper& bhBufInfo,
								  std::vector<double>& vdWeightsBuffer, const std::vector<BASE::iterator>&
								  vplwToBeProcessed, double& dApproxStandardDeviation)
{
	// weights first, computing COG on way
	wlCOG.displacement = 0;
	wlCOG.probability = 0;
	double dVariance(0);
	for (size_t iIndex = 0; bhBufInfo.isleaf(iIndex); iIndex ++)
	{
		wlCOG.probability += (vdWeightsBuffer[iIndex] = vplwToBeProcessed[iIndex]->probability);
		wlCOG.displacement += (vplwToBeProcessed[iIndex]->displacement * vplwToBeProcessed[iIndex]->probability);
		dVariance += (vplwToBeProcessed[iIndex]->displacement * vplwToBeProcessed[iIndex]->displacement *
					  vplwToBeProcessed[iIndex]->probability);
	}
	_ASSERT(bhBufInfo.isleaf(0) && wlCOG.probability != 0);
	wlCOG.displacement /= wlCOG.probability;
	dVariance /= wlCOG.probability;
	dVariance -= wlCOG.displacement * wlCOG.displacement;
	dApproxStandardDeviation = sqrt(dVariance);
}

void CCloudWeightLocation::LocationToVector(const size_t iDepthOfVector, std::valarray<double>& dPointsBuffer,
											const std::set<CWeightLocation, CWeightLocationCompare>::iterator&
											plwToBeProcessed, CWeightLocation& wlCOG, double& dApproxStandardDeviation)
{
	_ASSERT(iDepthOfVector <= dPointsBuffer.size());
	for (size_t j = 0; j < iDepthOfVector; j++)
		if (j < size_t(2))
		{
			dPointsBuffer[j] = (pow(plwToBeProcessed->displacement, (int)j));
		}
		else
		{
			dPointsBuffer[j] = (pow((plwToBeProcessed->displacement - wlCOG.displacement) /
								dApproxStandardDeviation, (int)j));
		}
}
//
//unsigned int CCloudWeightLocation::IdentifyLocationsRemainingAndTheirNewWeights( CTreeBufferHelper &bhBufInfo, std::map<size_t , size_t > &miTreePosition, std::vector<double> &vdWeightsBuffer, std::vector< std::valarray<double> > &vdPointsBuffer, std::vector < double >&weights, unsigned int &ICountCalls )
//{
//	weights.clear();
//	weights.resize(bhBufInfo.iNoTrees); 
//	std::vector < double > points (bhBufInfo.iNoTrees * (bhBufInfo.iNoTrees - 1)) ; 
//	std::vector < double > MassCog (bhBufInfo.iNoTrees);
//	std::vector < size_t > currentroots (bhBufInfo.iNoTrees);
//	std::vector <int> maxset ;
//	for (size_t iTreeIndexInFixedBuffer = 0; iTreeIndexInFixedBuffer < bhBufInfo.iNoTrees ; iTreeIndexInFixedBuffer ++ )
//	{
//		size_t currentroot = currentroots[iTreeIndexInFixedBuffer] = iTreeIndexInFixedBuffer + bhBufInfo.end() -  bhBufInfo.iNoTrees;
//		miTreePosition[currentroot] = iTreeIndexInFixedBuffer;
//		weights[iTreeIndexInFixedBuffer] = vdWeightsBuffer[currentroot];
//		for (size_t iM = 0; iM < bhBufInfo.iNoTrees -1 ; iM ++ )
//			points[iTreeIndexInFixedBuffer * ( bhBufInfo.iNoTrees - 1 ) + iM ] = (vdPointsBuffer[currentroot])[iM];
//	}
//
//	CLinearAlgebraReductionTool moLinearAlgebraReductionTool; 
//	_ASSERT((size_t)( (int) bhBufInfo.iNoTrees) == bhBufInfo.iNoTrees);
//	moLinearAlgebraReductionTool.INoCoords( (int) bhBufInfo.iNoTrees - 1) ;
//	moLinearAlgebraReductionTool.INoPoints( (int) bhBufInfo.iNoTrees) ;
//	while (!bhBufInfo.isleaf(miTreePosition.rbegin()->first))
//	{
//		moLinearAlgebraReductionTool.MoveMass(
//			weights , 
//			points , 
//			MassCog, 
//			maxset
//			);
//		while (maxset.size() &&! bhBufInfo.isleaf(miTreePosition.rbegin()->first))
//		{
//			size_t togoposition(maxset.back());
//			maxset.pop_back();
//			miTreePosition.erase(currentroots[togoposition]);
//			currentroots[togoposition] = -1;
//			size_t tosplit(miTreePosition.rbegin()->first);
//			if (! bhBufInfo.isleaf(tosplit))
//			{
//				size_t tosplitposition=miTreePosition[tosplit];
//				miTreePosition.erase(tosplit);
//				currentroots[tosplitposition] = -1;
//
//				weights[togoposition] = weights[tosplitposition] * vdWeightsBuffer[bhBufInfo.left(tosplit)]/vdWeightsBuffer[tosplit];
//				weights[tosplitposition] *= vdWeightsBuffer[bhBufInfo.right(tosplit)]/vdWeightsBuffer[tosplit];
//				for (size_t iM = 0; iM < bhBufInfo.iNoTrees -1 ; iM ++ )
//				{
//					points[togoposition * ( bhBufInfo.iNoTrees - 1 ) + iM ] = (vdPointsBuffer[bhBufInfo.left(tosplit)])[iM];
//					points[tosplitposition * ( bhBufInfo.iNoTrees - 1 ) + iM ] = (vdPointsBuffer[bhBufInfo.right(tosplit)])[iM];
//				}
//				currentroots[togoposition]=bhBufInfo.left(tosplit);
//				miTreePosition[bhBufInfo.left(tosplit)] = togoposition;
//				currentroots[tosplitposition]=bhBufInfo.right(tosplit);
//				miTreePosition[bhBufInfo.right(tosplit)] = tosplitposition;
//			}
//		}
//		ICountCalls = moLinearAlgebraReductionTool.INoCallsLinAlg();
//	}
//	while (maxset.size())
//	{
//		size_t togoposition(maxset.back());
//		maxset.pop_back();
//		miTreePosition.erase(currentroots[togoposition]);
//		currentroots[togoposition] = -1;
//	}	return ICountCalls;
//}

//unsigned int CCloudWeightLocation::IdentifyRebalancing( 
//	CTreeBufferHelper &bhBufInfo, 
//	std::vector<double> &vdWeightsBuffer, 
//	std::vector< BASE::iterator > &vplwToBeProcessed, 
//	std::vector< std::valarray<double> > &vdPointsBuffer, 
//	unsigned int &ICountCalls, 
//	std::map<size_t , size_t> &miTreePosition, 
//	std::vector< double > &weights 
//	)
//{
//	// insert the raw data as iNoLeaves leaves
//	ExpandRawLocationWeightsToConditionedVectorsOfDoubles(
//		bhBufInfo, 
//		vdWeightsBuffer, 
//		vplwToBeProcessed, 
//		vdPointsBuffer
//		);
//
//	// now fill the forest using the leafs, leaving exactly iNoTrees roots
//	ForestOfWeightedVectorsFromWeightedLeafVectors(
//		bhBufInfo, 
//		vdWeightsBuffer, 
//		vdPointsBuffer
//		);
//
//	ICountCalls = IdentifyLocationsRemainingAndTheirNewWeights(
//		bhBufInfo, 
//		miTreePosition, 
//		vdWeightsBuffer, 
//		vdPointsBuffer, 
//		weights, 
//		ICountCalls
//		);	
//	return ICountCalls;
//}
