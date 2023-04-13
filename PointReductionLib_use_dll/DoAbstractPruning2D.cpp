#include "stdafx.h"
//#include "SharedFunctionsAndTypes.h"
#include "recombine.h"
//#include "arrayofpdoublestovectorofpowers.h"
#include "cloudweightlocation2D.h"
//#include "topowersdata.h" // CConditioning
#include "RdToPowers.h" // CMultiDimensionalBufferHelper
#include "RdToPowers2.h" // RdToPowers
#include "EvaluateAllMonomials.h" //EvaluateAllMonomials::F
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <crtdbg.h>

namespace CWeightLocation2D {

	void CCloudWeightLocation::DoAbstractPruning
	(std::vector<std::set<CWeightLocation, CWeightLocationCompare>::iterator>& vplwToBeProcessed,
		//const unsigned int& iNoDimensionsToCubature, 
		const unsigned int& stCubatureDegree) //,
		//CMainProgramParameters& ipIntParams)
	{
		// create the input buffer of void pointers capturing the locations of points
		std::vector<const void*> locationbuffer;
		locationbuffer.resize(0);
		locationbuffer.reserve(vplwToBeProcessed.size());

		// create the input buffer of double pointers capturing the weights of points
		std::vector<double> weightbuffer;
		weightbuffer.resize(0);
		weightbuffer.reserve(vplwToBeProcessed.size());

		// extract data from vplwToBeProcessed to populate the buffers
		for (std::vector < std::set< CWeightLocation, CWeightLocationCompare >::iterator >::const_iterator
			it_in = vplwToBeProcessed.begin();
			it_in != vplwToBeProcessed.end();
			++it_in)
		{
			locationbuffer.push_back( &((*(*it_in)).displacement[0])  );
			weightbuffer.push_back((*(*it_in)).probability);
		}

		// set up the input structure for conditioning
		/*CConditioning sConditioning;
		sConditioning.dMean = dMean;
		sConditioning.dStdDev = dStdDev;*/
		CMultiDimensionalBufferHelper sConditioning;
		sConditioning.D = stCubatureDegree;
		sConditioning.L = dim__;

		// set up the input structure for data reduction "in"
		sCloud in;

		// chain optional extension information used to condition the data
		in.end = &sConditioning;

		// place the sizes of the buffers and their locations into the structure "in"
		in.NoActiveWeightsLocations = vplwToBeProcessed.size();
		in.LocationBuf = &locationbuffer[0];
		in.WeightBuf = &weightbuffer[0];
		//SHOW(weightbuffer);

		// set up the output structure for data reduction "out"
		sRCloudInfo out;
		out.end = 0;

		// in the new example, iNoDimensionsToCubature is calculated
		size_t iNoDimensionsToCubature = EvaluateAllMonomials::F(dim__, stCubatureDegree);
		//SHOW(iNoDimensionsToCubature);
		// setup a buffer of size iNoDimensionsToCubature to store indexes to the kept points
		std::vector<size_t> KeptLocations(iNoDimensionsToCubature);

		// setup a buffer of size iNoDimensionsToCubature to store the weights of the kept points
		std::vector<double> NewWeights(iNoDimensionsToCubature);

		// set the locations of these buffers into the structure "out"
		out.KeptLocations = &KeptLocations[0];
		out.NewWeightBuf = &NewWeights[0];

		// and the max dimension of the buffers
		out.No_KeptLocations = iNoDimensionsToCubature;

		//setup the Recombine Interface data which will join the input and output
		sRecombineInterface data;
		data.end = 0;

		// bind in and out together in data
		data.pInCloud = &in;
		data.pOutCloudInfo = &out;

		// add the degree of the vectors used and the callback function that expands 
		// the array of pointers to points into a long buffer of vectors
		data.degree = (size_t)iNoDimensionsToCubature;
		//data.fn = &ArrayOfpDoublesToVectorOfPowers1;
		data.fn = &RdToPowers;

		// CALL THE DLL THAT DOES THE HEAVYLIFTING
		Recombine(&data);
		//SHOW(NewWeights);
		//SHOW(out.No_KeptLocations);
		//SHOW(KeptLocations);

		// create a lookup service to determine which pointers should remain
		std::map<size_t, size_t> miIn2Out;
		for (size_t i = 0; i < out.No_KeptLocations; i++)
			miIn2Out[out.KeptLocations[i]] = i;

		// update weights and prune points from cloud
		for (size_t j = 0; j < in.NoActiveWeightsLocations; j++)
		{
			if (miIn2Out.find(j) != miIn2Out.end())
			{
				// change weight of the point
				std::set<CWeightLocation, CWeightLocationCompare>::value_type temp = *vplwToBeProcessed[j];//->probability = NewWeights[miIn2Out[j]];
				erase(vplwToBeProcessed[j]);
				temp.probability = NewWeights[miIn2Out[j]];
				vplwToBeProcessed[j] = insert(temp).first;
			}
			else
			{
				// remove point from base class
				erase(vplwToBeProcessed[j]);
				vplwToBeProcessed[j] = end();
			}
		}
	}
}

