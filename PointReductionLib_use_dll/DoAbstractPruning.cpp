#include "stdafx.h"
//#include "SharedFunctionsAndTypes.h"
#include "recombine.h"
#include "arrayofpdoublestovectorofpowers.h"
#include "cloudweightlocation.h"
#include "topowersdata.h"
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <crtdbg.h>

//************************************
// Method:    DoAbstractPruning
// FullName:  CCloudWeightLocation::DoAbstractPruning
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: std::vector < std::set< CWeightLocation, CWeightLocationCompare >::iterator > & vplwToBeProcessed will be modified 
// the targets of the iterators will be reweighted and or erased; in the latter case the iterator  will be set to point to this->end()
// Parameter: const unsigned int & iNoDimensionsToCubature = the depth of the vectors (including the one) (add 1 to get the number of trees)
// Parameter: CMainProgramParameters & ipIntParams
//************************************
void CCloudWeightLocation::DoAbstractPruning
	(std::vector<std::set<CWeightLocation, CWeightLocationCompare>::iterator>& vplwToBeProcessed,
	 const unsigned int& iNoDimensionsToCubature, CMainProgramParameters& ipIntParams)
{
	// create the input buffer of void pointers capturing the locations of points
	std::vector<const void*> locationbuffer;
	locationbuffer.resize(0);
	locationbuffer.reserve(vplwToBeProcessed.size());

	// create the input buffer of double pointers capturing the weights of points
	std::vector<double> weightbuffer;
	weightbuffer.resize(0);
	weightbuffer.reserve(vplwToBeProcessed.size());

	//populate the location and weight buffer and optionally collect conditioning data (an approximate Mean and Standard Deviation)
	double dMean(0);
	double dStdDev(1);
	{
		// It might be worth conditioning the cloud to have mean 0 and variance 1 before computing the moments
		// changing the line below to #if 1 will toggle between doing nothing and this scenario to compare
#if 0
		// extract data from vplwToBeProcessed
		for (std::vector < std::set< CWeightLocation, CWeightLocationCompare >::iterator >::const_iterator
			it_in = vplwToBeProcessed.begin ( ) ;
			it_in != vplwToBeProcessed.end ( );
			++it_in)
		{
			locationbuffer.push_back( &((*(*it_in)).displacement ));
			weightbuffer.push_back( (*(*it_in)).probability );
		}
#else
		double dWeight(0);
		// get the mean
		for (std::vector<std::set<CWeightLocation, CWeightLocationCompare>::iterator>::const_iterator
			 it_in = vplwToBeProcessed.begin();
			 it_in != vplwToBeProcessed.end();
			 ++it_in)
		{
			const double& displacement = (*(*it_in)).displacement;
			const double& probability = (*(*it_in)).probability;

			dWeight += probability;
			dMean += probability * displacement;
		}
		dMean /= dWeight;	// COM

		double dVar(0);
		// extract data from vplwToBeProcessed
		for (std::vector<std::set<CWeightLocation, CWeightLocationCompare>::iterator>::const_iterator
			 it_in = vplwToBeProcessed.begin();
			 it_in != vplwToBeProcessed.end();
			 ++it_in)
		{
			const double& displacement = (*(*it_in)).displacement;
			const double probability = (*(*it_in)).probability;
			dVar += probability * (displacement - dMean) * (displacement - dMean);

			locationbuffer.push_back(&displacement);
			weightbuffer.push_back(probability);
		}
		
		dVar /= dWeight;
		_ASSERT(dVar > 0. || vplwToBeProcessed.size() <= 1);
		dStdDev = sqrt(dVar);
#endif
	}

	// set up the input structure for conditioning
	CConditioning sConditioning;
	sConditioning.dMean = dMean;
	sConditioning.dStdDev = dStdDev;

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
	data.fn = &ArrayOfpDoublesToVectorOfPowers1;

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


