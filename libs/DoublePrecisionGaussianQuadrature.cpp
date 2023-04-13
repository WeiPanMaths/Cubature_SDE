/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#include ".\doubleprecisiongaussianquadrature.h"
#include <map>
#include "CloudWeightLocation.h"

#include <sstream>
#include <fstream>
#include <iostream>

#if 0
#include "GaussQuadratureSet.h"
#include <ios>

//#include "SharedFunctionsAndTypes.h"
typedef CHighPrecisionGaussianQuadratures<300> CFixedPrecisionGaussianQuadrature;

void Normal01QuadratureDoublePrecision(CVectorWeightLocation& vlwNormalApproxMean0Var1, const unsigned int iNoPoints)
{
	// create the high precision cubature (and remember it)
	static std::map<unsigned int, CFixedPrecisionGaussianQuadrature> vhpFixedPrecisionEmpiricalDistributions;
	// populate it as required
	if (vhpFixedPrecisionEmpiricalDistributions.find(iNoPoints) == vhpFixedPrecisionEmpiricalDistributions.end())
		vhpFixedPrecisionEmpiricalDistributions[iNoPoints] = CFixedPrecisionGaussianQuadrature(iNoPoints);		// here it was iNoPoints -1;
	
	// create the "double" quadrature
	vlwNormalApproxMean0Var1.clear();
	const CFixedPrecisionGaussianQuadrature& a = vhpFixedPrecisionEmpiricalDistributions[iNoPoints];	
	for (std::vector<CFixedPrecisionGaussianQuadrature::CWeightLocation>::const_iterator it = a().begin();
		 it != a().end(); ++ it)
	{
		CWeightLocation arg;
		const CFixedPrecisionGaussianQuadrature::CHPNumber& value = it->first;
		const CFixedPrecisionGaussianQuadrature::CHPNumber& weight = it->second;
		arg.displacement = value.toDouble();
		arg.probability = weight.toDouble();
		vlwNormalApproxMean0Var1.push_back(arg);
	}
}
#endif

void Normal01QuadratureDoublePrecision(CVectorWeightLocation& vlwNormalApproxMean0Var1, const unsigned int iNoPoints)
{
	// create the "double" quadrature
	vlwNormalApproxMean0Var1.clear();

	std::stringstream filename;
	filename << "QuadraturePoints"<< iNoPoints <<".txt";
	std::ifstream inFile(filename.str().c_str());
	if(!inFile) {
		std::cout<< "File does not exist!" << std::endl;
		exit(1);
	}

	double displacement(0), probability(0);
	unsigned int countInputSize(0);
	while(!inFile.eof()) {
		inFile >> std::hex >> displacement >> probability;
		countInputSize++;
		CWeightLocation arg;
		arg.displacement = displacement;
		arg.probability = probability;
		vlwNormalApproxMean0Var1.push_back(arg);
	}
	
	vlwNormalApproxMean0Var1.pop_back();	// remove duplicate --- how to get rid off this?
	countInputSize--;

	// Size check.
	if(countInputSize != iNoPoints) {
		std::cout << "Inconsistent number of QuadraturePoints" << std::endl;
		exit(1);
	}

	inFile.close();
}

CDoublePrecisionGaussianQuadrature::CDoublePrecisionGaussianQuadrature(void) {}

CDoublePrecisionGaussianQuadrature::~CDoublePrecisionGaussianQuadrature(void) {}
