
/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once
//#include "SharedFunctionsAndTypes.h"
#include "CloudWeightLocation.h"

void Normal01QuadratureDoublePrecision(CVectorWeightLocation& vlwNormalApproxMean0Var1, const unsigned int iNoPoints =
									   10);

class CDoublePrecisionGaussianQuadrature
{
public:
	CDoublePrecisionGaussianQuadrature(void);
	~CDoublePrecisionGaussianQuadrature(void);
};
