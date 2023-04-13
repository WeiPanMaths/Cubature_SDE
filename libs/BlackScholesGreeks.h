#pragma once
#include "CloudWeightLocation.h"
#include "MainProgramParameters.h"


double BSDeltaCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon);
double BSDeltaExactSolution(CMainProgramParameters& ipIntParams);
double BSDeltaFivePointCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon);
double BSDeltaFwdDifferenceCubatureApprox(CMainProgramParameters& ipIntParams,const CVectorWeightLocation & vlwNormalApproxMean0Var1, double epsilon);