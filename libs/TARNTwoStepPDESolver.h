#pragma once

#ifndef __TARNTwoStepPDESolver_h__
#define __TARNTwoStepPDESolver_h__
//#include "SharedFunctionsAndTypes.h" // CVectorWeightLocation
#include "MainProgramParameters.h"
#include "CloudWeightLocation.h"
 
//solves the AtExpOpt
double TARNOneStepPDESolver(double spot, const CVectorWeightLocation & vlwNormalApproxMean0Var1, CMainProgramParameters & ipIntParams, double e );
double TestApproxAtExpiry(double x, unsigned int d, double int_err, CMainProgramParameters & ipIntParams );
double TestApproximationClass(double x, unsigned int d, double int_err);


#endif