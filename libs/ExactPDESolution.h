/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once 

#ifndef ExactPDESolution_h__
#define ExactPDESolution_h__

double ExactPDESolution(double y, double t);							// BM
double ExactPDESolution(double x, double t, double r, double sigma);	// overloaded BM with drift
double ExactPDESolution(double spot, double strike, double timeToMat, double rate, double vol);	// BS call
double ExactPDESolution(double spot, double strike, double barrier, double timeToMat, double rate, double vol); // BS atExpCall

#endif // ExactPDESolution_h__
