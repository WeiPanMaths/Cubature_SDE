/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#include <vector>
#include <math.h>
//************************************
// Method:    ErrorInCubature
// FullName:  ErrorInCubature
// Access:    public 
// Returns:   double
// Qualifier:
// Parameter: double dNumberOfCubaturePoints
// Parameter: double dTimeToBoundaryOfCubaturePoints
// Parameter: double dTimeToCubaturePoints
//************************************
double ErrorInCubature(double dNumberOfCubaturePoints, double dTimeToBoundaryOfCubaturePoints,
					   double dTimeToCubaturePoints)
{
	//f(t) = (f(0) + f1(0) h + f2(0)/2! h^2 + ... + \int 0_t fn(u)(t-u)^n-1 /(n-1)! du
	//So p(t)f - qf <= sup fn ^n...
	// in fact the  n-poly degree cubature over an interval of length s against P_t(g) where  g=(1-e^x)+ is at most
	// error(n, s, t) <= 2 * (e/dPi)**(1/4) * ((8 * n + 1) / (8 * e * t))**((2 * n + 1)/4) * (s/2) **((n+1)/2) /((n+1)/2)! 
	const double dPi = asin(1.);
	const double dE = exp(1.);
	double error = 2 *
		pow(dE / dPi, 1 / 4) *
		pow((8 * dNumberOfCubaturePoints + 1) / (8 * dE * dTimeToBoundaryOfCubaturePoints),
		((2 * dNumberOfCubaturePoints + 1) / 4)) *
		pow(dTimeToCubaturePoints / 2,
		(dNumberOfCubaturePoints + 1) / 2) /
		((dNumberOfCubaturePoints + 1) / 2);	
	return error;
}

//************************************
// Method:    TimeStepForAcceptableError
// FullName:  TimeStepForAcceptableError
// Access:    public 
// Returns:   double
// Qualifier:
// Parameter: const double & dCurrentTime
// Parameter: const double & dAcceptableError
// Parameter: const int & degree
//************************************
double TimeStepForAcceptableError(const double& dCurrentTime, const double& dAcceptableError, const int& degree)
{
	double t(dCurrentTime);
	double s = t / 2;
	t -= s;
	double n = degree;
	double dEstimatedError = ErrorInCubature(n, t, s);
	while (dEstimatedError > dAcceptableError)
	{
		double temp = s / 2;
		s -= temp;
		t += temp;
		dEstimatedError = ErrorInCubature(n, t, s);
	}
	return s;
}
