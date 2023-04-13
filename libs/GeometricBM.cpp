#include "GeometricBM.h"
#include <cmath>
/*
double GeometricBM(double s, double r, double sigma, double t, double bm)
{
	return s * exp( (r - .5*sigma*sigma)*t + sigma * bm );
}

double GeometricBM(double s, double drift, double bm)
{
	return s * exp( drift + bm);
}
*/
double GeometricBM(double s, double r, double sigma, double t, double bm)
/* s is log spot */
{
	//return exp(s + (r - .5*sigma*sigma)*t + sigma * bm);
	return s * exp((r - .5*sigma*sigma)*t + sigma * bm );
//	return s + bm;
}

double GeometricBM(double s, double drift, double bm)
/* s is log spot */
{
	//return exp(s + drift + bm);
	return s* exp(drift + bm);
//	return s + bm;
}