#ifndef recombine_h__
#define recombine_h__

#ifdef RECOMBINE_EXPORTS
#define RECOMBINE_API __declspec(dllexport)
#else
#define RECOMBINE_API __declspec(dllimport)
#endif

#include "BufferConstructor.h"

/*
This header file defines the interface for programs seeking to do data reduction. 

The key input data are 
1) a set of N1 distinct points
	identified only by pointers to structures of type void
	and associated weights (positive doubles); 
2) d
3) a function that converts a point 
into a vector of d doubles and which will know what type to convert the void pointers to
it can deduce d from the size of the offered buffer

The key output data
1) a set of N2 <= d+1 distinct non-negative integers < N1, marking the points we retain and 
2) an equal number of revised weights so that the expectations of the marked points and new 
weights agree with the old expectation

We offer one interface
 a pure C interface
*/


#ifdef __cplusplus
extern "C"
{
#endif
	void RECOMBINE_API Recombine(void* recombine_interface);
#ifdef __cplusplus
} // extern "C" 
#endif

// the current interface uses
struct sCloud;
struct sRCloudInfo;
struct sRecombineInterface
{
	sCloud* pInCloud;
	sRCloudInfo* pOutCloudInfo;
	size_t degree;
	expander fn;
	void* end;
};
// where end must be null
// a C structure that points to the locations and weights in the cloud
struct sCloud
{
	size_t NoActiveWeightsLocations;
	double* WeightBuf;
	void* LocationBuf;
	void* end;
};

// a C structure for the returned information used for modifying the cloud
struct sRCloudInfo
{
	size_t No_KeptLocations; // number actually returned
	double* NewWeightBuf;  // a buffer containing the weights of the kept Locations // capacity must be at least degree + 1
	size_t* KeptLocations; // a buffer containing the offsets of the kept Locations // capacity must be at least degree + 1
	void* end;
};
// and in each case end must be null;
//



#endif // recombine_h__
