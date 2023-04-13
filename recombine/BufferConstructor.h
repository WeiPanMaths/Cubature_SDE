#ifndef BufferConstructor_h__
#define BufferConstructor_h__

#ifdef __cplusplus
extern "C"
{
#endif

	// Call Back Function that builds the vector co-ordinates of points
	// ARG1 points (in some way) to an indexed list of NoPointsToBeprocessed points to be reduced, and is const and typically NoPointsToBeprocessed pointers
	// ARG2 points to buffer of doubles length NoPointsToBeprocessed *(SmallestReducibleSetSize - 1)
	// ARG3 contains the size information and, if appropriate, preconditioning information that the expander function understands
	// the first two elements in ARG3 and the implied dimensioned buffer in ARG2 have a fixed interpretation
	typedef void (*expander)(void*, double*, void*);

	// the structure used to inform the callback function
	struct CConditionedBufferHelper
	{
		size_t SmallestReducibleSetSize;
		size_t NoPointsToBeprocessed;
		void* pvCConditioning;
	};

	typedef CConditionedBufferHelper CBufferHelper;

	// an example of a conditioning that might be used in a given callback function
	//struct CConditioning
	//{
	//	double dMean;
	//	double dStdDev;
	//};

	//an example callback function
	// void ArrayOfpDoublesToVectorOfPowers1( void * pIn , double * pOut , void* vpCConditionedBufferHelper );

#ifdef __cplusplus
}
#endif


#endif // BufferConstructor_h__
