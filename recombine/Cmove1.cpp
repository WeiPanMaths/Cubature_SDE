/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#include "stdafx.h"
#include "Cmove1.h" //CLinearAlgebraReductionTool
#include <algorithm>
//#include <vector>
//#include <iostream>
#include <cmath>
#include <crtdbg.h>
#include <limits>

#ifdef MKL
#pragma comment(lib, "mkl_intel_c.lib")
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "mkl_core.lib") 
#pragma comment(lib, "libiomp5md.lib")

#include "aligned_allocator.h"

typedef std::vector<double, TJL_alloc::aligned_allocator<double, 16> > dVECTOR;
typedef std::vector<int, TJL_alloc::aligned_allocator<int, 16> > iVECTOR;

#define daxpy_ DAXPY
#define ddot_ DDOT
#define dgelsd_ DGELSD
#define dgelss_ DGELSS
#define dnrm2_ DNRM2
#define dgesv_ DGESV
#define dcopy_ DCOPY
#define dgemm_ DGEMM

#elif defined(MKL64)

#pragma comment(lib, "mkl_intel_lp64.lib")
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "mkl_core.lib") 
#pragma comment(lib, "libiomp5md.lib")

//#include "aligned_allocator.h"

//typedef std::vector<double, TJL_alloc::aligned_allocator<double, 16> > dVECTOR;
//typedef std::vector<int, TJL_alloc::aligned_allocator<int, 16> > iVECTOR;

typedef std::vector<double> dVECTOR;
typedef std::vector<int> iVECTOR;


#define daxpy_ DAXPY
#define ddot_ DDOT
#define dgelsd_ DGELSD
#define dgelss_ DGELSS
#define dnrm2_ DNRM2
#define dgesv_ DGESV
#define dcopy_ DCOPY
#define dgemm_ DGEMM

#else

#pragma comment(lib,  "liblapack.a")
#pragma comment(lib,  "libf77blas.a")
#pragma comment(lib,  "libcblas.a")
#pragma comment(lib,  "libatlas.a")

#pragma comment(lib,  "libg2c.a")
#pragma comment(lib,  "libgcc.a")

typedef std::vector< double > dVECTOR;
typedef std::vector< int > iVECTOR;

#endif

inline double dUniformRandomVariable(const double lower=double(-1.), const double upper = double(1.))
{
	// returns a uniform random variable on the interval [lower, upper]
	return lower + (upper - lower) * ((double)rand()) / ((double)RAND_MAX);
}

template <class T>
std::ostream& operator <<(std::ostream& os, const std::vector<T>& in)
{
	os << "{";
	for (typename std::vector<T>::const_iterator it(in.begin()); it != in.end(); ++it)
	{
		os << *it;
		if (it + 1 != in.end())
			os << ", ";
	}
	os << "}";
	return os;
}

using std::max;
using std::min;

typedef double doublereal;
typedef int integer;

// Y=a*X+Y
extern "C"
void daxpy_(integer* N, doublereal* A, doublereal* X, integer* INCX, doublereal* Y, integer* INCY);

// returns X.Y
extern "C"
double ddot_(integer* N, doublereal* X, integer* INCX, doublereal* Y, integer* INCY);

//minimize 2-norm(| b - A*x |)
extern "C"
void dgelsd_(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
			 doublereal*
			 s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* iwork, integer* info);
extern "C"
void dgelss_(integer* m, integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
			 doublereal*
			 s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork, integer* info);

extern "C"
double dnrm2_(int* n, double* x, int* incx);
extern "C"
void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
extern "C"
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
extern "C"
void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int
			* lda, double* b, int* ldb, double* beta, double* c, int* ldc);



//int iNoCoords=2; //M 
//int iNoPoints=3; //N
int iNoRhs = 1; //NRHS
// following not constant because the lapack programmes can change them (but dont....)
int iOne(1);
int iZero(0);
double dOne(1);
double dMOne(-1);
double dZero(0);
char charN('N');




CLinearAlgebraReductionTool::CLinearAlgebraReductionTool(void)
	: iNoCoords(1),
	  iNoPoints(1),
	  iNoRhs(1),
	  iNoCallsLinAlg(0) {}

CLinearAlgebraReductionTool::~CLinearAlgebraReductionTool(void) {}

//#define  _ORTHOG
//#define  _ROTATE

void CLinearAlgebraReductionTool::MoveMass
	(std::vector<double>& eweights, const std::vector<double>& epoints, std::vector<double>& eMassCog,
	 std::vector<int>& maxset)
{
	//	CDebugTools& debugtool = ipIntParams.dbgLogFile;
#ifdef _ROTATE
	dVECTOR dDiagonalMatrix;
#endif
	dVECTOR weights(eweights.begin(),eweights.end());
	dVECTOR points;
	for (int i1 = 0; i1 < iNoPoints; i1 ++)
	{
#ifdef _ROTATE
		double dFactor = dUniformRandomVariable(.5,1.5);
		weights[i1]/=dFactor;
		dDiagonalMatrix.push_back(dFactor);
#endif		
		for (int i2 = 0; i2 < iNoCoords; i2 ++)
		{
#ifdef _ROTATE
			points.push_back(dFactor*epoints[i2+i1*iNoCoords]);
#else
			points.push_back(epoints[i2 + i1 * iNoCoords]);
#endif
		}

#ifdef _ORTHOG
		points.push_back(weights[i1]);
#endif	

	}
	dVECTOR MassCog(eMassCog.begin(),eMassCog.end());

#ifdef _ORTHOG
	MassCog.push_back(0);
#endif	



	// need to solve Ax=Aw; then we use x-w, using geslv to find x is a stable process
	// and this works well if the matrix is sizeable, because the chances of x==w are very small;
	// one might think that taking x orthogonal to w would be a good idea, but this seems very unstable.
	// To test this w are the original weights, there are more weights than dimensions to points by 1.
	// need to append w as a row to A and add 0 to Aw at the same height then use dgeslv.
	//
	// better apply a random rescaling of each row
	dVECTOR& dvecInitWeights(weights);
	dVECTOR dvecPointSet(points);//A M=#rows  N=#cols
	dVECTOR& dvecMassCog = MassCog;

#ifdef _ORTHOG
	int iNoExtendedCoords = iNoCoords + 1;
#else
	int iNoExtendedCoords = iNoCoords;
#endif	


	dgemm_(&charN, &charN, &iNoExtendedCoords, &iOne, &iNoPoints, &dOne,
		&dvecPointSet[0], &iNoExtendedCoords, &dvecInitWeights[0], &iNoPoints,
		&dZero, &dvecMassCog[0], &iNoExtendedCoords);
	dVECTOR dvecRHStoLHS(dvecMassCog);
	dvecRHStoLHS.resize(max(max(iNoExtendedCoords, iNoPoints), 1)); //starts as b finishes as x nRHS times
	//DOUBLE PRECISION array, dimension (min(M,N)) 
	//containing The singular values of A in decreasing order
	//The condition number of A in the 2-norm = S(1)/S(min(m,n))
	dVECTOR dvecS(min(iNoExtendedCoords, iNoPoints));
	//Singular values S(i) <= RCOND*S(1) are treated as zero.
	//If RCOND < 0, machine precision is used instead
	double RCOND (-1);
	//double RCOND (1.e-20L);

	int nLDB1(max(iNoExtendedCoords, iNoPoints));
	int iRank;
	//initially an array of size Lwork (>=1), should be resized using the returned value of *work
	dVECTOR dvecWork(100);
	//set to query optimal for lenght of dvecWork
	int iLWork(-1);

	//dim >= 3 * MINMN * MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) + 11 * MINMN
	// dim >= MINMN * (11 + 3* MAX( 0, INT( LOG_2( MIN( M,N )) ) + 1 ))
	//fill with symbols to check it does not overflow
	//	iVECTOR iWork (10000);
	iVECTOR iWork (min(iNoExtendedCoords, iNoPoints) * (11 + 3 * max(0, (int)(log((double)min(iNoExtendedCoords,
																			  iNoPoints)) * (1 / log((double)2))) +
														1)),787654);
	
	//_ASSERT( ( (int)&dvecRHStoLHS[0]  %16 ) == 0 );
	//_ASSERT( ( (int)&dvecS[0]         %16 ) == 0 );
	//_ASSERT( ( (int)&dvecWork[0]      %16 ) == 0 );
	//_ASSERT( ( (int)&iWork[0]         %16 ) == 0 );
	//_ASSERT( ( (int)&dvecPointSet[0]  %16 ) == 0 );

	int iInfo;
	dgelsd_(&iNoExtendedCoords,
		&iNoPoints,
		&iNoRhs,
		&dvecPointSet[0],
		&iNoExtendedCoords,
		&dvecRHStoLHS[0],
		&nLDB1,
		&dvecS[0],
		&RCOND,
		&iRank,
		&dvecWork[0],
		&iLWork,
		&iWork[0],//put this back to change between dgelsd_ dgelss_ 
		&iInfo);
	iLWork = (int)dvecWork[0];//get optimal work size
	dvecWork.resize(iLWork);
	dgelsd_(&iNoExtendedCoords,
		&iNoPoints,
		&iNoRhs,
		&dvecPointSet[0],
		&iNoExtendedCoords,
		&dvecRHStoLHS[0],
		&nLDB1,
		&dvecS[0],
		&RCOND,
		&iRank,
		&dvecWork[0], //must not be too small even at the enquiry stage
		&iLWork,
		&iWork[0], //not used in small example////put this back to change between dgelsd_ dgelss_ 
		&iInfo);
	iNoCallsLinAlg++;

	dvecPointSet.assign(points.begin(), points.end());
	dVECTOR dvecMassMean1(max(iNoExtendedCoords, 1));
	dgemm_(&charN, &charN, &iNoExtendedCoords, &iOne, &iNoPoints, &dOne,
		&dvecPointSet[0], &iNoExtendedCoords, &dvecRHStoLHS[0], &iNoPoints,
		&dZero, &dvecMassMean1[0], &iNoExtendedCoords);
	// take the difference of RHS: mass and cogs
	daxpy_(&iNoExtendedCoords, &dMOne, &dvecMassCog[0], &iOne, &dvecMassMean1[0], &iOne);
	// take the difference of the LHS: weights and check orthogonality
	daxpy_(&iNoPoints, &dMOne, &dvecInitWeights[0], &iOne, &dvecRHStoLHS[0], &iOne);

	double ratio (0);
	double dRatio (- std::numeric_limits<double>:: infinity());
	maxset.clear();
	for (int i = 0; i != iNoPoints; i++)
	{
		// solve 0 <= dvecInitWeights + dRatio * dvecRHStoLHS 
		// for the smallest version attaining 0 at at least one i and recording all values i
		if (dvecInitWeights[i] + dRatio * dvecRHStoLHS[i] <= 0.)
		{
			if (dvecInitWeights[i] + dRatio * dvecRHStoLHS[i] == 0.)
				maxset.push_back(i);
			else
			{
				maxset.clear();
				maxset.push_back(i);
				_ASSERT(dRatio < - dvecInitWeights[i] / dvecRHStoLHS[i]);
				dRatio = - dvecInitWeights[i] / dvecRHStoLHS[i];
			}
		}

		//double temp;
		//if ((temp=(dvecRHStoLHS[i]/dvecInitWeights[i])) == ratio)
		//{
		//	maxset.push_back(i);
		//}
		//else
		//	if (temp > ratio)
		//	{
		//		maxset.clear();
		//		maxset.push_back(i);
		//		ratio = temp;
		//	}
	}

	//double dRatio = -1/ ratio;
	//std::cout << "The ratio is " << ratio << " and the max points are "<< maxset << std::endl<< std::endl;
	daxpy_(&iNoPoints, &dRatio, &dvecRHStoLHS[0], &iOne, &dvecInitWeights[0], &iOne);
#ifdef _ROTATE
	for (int i1 = 0; i1 < iNoPoints ; i1 ++ )
	{
		weights[i1]*=dDiagonalMatrix[i1];
	}
#endif

#ifdef _ORTHOG
	MassCog.pop_back();
#endif
	// return data
#ifdef MKL
	// using custom allocator internally so in general cannot swap
	eMassCog.assign(MassCog.begin(), MassCog.end());
	eweights.assign(weights.begin(), weights.end());

#elif defined(MKL64)

	// using custom allocator internally so in general cannot swap
	eMassCog.assign(MassCog.begin(), MassCog.end());
	eweights.assign(weights.begin(), weights.end());

#else
	eMassCog.swap(MassCog);
	eweights.swap(weights);
#endif
	
	for (std::vector<int>::const_iterator it = maxset.begin(); it != maxset.end(); ++it)
		eweights[*it] = 0;
}

#ifdef _ORTHOG
#undef _ORTHOG
#endif
#ifdef _ROTATE
#undef _ROTATE
#endif
