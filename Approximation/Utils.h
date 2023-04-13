#pragma once

#ifndef Utils_h__
#define Utils_h__

#include <vector>
#include <array>
#include "wpUtilities.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <exception>

#define RMSE 0
#define _normalisation_ 1    // controls if we need to renormalise points using sigmoid in adaptive approximation

class myexception : public std::exception
{};

using boost::math::erfc;
using boost::math::erf;
//#define EXPONENTIATE				/// use exponential to transform all coordinates to the positive axis

const bool LEFTNBR  = true;
const bool RIGHTNBR = false;
const bool NOTENLARGEABLE = false;
const bool ENLARGEABLE = true;

const size_t DEFAULTDEPTH = 52;   // 52 TREATS [0.5, 1)^d as a single element // this is the "highest" fraction bit
const size_t DEFAULT_RES = 52; // sizeof(double)* CHAR_BIT - 1;  // == 63
const size_t DEFAULTMAXDEPTH = 35;  // do not subdivde beyond this

typedef bool		RejectionFlag;
typedef bool		LeftOrRight;

typedef std::vector<double>	vectordouble;

#define dim_ 2

//typedef double Point_POD[dim_];							// for testing purpose change to 1d
//typedef morton::block<Point_POD>::BLOCK BLOCK;
//typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
//typedef morton::block<Point_POD>::compare COMPARE;
//typedef		BLOCK		MortonTreeNode;

typedef wpUtilities::Point_POD<dim_> Point_POD;
//typedef wpUtilities::BLOCK<dim_> BLOCK;
//typedef wpUtilities::UNBLOCK<dim_> UNBLOCK;
//typedef wpUtilities::COMPARE<dim_> COMPARE;
typedef morton::block<Point_POD>::BLOCK BLOCK;
typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
typedef morton::block<Point_POD>::compare COMPARE;
typedef wpUtilities::Point<dim_> Point;

#if dim_ == 2
typedef Point_POD point_2d_POD;
#endif

#if dim_ == 1
typedef Point_POD point_POD; // for testing purpose change to 1d
#endif

typedef std::array<double, 1>	point;				// for interpolation
typedef std::array<double, 2>	point_2d;			// for interpolation
typedef		BLOCK		MortonTreeNode;

//#ifdef MKL64
//
//#pragma comment(lib, "mkl_intel_lp64.lib")
//#pragma comment(lib, "mkl_intel_thread.lib")
//#pragma comment(lib, "mkl_core.lib") 
//#pragma comment(lib, "libiomp5md.lib")
//
//#define dgels_ DGELS
//
//#endif

typedef double doublereal;
typedef int integer;

//minimize 2-norm(| b - A*x |)  least squares
//extern "C"
//void dgels_(char* trans, integer* m, integer* n, integer* nrhs, doublereal* a, 
//			integer* lda, doublereal* b, integer* ldb,
//			doublereal* work, integer* lwork, integer* info);


template<typename InputIterator>
double bitstring_to_double(InputIterator begin, InputIterator end)
{
	uint64_t x = 0;
	for (; begin != end; ++begin)
	{
		x = (x << 1) + (*begin - '0');
	}
	double d;
	memcpy(&d, &x, 8);
	return d;
}

struct FunctionBase {

	template<typename D, typename PointType>
	double operator() (const PointType& pt, D & data) const
	{ return this->operator()(pt); }

	virtual double operator()(const point_2d& pt) const
	{
		////return pt[0]*pt[0]*pt[0]*pt[0]+pt[1]*pt[1] + pt[0]*pt[1];
		////return (pt[0]<0.1)? pt[0]*pt[0]+pt[1]*pt[1] + pt[0]*pt[1] : -1.;
		////return 1.; //pt[0]*pt[0]; 

		//double s = pt[0];
		//double m = pt[1];
		//double k = 1.;			// TEST value: strike K = 1
		//if (m <0) return 0.;
		//else 
		//{
		//	double log_k_m = std::log(k+m);
		//	double one_over_root_2 = 1./std::sqrt(2.);
		//	double ans = std::exp(s + .5) * (1. + erf( (log_k_m - s - 1.)*one_over_root_2 ) );
		//	ans	-= k * (1. + erf( (log_k_m - s )*one_over_root_2 ) );
		//	ans += m * erfc( (log_k_m - s) * one_over_root_2 );
		//	ans *= .5;
		//	return ans;
		//}
		return 0.;
	}
	virtual double operator()(const point& pt) const
	{
		return 0.;
	}
};

/// use this so we don't need to template MortonTreeApproximation class
template <typename Function, typename PointType>
struct FunctionWrapper : public FunctionBase
{
	Function* pFunc;
	FunctionWrapper(Function* func) : pFunc(func) {}
	double operator()(const PointType& pt) const { return (*pFunc)(pt); }
};

template <typename Function>
struct FuncWrapper1D
{
	Function* pFunc;
	FuncWrapper1D(Function* func) : pFunc(func) {}
	double operator()(const point& pt) const { return (*pFunc)(pt); }
};


#define dlength sizeof(double)*CHAR_BIT

double getIntervalLength_(size_t bit);

// maps R to [0, 1]
inline double sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

inline double sigmoid_inv(double x)
{
	return log(x / (1.0 - x));
}

// transforms to [1, 2] to avoid near 0 values
inline double sigmoid_transform(double x)
{
	return 1.0 / (1.0 + exp(-x)) + 1.0;
}

inline double sigmoid_inv_transform(double x)
{
	double p = x - 1.0;
	return log(p / (1.0 - p));
}

inline double chebyshev_interval_shift(double x, double low, size_t depth)
{
	return  (x - low)* (1 << (DEFAULT_RES - depth + 1)) - 1.;
}


//template <typename F>
//class InterpolationMorton;
//class InterpolationMorton				// defined this for testing purposes
//{
//	BLOCK pt;
//public:
//	InterpolationMorton(BLOCK _pt) : pt(_pt) {}
//	~InterpolationMorton(){}
//	BLOCK getPT() const { return pt; }
//	double approximateFunction( const point_2d & x ) {return 1.;}
//};

//struct InterpolationMortonCompare
//{
//	bool operator()(const InterpolationMorton & l , const InterpolationMorton & r)
//	{
//		if (l.getPT() < r.getPT())
//			return true;
//		else
//			return false;
//	}
//};

//struct TestFunction{
//	double operator() ( const point_2d & pt )
//	{
//		double x = std::log(pt[0]);
//		double y = std::log(pt[1]);
//		double ans = -0.834 -0.25*x + 0.25* y  + 0.75*x*x + 0.375*x*y +0.75*y*y;
//		return ans;
//	}
//};

#endif // Utils_h__