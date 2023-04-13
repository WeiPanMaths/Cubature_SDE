#pragma once
#ifndef wpUtilities_h__
#define wpUtilities_h__

#include "morton.h"
#include <array>

namespace wpUtilities
{
	template<size_t Dimension>
	using Point_POD = double[Dimension];

	//template<size_t Dimension>
	//using BLOCK = morton::block<Point_POD<Dimension>>::BLOCK;

	//template<size_t Dimension>
	//using UNBLOCK = morton::block<Point_POD<Dimension>>::UNBLOCK;

	//template<size_t Dimension>
	//using COMPARE = morton::block<Point_POD<Dimension>>::compare;

	template<size_t Dimension>
	using Point = std::array<double, Dimension>; 
}

//namespace wpUtilities2D
//{
//	typedef double Point_POD[2];
//	typedef morton::block<Point_POD>::BLOCK BLOCK;
//	typedef morton::block<Point_POD>::UNBLOCK UNBLOCK;
//	typedef morton::block<Point_POD>::compare COMPARE;
//	typedef std::array<double, 2> Point;
//}

#endif // wpUtilities_h__

