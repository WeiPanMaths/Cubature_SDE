/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once
#ifndef Cmove1_h__
#define Cmove1_h__

#include <vector>

class CLinearAlgebraReductionTool
{
public:
	CLinearAlgebraReductionTool(void);
	~CLinearAlgebraReductionTool(void);

	inline int INoCoords() const
	{
		return iNoCoords;
	}
	inline const int& INoCoords(int val)
	{
		return iNoCoords = val;
	}
	inline int INoPoints() const
	{
		return iNoPoints;
	}
	inline const int& INoPoints(int val)
	{
		return iNoPoints = val;
	}
	inline unsigned int INoCallsLinAlg() const
	{
		return iNoCallsLinAlg;
	}
	void MoveMass(std::vector<double>& weights, const std::vector<double>& points, std::vector<double>& MassCog,
				  std::vector<int>& maxset);

private:
	int iNoCoords;
	int iNoPoints;
	int iNoRhs;

	// counts the number of calls to the linear reduction package
	unsigned int iNoCallsLinAlg;
	
};
#endif // Cmove1_h__
