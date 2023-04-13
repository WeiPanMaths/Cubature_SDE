#include "MortonTreeOccupant.h"
#include "MortonTreeApproximation.h"

MortonTreeOccupant::MortonTreeOccupant(size_t depth = DEFAULT_RES)
	: mDepth(depth), pInterpolation(nullptr), pTreeNode(),pTree(0)// ,_TEST_NODE_X(0),_TEST_NODE_Y(0)
{
}

MortonTreeOccupant::~MortonTreeOccupant(void)
{
}

bool MortonTreeOccupant::getLocalApprox(const BLOCK& pt)		
{
	if(pInterpolation) return true;
	else return false;
	/*else  {
		auto lower = pTree->pApprox->getLowerbound(pt,COMPARE(mDepth)); 
		auto upper = pTree->pApprox->getUpperbound(pt,COMPARE(mDepth));
		auto deg = pTree->pApprox->deg;
		auto value = distance(lower,upper);
		if(value>=(2+deg)*(1+deg)) 	return pTree->pApprox->tryInterpolation(lower,upper,this, pt);
		else						return false;
	}*/
}


bool MortonTreeOccupant::attemptApprox(const BLOCK& pt)
{
	auto lower = pTree->pApprox->getLowerbound(pt,COMPARE(mDepth)); 
	auto upper = pTree->pApprox->getUpperbound(pt,COMPARE(mDepth));
	auto deg = pTree->pApprox->deg;
	auto value = distance(lower,upper);

	/*if(value == 239)
	{
		auto pod_lower = point_2d_POD(UNBLOCK(lower->first));
		auto pod_upper = point_2d_POD(UNBLOCK(upper->first));
		double pod_lower_x = pod_lower[0];
		double pod_lower_y = pod_lower[1];
		double pod_upper_x = pod_upper[0];
		double pod_upper_y = pod_upper[1];
		std::cout << pod_lower_x << " , " << pod_lower_y << " ||| " << pod_upper_x << " , " << pod_upper_y << "\n";
	}*/

	if(value>=(2+deg)*(1+deg)) 	return pTree->pApprox->tryInterpolation(lower,upper,this, pt);
	else						return false;
}

//
//bool MortonTreeOccupant::getLocalApprox(const BLOCK& pt, bool& _badrange) // for testing
//{
//	if(pInterpolation) return true;
//	else  {
//		auto lower = pTree->pApprox->getLowerbound(pt,COMPARE(mDepth)); 
//		auto upper = pTree->pApprox->getUpperbound(pt,COMPARE(mDepth));
//		auto deg = pTree->pApprox->deg;
//		auto value = distance(lower,upper);
//		if(value>=(2+deg)*(1+deg))
//		{
//			if(_badrange) 
//				_badrange = true;
//			return pTree->pApprox->tryInterpolation(lower,upper,this, pt);
//		} else 
//			return false;
//	}
//
//}

double MortonTreeOccupant::getInterpApprox(const point_2d& pt) const
{
	//TEST
	pTree->pApprox->TEST_INC_NUM_APPROX_CALLS();
	return pInterpolation->approximateFunction(pt);
}	

size_t MortonTreeOccupant::getMortonLevel() const
{
	return mDepth;
}
