#include "MortonTree.h"
#include "MortonTreeApproximation.h"
#include "MortonTreeOccupant.h"

MortonTree::MortonTree(void) : pApprox(0)
{
}

MortonTree::~MortonTree(void)
{
}

void MortonTree::subdivideNode(MortonTree::iterator& node)
{
	auto olddepth = node->second->getMortonLevel();
	auto newdepth = olddepth-1;
	auto left = pApprox->getLowerbound(node->first,COMPARE(olddepth));
	auto rightmost = pApprox->getUpperbound(node->first,COMPARE(olddepth));
	//point_2d_POD left_pod = { point_2d_POD(UNBLOCK(  left->first  ))[0], point_2d_POD(UNBLOCK(left->first))[1]};
	//point_2d_POD right_pod = { point_2d_POD(UNBLOCK(rightmost->first))[0], point_2d_POD(UNBLOCK(rightmost->first))[1]};
	erase(node);
	//node->second->mDepth--;
	while(left != rightmost)
	{
		addNode(left->first,newdepth);
		left = pApprox->getUpperbound(left->first,COMPARE(newdepth));
		//left_pod[0] = point_2d_POD(UNBLOCK(  left->first  ))[0];
		//left_pod[1] = point_2d_POD(UNBLOCK(  left->first  ))[1];
	}
}

void MortonTree::subdivideNode(MortonTreeNode& node)
{
	subdivideNode(find(node));
}

auto MortonTree::addNode(const MortonTreeNode& node, const_iterator &pos, size_t& depth)->iterator
{
	auto it = insert(pos, std::pair<MortonTreeNode, std::shared_ptr<MortonTreeOccupant>>(node, std::make_shared<MortonTreeOccupant>(depth)));
	it->second->pTree = this;
	it->second->pTreeNode = it;
	
	/// testing ************************************
	/*it->second->_TEST_NODE_X = point_2d_POD(UNBLOCK(node))[0];
	it->second->_TEST_NODE_Y = point_2d_POD(UNBLOCK(node))[1];*/
	///************** end testing
	return it;
}

auto MortonTree::addNode(const MortonTreeNode& node, size_t& depth)->iterator
{
	return addNode(node, begin(), depth);
}

std::shared_ptr<MortonTreeOccupant> MortonTree::query(const MortonTreeNode& pt)
{
	size_t new_depth = DEFAULT_RES;
	const auto itpt = addNode(pt, new_depth);
	auto rightnbr = itpt; ++rightnbr;
	auto leftnbr = itpt; --leftnbr;
	auto itend = end();
	size_t lnbr_LCAdepth = DEFAULT_RES+1;	// root depth definition for determining forest. in this case it's 64
	size_t rnbr_LCAdepth = DEFAULT_RES+1;	// root depth definition
	if (leftnbr!=itend) 
	{
		auto ldepth = leftnbr->second->getMortonLevel();
		COMPARE lless(ldepth);
		if (!lless(leftnbr->first,pt))
		{
			erase(itpt);
			return leftnbr->second;
		}
		lnbr_LCAdepth = getLowestCommonAncestor(leftnbr->first,pt,ldepth);		
	} 
	if (rightnbr!=itend)
	{
		auto rdepth = rightnbr->second->getMortonLevel();
		COMPARE rless(rdepth);
		if (!rless(pt,rightnbr->first))
		{
			erase(itpt);
			return rightnbr->second;
		}
		rnbr_LCAdepth = getLowestCommonAncestor(pt,rightnbr->first,rdepth);
	}
	itpt->second->mDepth = std::min(lnbr_LCAdepth, rnbr_LCAdepth)-1; // wud there be a problem if both are 64
	return itpt->second;
}

// return the lowest level in (refDepth, DEFAULTDEAPTH] such that n1 and n2 are different but not bigger than refDepth. 
// if return 64 then that means n1 and n2 have no common root
size_t MortonTree::getLowestCommonAncestor(const MortonTreeNode& n1,const MortonTreeNode& n2,size_t refDepth) const
{
	//size_t max_res = sizeof(double)* CHAR_BIT - 1 ;		// for determining whether it's possible to have common ancestor; -1 is to discount sgn bit
	size_t max_res = DEFAULT_RES;
	bool notbigenough(true);
	COMPARE maxless(max_res);
	if ( maxless(n1, n2)) return max_res+1;			// no common ancestors, eg, in 1 D: 0.1 and -0.1 won't be in the same box
	--refDepth;
	while ( notbigenough && (++refDepth < max_res) )
	{
		COMPARE less(refDepth);
		notbigenough = less(n1, n2);		// true if A < B, hence A and B are not in the same box according to resolution temp_depth
	}
	return refDepth;
}