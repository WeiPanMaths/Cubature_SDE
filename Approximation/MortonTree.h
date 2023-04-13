#pragma once

#include<map>
#include<memory>
#include"Utils.h"

class MortonTreeApproximation;
class MortonTreeOccupant;

class MortonTree : private std::map<MortonTreeNode, std::shared_ptr<MortonTreeOccupant>>
{
	MortonTreeApproximation* pApprox;

public:
	MortonTree(void);
	~MortonTree(void);

	void subdivideNode(MortonTree::iterator&);
	void subdivideNode(MortonTreeNode&);
	iterator addNode(const MortonTreeNode&, const_iterator&, size_t&);
	iterator addNode(const MortonTreeNode&, size_t&);
	std::shared_ptr<MortonTreeOccupant> query(const MortonTreeNode&);
	
	friend class MortonTreeOccupant;
	friend class MortonTreeApproximation;

private:
	size_t getLowestCommonAncestor(const MortonTreeNode&, const MortonTreeNode&,size_t) const;
};