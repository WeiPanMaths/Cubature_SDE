#pragma once
#include "LeastSquaresInterpolation.h"
#include <memory>
#include "MortonTree.h"

class MortonTreeOccupant
{
private:
	MortonTree::iterator pTreeNode;
	size_t mDepth;
	std::shared_ptr<InterpolationMorton> pInterpolation;
	MortonTree* pTree;

	/*double _TEST_NODE_X;
	double _TEST_NODE_Y;*/

public:
	MortonTreeOccupant(size_t);
	~MortonTreeOccupant(void);

	bool	getLocalApprox(const BLOCK&);
	//bool	getLocalApprox(const BLOCK&, bool&);
	double	getInterpApprox(const point_2d&) const;
	size_t getMortonLevel() const;
	bool attemptApprox(const BLOCK&);

	//double getTestNodeX() const { return _TEST_NODE_X;}
	//double getTestNodeY() const { return _TEST_NODE_Y; }
	friend class MortonTree;
	friend class MortonTreeApproximation;
};