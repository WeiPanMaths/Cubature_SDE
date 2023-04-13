/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#include "stdafx.h"
#include ".\TreeBufferHelper.h"
#include <algorithm> //min
#include <crtdbg.h> //_ASSERT

CTreeBufferHelper::CTreeBufferHelper(size_t SmallestReducibleSetSize, size_t NoPointsToBeprocessed)
	: iNoTrees(SmallestReducibleSetSize),
	  iInitialNoLeaves(NoPointsToBeprocessed)
{
	_ASSERT(iInitialNoLeaves >= iNoTrees && iNoTrees > 0);
}

bool CTreeBufferHelper::isleaf(const size_t& node) const
{
	return (node < iInitialNoLeaves && node >= 0);
}

size_t CTreeBufferHelper::end() const
{
	return 2 * iInitialNoLeaves - iNoTrees;
}

bool CTreeBufferHelper::isnode(const size_t& node) const
{
	return node >= 0 && node < end();
}

size_t CTreeBufferHelper::parent(const size_t& node) const
{
	_ASSERT(isnode(node));
	return std::min(iInitialNoLeaves + (node / 2), end());
}

bool CTreeBufferHelper::isroot(const size_t& node) const
{
	_ASSERT(isnode(node));
	return parent(node) == end();
}

size_t CTreeBufferHelper::left(const size_t& node) const
{
	_ASSERT(isnode(node));
	// returns negative if leaf
	return (node - iInitialNoLeaves) * 2;
}

size_t CTreeBufferHelper::right(const size_t& node) const
{
	return (left(node) < 0) ? -1 : left(node) + 1;
}
