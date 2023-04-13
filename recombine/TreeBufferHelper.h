/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once

struct CTreeBufferHelper
{
	// the number of trees in the initial forest
	size_t iNoTrees;
	// the number of leaves in the initial forest
	size_t iInitialNoLeaves;
	// vdBuffer[iIndex +  iNoPointsToBeprocessed]
	// = vdBuffer[ 2 * iIndex ] + vdBuffer[ 2 * iIndex + 1 ] ;

	CTreeBufferHelper(size_t SmallestReducibleSetSize, size_t NoPointsToBeprocessed);				
	bool isleaf(const size_t& node)const;
	size_t end() const;
	bool isnode(const size_t& node) const;
	size_t parent(const size_t& node) const;
	bool isroot(const size_t& node) const;
	size_t left(const size_t& node) const;
	size_t right(const size_t& node) const;
};
