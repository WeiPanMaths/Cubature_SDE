////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#include "stdafx.h"
unsigned short splice(const unsigned char lhs, const unsigned char rhs)
{
	return
		(unsigned short (rhs & (unsigned char)1) << 0)
		| (unsigned short (lhs & (unsigned char)1) << 1)
		| (unsigned short (rhs & (unsigned char)2) << 1)
		| (unsigned short (lhs & (unsigned char)2) << 2)
		| (unsigned short (rhs & (unsigned char)4) << 2)
		| (unsigned short (lhs & (unsigned char)4) << 3)
		| (unsigned short (rhs & (unsigned char)8) << 3)
		| (unsigned short (lhs & (unsigned char)8) << 4)
		| (unsigned short (rhs & (unsigned char)16) << 4)
		| (unsigned short (lhs & (unsigned char)16) << 5)
		| (unsigned short (rhs & (unsigned char)32) << 5)
		| (unsigned short (lhs & (unsigned char)32) << 6)
		| (unsigned short (rhs & (unsigned char)64) << 6)
		| (unsigned short (lhs & (unsigned char)64) << 7)
		| (unsigned short (rhs & (unsigned char)128)<< 7)
		| (unsigned short (lhs & (unsigned char)128)<< 8);
}

unsigned char lunsplice(const unsigned short arg)
{
	return
		(unsigned char)(arg >> 1) & (unsigned char)1
		| (unsigned char)(arg >> 2) & (unsigned char)2
		| (unsigned char)(arg >> 3) & (unsigned char)4
		| (unsigned char)(arg >> 4) & (unsigned char)8
		| (unsigned char)(arg >> 5) & (unsigned char)16
		| (unsigned char)(arg >> 6) & (unsigned char)32
		| (unsigned char)(arg >> 7) & (unsigned char)64
		| (unsigned char)(arg >> 8) & (unsigned char)128;
}

unsigned char runsplice(const unsigned short arg)
{
	return
		(unsigned char)(arg >> 0) & (unsigned char)1
		| (unsigned char)(arg >> 1) & (unsigned char)2
		| (unsigned char)(arg >> 2) & (unsigned char)4
		| (unsigned char)(arg >> 3) & (unsigned char)8
		| (unsigned char)(arg >> 4) & (unsigned char)16
		| (unsigned char)(arg >> 5) & (unsigned char)32
		| (unsigned char)(arg >> 6) & (unsigned char)64
		| (unsigned char)(arg >> 7) & (unsigned char)128;
}

bool testsplice()
{
	unsigned char i, j;
	bool ans = true;
	for (i = 0; i < 256; ++i)
		for (j = 0; j < 256; ++j) {
			ans = ans && ((i == lunsplice(splice(i,
				j))) && (j == runsplice(splice(i, j))));
		}
	return ans;
}
