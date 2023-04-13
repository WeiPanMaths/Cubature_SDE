////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef lsplice_h__
#define lsplice_h__
#include <type_traits>
#include "bitutilities.h"
#include <emmintrin.h>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned short splice(const unsigned char lhs, const unsigned char rhs);
unsigned char lunsplice(const unsigned short arg);
unsigned char runsplice(const unsigned short arg);

template <class tt, class T>
	T& lsplice(T& Ans, const tt& Arg)
	{
#define INTCHOICE 2

#if INTCHOICE==3
		typedef __m128i LONGT;//short LONGT; //
		typedef unsigned __int64 SHORTT; //char SHORTT; //
#define SPLICE nsplice<LONGT, SHORTT>
#endif
#if INTCHOICE==2
		typedef unsigned long long LONGT;//short LONGT; //
		typedef unsigned long SHORTT; //char SHORTT; //
#define SPLICE nsplice<LONGT, SHORTT>
#endif
#if INTCHOICE==1
		// the fastest in empirical tests
		typedef unsigned long LONGT; //short LONGT;
		typedef unsigned short SHORTT; //char SHORTT;
#define SPLICE nsplice<LONGT, SHORTT>
#endif
#if INTCHOICE==0
		// the fastest in empirical tests
		typedef unsigned short LONGT; //short LONGT;
		typedef unsigned char SHORTT; //char SHORTT;
#define SPLICE nsplice<LONGT, SHORTT>
//#define SPLICE splice
#endif

	//reinterpret types to arrays of bytes and shorts 
		static_assert((sizeof(Arg) == sizeof(Ans)) && (sizeof(Arg) % sizeof(LONGT) == 0), "lsplice: assumptions about data structure sizes are not true");
		typedef LONGT TBlock [sizeof(Arg)/sizeof(LONGT)];
		typedef SHORTT  tBlock [sizeof(Arg)/sizeof(LONGT)];
		typedef tBlock         ttBlock[2];

		TBlock& Zipped = reinterpret_cast <TBlock&> (Ans);
		const ttBlock& Split = reinterpret_cast <const ttBlock&> (Arg);

		std::transform (std::begin(Split[0]), std::end(Split[0]), std::begin(Split[1]), std::begin(Zipped), SPLICE);
		return reinterpret_cast<T&>(Zipped);
	};

template <class tt, class T>
	tt& lusplice(tt& Ans, const T& Arg)
	{
	//reinterpret types to arrays of bytes and shorts 
		static_assert((sizeof(Arg) == sizeof(Ans)) && (sizeof(Arg) % sizeof(unsigned short) == 0), "lsplice: assumptions about data structure sizes are not true");
		typedef unsigned short TBlock [sizeof(Arg)/sizeof(unsigned short)];
		typedef unsigned char  tBlock [sizeof(Arg)/sizeof(unsigned short)];
		typedef tBlock         ttBlock[2];

		const TBlock& Zipped = reinterpret_cast <const TBlock&> (Arg);
		ttBlock& Split = reinterpret_cast <ttBlock&> (Ans);

		typedef decltype(std::_Unchecked(std::begin(Zipped))) T_iterator;
		typedef decltype(std::_Unchecked(std::begin(Split[0]))) t_iterator;

		T_iterator NextIn(std::begin(Zipped)), Last(std::end(Zipped));
		t_iterator Dest1(std::begin(Split[0])), Dest2(std::begin(Split[1]));

		for (; NextIn != Last; ++NextIn, ++Dest1, ++Dest2)
		{
			*Dest1 = lunsplice(*NextIn);
			*Dest2 = runsplice(*NextIn);
		}
		return reinterpret_cast<tt&>(Split);
	};

#endif // lsplice_h__
