////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef bitutilities_h__
#define bitutilities_h__
#include "emmintrin.h"

// Some templates to construct various complex compile time integer constants

// Ones<TargetType, NoOnes>::ans is a const of type TargetType whose lowest NoOnes bits are 1 and the rest are zeros
template <class TargetType, unsigned NoOnes>
struct Ones {
	enum : TargetType
	{
		ans = (TargetType(Ones<TargetType, NoOnes - 1>::ans) << 1) | TargetType(1)
	};
};

// Ones<TargetType, NoOnes>::ans is a const of type TargetType whose lowest NoOnes bits are 1 and the rest are zeros
template <class TargetType>
struct Ones<TargetType, 0> {
	enum : TargetType
	{
		ans = TargetType(0)
	};
};

// repeat<TargetType, NoRepetitions, BitWidth, BitPattern>::ans is a const of type TargetType 
// which replicates the lowest BitWidth bits of BitPattern NoRepetitions times into a zero intialized variable
template <class TargetType, unsigned NoRepetitions, unsigned BitWidth,
	TargetType BitPattern>
struct repeat {
	enum : TargetType
	{
		// shift the intermediate ans  by BitWidth bits, zero BitPattern off BitWidth bits
		// insert BitPattern into the lowest bits of the ans and pass ans up for further modification 
		ans = ((TargetType(repeat <TargetType, NoRepetitions - 1, BitWidth, BitPattern>::ans) << BitWidth) | (TargetType(BitPattern) & TargetType(Ones<TargetType, BitWidth>::ans)))
	};
};

template <class TargetType, unsigned BitWidth, TargetType BitPattern>
struct repeat<TargetType, 1, BitWidth, BitPattern> {
	enum : TargetType
	{
		// zero BitPattern off BitWidth bits and return it  
		ans = (TargetType(BitPattern) & TargetType(Ones<TargetType, BitWidth>::ans))
	};
};

template <class TargetType, unsigned BitWidth, TargetType BitPattern>
struct repeat<TargetType, 0, BitWidth, BitPattern> {
	enum : TargetType
	{
		ans = TargetType(0)
	};
};

// expand<TargetType, NoRepetitions, BitWidth, BitPattern>::ans is a const of type TargetType 
// takes BitPattern truncated after BitWidth and replaces each bit in [0,...,BitWidth) 
// by NoRepetitions identical bits so that ans is a const TargetType bit pattern 
// with BitWidth*NoRepetitions bits extended by zeros
template <class TargetType, unsigned NoRepetitions, unsigned BitWidth,
	TargetType BitPattern>
struct expand {
	enum : TargetType
	{
		field = (BitPattern%2) ? TargetType(Ones<TargetType, NoRepetitions>::ans) : TargetType(0),
		ans = ((TargetType(expand <TargetType, NoRepetitions, BitWidth - 1, (BitPattern >> 1)>::ans) << NoRepetitions) | TargetType(field))
	};
};

template <class TargetType, unsigned NoRepetitions, TargetType BitPattern>
struct expand<TargetType, NoRepetitions, 0, BitPattern> {
	enum : TargetType
	{
		ans = TargetType(0)
	};
};

template <class TargetType, unsigned PatternWidth>
struct splice_int_data 
{
	enum : TargetType
	{
		Down = repeat<TargetType, (sizeof(TargetType) * 8) / PatternWidth, PatternWidth, expand<TargetType, PatternWidth / 4, 4, 4>::ans>::ans,
		Up =  repeat<TargetType, (sizeof(TargetType) * 8) / PatternWidth, PatternWidth, expand<TargetType, PatternWidth / 4, 4, 2>::ans>::ans,
		Static =  repeat<TargetType, (sizeof(TargetType) * 8) / PatternWidth, PatternWidth, expand<TargetType, PatternWidth / 4, 4, 9>::ans>::ans
	};	
};

template <unsigned PatternWidth>
struct shuffle{
	// interleaves the bits in the top and bottom parts of successive fields of length PatternWidth abAB->aAbB in arg 
	// leaving the top bit unmoved in each field
	template <class TargetType > //
	inline TargetType operator() (TargetType arg) const
	{
		static_assert((sizeof(TargetType)*8)%PatternWidth == 0 && (sizeof(TargetType)*8) / PatternWidth > 0,  "arg not a positive exact number of patterns");
		const static unsigned bits = (PatternWidth>>2);
		const TargetType & Down = splice_int_data<TargetType, PatternWidth>::Down;
		const TargetType & Up = splice_int_data<TargetType, PatternWidth>::Up;
		const TargetType & Static = splice_int_data<TargetType, PatternWidth>::Static;
		return shuffle<(PatternWidth >> 1)>() ((arg & Static) | ((arg << bits) & Down) | ((arg >> bits) & Up));
	}
};


template <>
struct shuffle <unsigned(2)> {
	template <class TargetType > //
	inline TargetType operator() (TargetType arg) const
	{
		return arg;
	}
};

template <class LongType, class ShortType >
LongType nsplice(ShortType lhs, ShortType rhs)
{
	static_assert((sizeof(ShortType)<<1) == sizeof(LongType), "size of type to splice must be exactly half of size of the type output");
	LongType ans(lhs);
	ans = (ans << (sizeof(ShortType)*8)) | LongType(rhs);
	return shuffle<(sizeof(LongType)*8)>()(ans);
}

inline __m128i nsplice(unsigned __int64 lhs, unsigned __int64 rhs);
//__m128 _mm_shuffle_ps(__m128 a , __m128 b , int i );



#endif // bitutilities_h__

