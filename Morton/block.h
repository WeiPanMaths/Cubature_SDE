////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef block_h__
#define block_h__

#include <utility>
#include "lsplice.h"

namespace tjlblockdetail {

	template<class V_T, unsigned ZippedLogDim, unsigned UnZippedLogDim>
	union Block;

	template<class V_T>
	union Block<V_T, 0, 0>
	{
		// a one dimensional array needs no zipping
		V_T internal_V_T_array[1];

		bool operator < (const Block rhs) const {return internal_V_T_array[0] < rhs.internal_V_T_array[0];}
		bool operator == (const Block rhs) const{return internal_V_T_array[0] == rhs.internal_V_T_array[0];}

		// The format for the transformed block we will zip to
		V_T reblocked;

		// the code that splices pairs of spliced data.... until there is only one data block
		V_T& flip() {return reblocked;}
	};

	template<class V_T, unsigned UnZippedLogDim>
	union Block<V_T, 0, UnZippedLogDim>
	{
		// first member of Block must bit fill Block without any padding
		// since a union is (by definition) zero initialized by the first member being zero initialized
		V_T internal_V_T_array[1 << UnZippedLogDim];

		bool operator < (const Block rhs) const {
			typedef std::reverse_iterator<decltype(std::
				begin (internal_V_T_array))> R_It_t;
			R_It_t lhsREnd( std::begin (internal_V_T_array) );
			R_It_t lhsRBegin(std::end (internal_V_T_array));
			R_It_t rhsREnd(std::begin(rhs.internal_V_T_array));
			R_It_t rhsRBegin(std::end(rhs.internal_V_T_array));
			return std::lexicographical_compare(lhsRBegin, lhsREnd,
				rhsRBegin, rhsREnd);
		}

		bool operator == (const Block rhs) const
		{
			typedef std::reverse_iterator<decltype(std::
				begin (internal_V_T_array))> R_It_t;
			R_It_t lhsREnd( std::begin (internal_V_T_array) );
			R_It_t lhsRBegin(std::end (internal_V_T_array));
			R_It_t rhsRBegin(std::end(rhs.internal_V_T_array));
			return std::equal(lhsRBegin, lhsREnd, rhsRBegin);
		}

		// The format for the transformed block we will zip to
		Block<V_T, 1, UnZippedLogDim - 1 > reblocked;

		// the code that splices pairs of spliced data.... until there is only one data block
		Block<V_T, UnZippedLogDim, 0>& flip()
		{
			{
				auto itNextIn = std::begin(reblocked.data);
				auto itLastIn = std::end(reblocked.data);
				auto itNextOut = std::begin(reblocked.data);
				// working space for zipped output
				Block<V_T, 1, UnZippedLogDim -
					1 >::currently_zipped_as_1dim_t temp;
				for (; itNextIn != itLastIn; ++itNextIn, ++itNextOut) {
					lsplice(temp, *itNextIn);
					// in and out point to the same space
					auto itOut = std::begin(*itNextOut);
					std::copy(std::begin(temp), std::end(temp), itOut);
				}
			}
			return reblocked.flip();
		}
	};

	template<class V_T, unsigned ZippedLogDim, unsigned UnZippedLogDim>
	union Block
	{

		typedef V_T currently_zipped_as_1dim_t[1<<ZippedLogDim];
		currently_zipped_as_1dim_t data[1<<UnZippedLogDim];
		Block<V_T, ZippedLogDim + 1, UnZippedLogDim - 1> reblocked;

		Block<V_T, ZippedLogDim + UnZippedLogDim,  0 >& flip()
		{
			{
				auto itNextIn = std::begin(reblocked.data);
				auto itLastIn = std::end(reblocked.data);
				auto itNextOut = std::begin(reblocked.data);
				Block<V_T, ZippedLogDim + 1 , UnZippedLogDim -
					1>::currently_zipped_as_1dim_t temp;
				for (; itNextIn != itLastIn; ++itNextIn, ++itNextOut) {
					lsplice(temp, *itNextIn);
					// in and out point to the same space
					std::copy(std::begin(temp),
						std::end(temp), (*itNextOut));
				}
			}
			return reblocked.flip();
		}
	};

	template<class V_T, unsigned ZippedLogDim>
	union
		Block<V_T, ZippedLogDim, 0> {

			typedef V_T currently_zipped_as_1dim_t[1<<ZippedLogDim];
			// a 2 dimensional array
			currently_zipped_as_1dim_t data[1];
			// flip completed
			Block& flip()
			{
				return *this;
			}
	};
	// a function object chaining a unary operator that sets the last "bits" bits of each co-ordinate in a block to zero 
	// via an operator that sets the last "bits" bits in a monolithic int data stream to zero and is agnostic to type
	struct round_down_block {
		// the agnostic rounding tool
		tjlUtility::round_down_block_tjlU round;
		
		// the constructor for the function object
		round_down_block(const size_t& bits) : round(bits)
		{
		};

		template<class V_T, unsigned UnZippedLogDim>
		Block<V_T, 0, UnZippedLogDim> operator()(const Block<V_T, 0, UnZippedLogDim>& lhs)
		{
			// this conversion is messy as one side is a union containing the array 
			// the other side is the array
			static_assert(sizeof(tjlUtility::smallest_int_type<sizeof(V_T) << UnZippedLogDim>::int_t) == sizeof(lhs.internal_V_T_array), "size mismatch");
			auto& lh = (const tjlUtility::smallest_int_type< sizeof(V_T) << UnZippedLogDim >::int_t&) lhs.internal_V_T_array;
			Block<V_T, 0, UnZippedLogDim> ans;
			reinterpret_cast<tjlUtility::smallest_int_type<sizeof(lhs.internal_V_T_array)>::int_t&>(ans) = round.operator() <sizeof(V_T) << UnZippedLogDim> (lh);
			return ans;
		}
	};
	// a binary comparison operator that ignores the last few bits
	struct compare_less_block {
		tjlUtility::compare_less_tjlU less;
		compare_less_block(const size_t& bits) : less(bits)
		{
		};

		template<class V_T, unsigned UnZippedLogDim>
		bool operator()(const Block<V_T,
			0, UnZippedLogDim>& lhs, const Block<V_T,
			0, UnZippedLogDim>& rhs)
		{
			// this conversion is messy as one side is a union containing the array 
			// the other side is the array
			static_assert(sizeof(tjlUtility::smallest_int_type<sizeof(V_T) <<
				UnZippedLogDim>::int_t) == sizeof(lhs.
				internal_V_T_array), "size mismatch");
			auto& lh = (const tjlUtility::smallest_int_type<sizeof(V_T) <<
				UnZippedLogDim>::int_t&)lhs.internal_V_T_array;
			auto& rh = (const tjlUtility::smallest_int_type<sizeof(V_T) <<
				UnZippedLogDim>::int_t&)rhs.internal_V_T_array;
			return less.operator() <sizeof(V_T) << UnZippedLogDim> (lh, rh);
		}
	};
}
#endif // block_h__
