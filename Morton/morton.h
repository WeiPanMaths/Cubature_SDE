////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef morton_h__
#define morton_h__

#include "log_base2.h"
#include "smallest_int_type.h"
#include "block.h"
#include "unblock.h"
#include <utility>
#include <algorithm>
#include <iterator>

//Copyright January 1st 2011 Terry Lyons

namespace morton {

	using namespace tjlUtility;
	using namespace tjlblockdetail;

	// block<point> holds no data; it defines in a consistent way four public objects
	// block<point>::BLOCK a union that holds multidimensional data in a binary form useful for localising and categorising
	// block<point>::UNBLOCK a union that holds data in conventional form and has a constructor that accepts a BLOCK so can be used to roundtrip points
	// block<point>::compare provides a binary function object for comparing BLOCK objects at a given resolution
	// block<point>::round_down provides a unary function object for coarsening BLOCK objects to a given resolution
	template <class POINT>
	struct block;

	// the only implementation of block provided requires that point is an array type V[D] where V is POD and D is size_t
	// this implementation introduces the types P_T for the array V[D] and V_T for the POD type V and an integer constant D of type size_t for the dimension
	template <class V, size_t D>
	struct block<V[D]> {
		typedef V P_T[D];
		typedef V V_T;
		static const size_t DIM = D;
	private:
		friend union BLOCK;
		friend union UNBLOCK;
		static const unsigned LOGDIMPLUS = log_base2<DIM>::ans + ((DIM > (size_t(1) << log_base2<DIM>::ans))? 1 : 0);
		static const unsigned LOGVWIDTH =log_base2<sizeof(V_T)>::ans;
		static const size_t VWIDTHPADDED = (1 << LOGVWIDTH);
		typedef typename smallest_int_type<VWIDTHPADDED>::int_t INTVPADDED;
		static_assert((sizeof(V_T) <= sizeof(INTVPADDED)) && (sizeof(INTVPADDED) < (sizeof(V_T) << 1)), "INTVTPADDED is not the smallest integer cover of V_T");
	    static_assert((DIM >= (size_t(1) << log_base2<DIM>::ans)) && (DIM <= (size_t(1) << LOGDIMPLUS)), "LOGDIMPLUS not capturing length");
	public:

		union BLOCK
		{
			friend union UNBLOCK;
			friend struct compare;
			friend struct round_down;
		public:
			// default constructor
			BLOCK() : pmBlock()
			{
				memset(this,0,sizeof(BLOCK));
			}
			// copy constructor
			BLOCK(const BLOCK& arg) : pmBlock (arg.pmBlock)
			{
			}
			// move constructor
			BLOCK(BLOCK&& arg) : pmBlock (std::move(arg.pmBlock))
			{
			}
			// inversion constructor from UNBLOCK
			BLOCK(const UNBLOCK& arg) : pmUnBlock (arg.pmUnBlock)
			{
				pmBlock.flip();
			}

			template <class POINTLIKE>
			BLOCK(const POINTLIKE& arg): pmBlock() //essential if there is padding: including pmBlock() in the initializer list zeros it (and any padding) by the standard 
			{
				memset(this,0,sizeof(BLOCK));
				std::copy_n(std::begin(arg), DIM, external_data);
				pmBlock.flip();
			}

			bool operator <(const BLOCK& rhs) const
			{
				return pmBlock < rhs.pmBlock;
			}
			bool operator ==(const BLOCK& rhs) const
			{
				return pmBlock == rhs.pmBlock;
			}
		private:
			UnBlock<INTVPADDED, LOGDIMPLUS, 0> pmUnBlock;
			Block<INTVPADDED, 0, LOGDIMPLUS> pmBlock;
			P_T external_data;
		};


		// a function object providing a binary comparison operator that ignores the last few bits when making the comparison
		struct compare {
		private:
			tjlblockdetail::compare_less_block less;
		public:
			// the constructor sets the number of lower order bits to ignore
			compare(const size_t& bits) : less(bits << LOGDIMPLUS)
			{
			}

			bool operator()(const BLOCK& lhs, const BLOCK& rhs)
			{
				return less(lhs.pmBlock, rhs.pmBlock);
			}
		};

		//round_down_block
		// a function object providing a rounding operator that zeros the last few bits
		struct round_down {
		private:
			tjlblockdetail::round_down_block round;
		public:
			// the constructor sets the number of lower order bits to ignore
			round_down(const size_t& bits) : round(bits << LOGDIMPLUS)
			{
			}

			BLOCK operator()(const BLOCK& lhs)
			{
				BLOCK ans;
				ans.pmBlock = round(lhs.pmBlock);
				return ans;
			}
		};

		union UNBLOCK
		{
			friend union BLOCK;
			friend bool operator == (const UNBLOCK& lhs, const UNBLOCK& rhs) {
				return std::equal( std::begin((const P_T&)lhs), std::end((const P_T&)lhs),std::begin((const P_T&)rhs));
			}
		public:
			// default constructor
			UNBLOCK() : pmUnBlock()
			{
				memset(this,0,sizeof(BLOCK));
			}
			// copy constructor
			UNBLOCK(const UNBLOCK& arg) : pmUnBlock (arg.pmUnBlock)
			{
			}
			// move constructor
			UNBLOCK(UNBLOCK&& arg) : pmUnBlock (std::move(arg.pmUnBlock))
			{
			}
			// inversion constructor from BLOCK
			UNBLOCK(const BLOCK& arg) : pmUnBlock(arg.pmUnBlock)
			{
				pmUnBlock.flop();
			}
			// public conversion operator giving access to the coordinate values
			operator P_T&()
			{
				return external_data;
			}

			operator const P_T&()
			{
				return external_data;
			}

			//bool operator == (const BLOCK& rhs) const {return pmBlock == rhs.pmBlock;}
		private:
			UnBlock<INTVPADDED, LOGDIMPLUS, 0> pmUnBlock;
			Block<INTVPADDED, LOGDIMPLUS, 0> pmBlock;
			P_T external_data;
		};
	};
}

#endif // morton_h__
