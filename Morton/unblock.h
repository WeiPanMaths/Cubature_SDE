////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef unblock_h__
#define unblock_h__

#include <utility>
#include "lsplice.h"

namespace tjlblockdetail {

	template<class V_T, unsigned ZippedLogDim, unsigned UnZippedLogDim>
	union UnBlock;

	template<class V_T>
	union UnBlock<V_T, 0, 0>
	{
		typedef V_T currently_zipped_as_1dim_t[1];
		currently_zipped_as_1dim_t data[1];

		// The format for the transformed block we will unzip from
		V_T reblocked;

		// the code that splices pairs of spliced data.... until there is only one data block
		V_T& flop()
		{
		return reblocked;
		}
	};

	template<class V_T, unsigned ZippedLogDim>
	union UnBlock<V_T, ZippedLogDim, 0>
		  {

			// data as a 2 dimensional array 2^ZippedLogDim by 2^0
			typedef V_T currently_zipped_as_1dim_t[1<<ZippedLogDim];
			currently_zipped_as_1dim_t data[1];

			// The format for the transformed block we will unzip from
			UnBlock<V_T, ZippedLogDim - 1, 1> reblocked;

			// the code that splices pairs of spliced data.... until there is only one data block
		UnBlock<V_T, 0, ZippedLogDim>& flop()
			{
				{
				auto itNextIn = std::begin(data);
				auto itLastIn = std::end(data);
				auto itNextOut = std::begin(data);
					// working space for unzipped output
					currently_zipped_as_1dim_t temp;
					for (; itNextIn != itLastIn; ++itNextIn, ++itNextOut) {
						lusplice(temp, *itNextIn);
						// in and out point to the same space
					auto itOut = std::begin(*itNextOut);
						std::copy(std::begin(temp), std::end(temp), itOut);
					}
				}
				return reblocked.flop();
			}
		  };

	template<class V_T, unsigned ZippedLogDim, unsigned UnZippedLogDim>
	union UnBlock
		  {

			typedef V_T currently_zipped_as_1dim_t[1<<ZippedLogDim];
			currently_zipped_as_1dim_t data[1<<UnZippedLogDim];

			// The format for the transformed block we will unzip from
			UnBlock<V_T, ZippedLogDim - 1, UnZippedLogDim + 1> reblocked;

		UnBlock<V_T, 0, ZippedLogDim + UnZippedLogDim>& flop()
			{
				{
				auto itNextIn = std::begin(data);
				auto itLastIn = std::end(data);
				auto itNextOut = std::begin(data);
					currently_zipped_as_1dim_t temp;
					for (; itNextIn != itLastIn; ++itNextIn, ++itNextOut) {
						lusplice(temp, *itNextIn);
						// in and out point to the same space
						std::copy(std::begin(temp),
							std::end(temp), (*itNextOut));
					}
				}
				return reblocked.flop();
			}
		  };

	template<class V_T, unsigned UnZippedLogDim>
	union
		UnBlock<V_T, 0, UnZippedLogDim> {
		// flop completed
		UnBlock& flop()
		{
			return *this;
		}
	};
}
#endif // unblock_h__
