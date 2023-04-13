////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef smallest_int_type_h__
#define smallest_int_type_h__
#include <limits>
#include <xutility>


namespace tjlUtility {
	// smallest_int_type<W>::int_t 
	// provides a template POD type for dealing with a block of binary data of size W (which must be a power of two) 
	// it also provides limited methods and templated functions objects (which bind the number of bits) for
	// comparison, ignoring low order bits, and zeroing of the low order bits for any arguments of type smallest_int_type<W>::int_t
	template <size_t W>
	struct
		smallest_int_type;

	// bind "bits" in smallest_int_type<W>::round_down_sit(lhs, bits)
	struct round_down_block_tjlU{
		size_t bits;
		round_down_block_tjlU(const size_t& Bits) : bits(Bits)
		{
		};
		template <size_t W>
		typename smallest_int_type<W>::int_t operator()(const typename smallest_int_type<W>::int_t& lhs)
		{
			return smallest_int_type<W>::round_down_sit(lhs, bits);
		}
	};

	// bind "bits" in smallest_int_type<W>::compare_less_sit(lhs, rhs, bits)
	struct compare_less_tjlU {
		size_t bits;
		compare_less_tjlU(const size_t& Bits) : bits(Bits)
		{
		};

		template <size_t W>
		bool operator()(const typename smallest_int_type<W>::int_t& lhs,
			const typename smallest_int_type<W>::int_t& rhs)
		{
			return smallest_int_type<W>::compare_less_sit(lhs, rhs, bits);
		}
	};

	// smallest_int_type<W>::int_t is a POD integer type structure whose size is the power of two that equals W
	// smallest_int_type<W>::INTTYPE is the integer type used to build the structure int_t and will be equal to it if there is a standard int type that is large enough
	// static bool smallest_int_type<W>::compare_less_sit(const int_t& lhs, const int_t& rhs, size_t bits) compares lhs and rhs of type int_t ignoring the bottom bits bits
	// static int_t smallest_int_type<W>::round_down_sit(const int_t& lhs, size_t bits) returns an object of type int_t equating to lhs but with the bottom bits bits zeroed

	template <size_t W>
	struct
		smallest_int_type {
			static_assert(W == (1 << log_base2<W>::ans), "Power of 2 expected");
		typedef unsigned long long INTTYPE;
		enum {
			ARRAYDIM = ( W >> log_base2<sizeof(INTTYPE)>::ans)
		};
		union int_t {
			INTTYPE m_array[ARRAYDIM];
		} ;

		static int_t round_down_sit(const int_t& lhs, size_t bits)
		{
			size_t bits_to_ignore = bits%64; // bits and bits to be ignored
			size_t elements_to_ignore = bits/64;
			int_t ans = {0};
			std::copy(std::begin(lhs.m_array) + elements_to_ignore + 1, std::end(lhs.m_array), std::begin(ans.m_array) + elements_to_ignore + 1);
			(*(std::begin(ans.m_array) + elements_to_ignore)) = ((*(std::begin(lhs.m_array) + elements_to_ignore)) >>
				bits_to_ignore) << bits_to_ignore;
			return ans;
		}

		static bool compare_less_sit(const int_t& lhs, const int_t& rhs,
			size_t bits)
		{

			size_t bits_to_ignore = bits%64; // bits and bits to be ignored
			size_t elements_to_ignore = bits/64;

			bool temp1 = std::lexicographical_compare(lhs.m_array +
				elements_to_ignore + 1, std::end(lhs.m_array),
				rhs.m_array + elements_to_ignore + 1, std::end(rhs.m_array));
			return (temp1)
				?true
				:
				((((*(lhs.m_array + elements_to_ignore)) >>
				bits_to_ignore) < ((*(rhs.m_array + elements_to_ignore)) >>
				bits_to_ignore)) &&
				(std::equal(lhs.m_array + elements_to_ignore + 1,
				std::end(lhs.m_array), rhs.m_array+ elements_to_ignore + 1)
				))? true
				:false;
		}
	};

	template <>
	struct
		smallest_int_type<8> {
		static_assert(sizeof(unsigned long long)==8,"unsigned long long size assumptions invalid");
		typedef unsigned long long INTTYPE;
		typedef unsigned long long int_t;
		static int_t round_down_sit(const int_t& lhs, size_t bits)
		{
			return (lhs >> bits) << bits;
		}
		static bool compare_less_sit(const int_t& lhs, const int_t& rhs,
			size_t bits)
		{
			return (lhs >> bits) < (rhs >> bits);
		}
	};

	template <>
	struct
		smallest_int_type<4> {
		static_assert(sizeof(unsigned long)==4,"unsigned long size assumptions invalid");
		typedef unsigned long INTTYPE;
		typedef unsigned long int_t;
		static int_t round_down_sit(const int_t& lhs, size_t bits)
		{
			return (lhs >> bits) << bits;
		}
		static bool compare_less_sit(const int_t& lhs, const int_t& rhs,
			size_t bits)
		{
			return (lhs >> bits) < (rhs >> bits);
		}
	};

	template <>
	struct
		smallest_int_type<2> {
		static_assert(sizeof(unsigned short)==2,"unsigned short size assumptions invalid");
		typedef unsigned short INTTYPE;
		typedef unsigned short int_t;
		static int_t round_down_sit(const int_t& lhs, size_t bits)
		{
			return (lhs >> bits) << bits;
		}
		static bool compare_less_sit(const int_t& lhs, const int_t& rhs,
			size_t bits)
		{
			return (lhs >> bits) < (rhs >> bits);
		}
	};

	template <>
	struct
		smallest_int_type<1> {
		typedef unsigned char INTTYPE;
		typedef unsigned char int_t;
		static int_t round_down_sit(const int_t& lhs, size_t bits)
		{
			return (lhs >> bits) << bits;
		}
		static bool compare_less_sit(const int_t& lhs, const int_t& rhs,
			size_t bits)
		{
			return (lhs >> bits) < (rhs >> bits);
		}
	};

}

#endif // smallest_int_type_h__
