////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
#ifndef log_base2_h__
#define log_base2_h__

namespace tjlUtility {
	template<unsigned x>
	struct
		log_base2 {
		enum {
			ans = 1 + log_base2<(x >> 1)>::ans
		};
	};

	template<>
	struct
		log_base2<1> {
		enum {
			ans = 0
		};
	};

	template<>
	struct
		log_base2<0> {
		//log_base2(){static_assert(false, "log can only be applied to strictly positive numbers");}
		enum {
			ans = 0
		};
	};
}

#endif // log_base2_h__
