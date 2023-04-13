////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////

#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include "morton.h"
#include <vector>
//
//Eg. for any ab…d in {0,1} the 6 points in 5d described by the columns:
//.1ab…d0               .1ab…d1               .1ab…d1               .1ab…d1               .1ab…d1               .1ab…d1
//	.1ab…d0               .1ab…d0               .1ab…d1               .1ab…d1               .1ab…d1               .1ab…d1
//	.1ab…d0               .1ab…d0               .1ab…d0               .1ab…d1               .1ab…d1               .1ab…d1
//	.1ab…d0               .1ab…d0               .1ab…d0               .1ab…d0               .1ab…d1               .1ab…d1               
//	.1ab…d0               .1ab…d0               .1ab…d0               .1ab…d0               .1ab…d0               .1ab…d1
//
//	should always be in the same order and changing the ab…d or adding extra digits after the 0 and 1 should not affect this. (It would be worth understanding the order that does emerge!)
//

const size_t dimn = 5;
typedef double point_t[dimn];
typedef morton::block<point_t>::BLOCK BLOCK;
typedef morton::block<point_t>::UNBLOCK UNBLOCK;


void test()
{
	const size_t max = 1<<5;
	std::vector<bool> answers;
	auto initfn = [](unsigned char a,unsigned char b, unsigned char c)->
		size_t
	{
		return ((((((1 << 1) + a)<<1) + b)<<1)+c)<<1;
	};
	
	std::vector<std::vector<double>> ds;
	for (size_t count = 0; count < 8 ; ++count)
	{
		size_t init = initfn(count%2,(count%4)/2,(count%8)/4);
		ds.clear();
		auto no_patterns = (1<<dimn);
		for (size_t bitcycler = 0; bitcycler < no_patterns; ++bitcycler)
		{
			std::vector<size_t> temp;
			size_t bits = bitcycler;
			for (size_t count = 0 ; count < dimn ; ++count){
				temp.push_back(init + bits % 2);
				bits/=2;
			}
			std::vector<double> tempd(temp.size());
			std::transform(begin(temp), end(temp), begin(tempd), [&max](const size_t arg)->double {return double(arg)/max;});
			ds.push_back(tempd);
		}

		if (answers.size()==0) for (auto it = begin(ds); it != end(ds); ++it) for (auto it1 = begin(ds); it1 != end(ds); ++it1)
			answers.push_back(BLOCK((point_t &)((*it)[0])) < BLOCK((point_t&)((*it1)[0])));
		std::vector<bool> temp_answers;
		temp_answers.clear();
		for (auto it = begin(ds); it != end(ds); ++it) for (auto it1 = begin(ds); it1 != end(ds); ++it1)
			temp_answers.push_back(BLOCK((point_t &)((*it)[0])) < BLOCK((point_t&)((*it1)[0])));
		if(temp_answers!=answers) std::cout << "Error with value " << init << "\n";
		std::cout << "test all done\n";
	}
	std::cout << "\n";
}
