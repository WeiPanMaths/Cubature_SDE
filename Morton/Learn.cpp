////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
// Learn.cpp : Defines main code for the console application.
//

#include "stdafx.h"
#include "morton.h"

#include <vector>
#include <iostream>
#include <random>
#include <functional>

// the point type
const size_t dimn = 9;
typedef double point_t[dimn];
typedef morton::block<point_t>::BLOCK BLOCK;
typedef morton::block<point_t>::UNBLOCK UNBLOCK; //nb UNBLOCK(BLOCK(UNBLOCK(a))) == UNBLOCK(a) 

// a demonstration of the techniques
void Learn()
{
	std::cout <<
		"A dataset comprises points distributed in" << dimn <<" dimensional cube with distribution \n\t(U1,U2^1/(1.125),...,U"<<dimn<<"^1/2).\n Other datasets are obtained by refection within the cube.\n\n" ;

	std::mt19937 Engine; // Mersenne twister MT19937

	std::uniform_real_distribution<double> UniformDistribution(.0, 1.);
	auto UniformGenerator = std::bind(UniformDistribution, Engine);

	std::uniform_int_distribution<size_t> ZeroOneDistribution(0, 1);
	auto ZeroOneGenerator = std::bind(ZeroOneDistribution, Engine);
	auto randombool = [&ZeroOneGenerator]()->
		bool {return (ZeroOneGenerator() > 0) ? true : false;};

	auto CreateNewScenario = [&UniformGenerator] (point_t& currentpoint, double Discrimination, bool flavour)->void
	{
		double myindex = 1.;
		for (auto it1 = std::begin(currentpoint);
			it1 != std::end(currentpoint) ; myindex += Discrimination, ++it1)
			*it1 = (flavour) ? .5 + .5 * pow(UniformGenerator(),
			1 / myindex) : 1 - .5 * pow(UniformGenerator(), 1 / myindex);
	};

	size_t TrainingSetSize = 3000000;
	double Discrimination = 0.125;

	std::vector<UNBLOCK> TemporaryNativeDataStore(TrainingSetSize);
	std::cout <<
		"Generating the 3000000 " << dimn <<" dimensional points\n";
	for(auto it = begin(TemporaryNativeDataStore);
		it != end(TemporaryNativeDataStore) ; ++it)
		CreateNewScenario( *it, Discrimination, false);
	std::vector<BLOCK> BlockedTrainingSetOne(TrainingSetSize);
	std::cout << "Splicing the coordinates points\n";
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedTrainingSetOne), [](const UNBLOCK &arg)->
		BLOCK {return BLOCK(arg);});
	std::cout << "Done! \n\nGenerating the 3000000 " << dimn <<" dimensional points\n";

	for(auto it = begin(TemporaryNativeDataStore);
		it != end(TemporaryNativeDataStore) ; ++it)
		CreateNewScenario( *it, Discrimination, true);
	std::vector<BLOCK> BlockedTrainingSetTwo(TrainingSetSize);
	std::cout << "Splicing the coordinates points\n";
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedTrainingSetTwo), [](const UNBLOCK &arg)->
		BLOCK {return BLOCK(arg);});
	std::cout << "Done! \n\n";
	std::cout << "Sorting points into order\n";
	std::sort(std::begin(BlockedTrainingSetOne),
		std::end(BlockedTrainingSetOne));
	std::sort(std::begin(BlockedTrainingSetTwo),
		std::end(BlockedTrainingSetTwo));
	std::cout << "Both datasets sorted\n\n";

	auto no_in_box = [] (const size_t BoxDepth, const std::vector<BLOCK>& Data, const BLOCK& Point)->size_t
	{
		morton::block<point_t>::compare less(BoxDepth);
		auto First = std::begin(Data), Last = std::end(Data);
		return std::upper_bound(First, Last, Point,
			less) - std::lower_bound(First, Last, Point, less);
	};

	auto d_stats_test = [](double x, double y)->double
	{
		//http://webcache.googleusercontent.com/search?q=cache:Sf2pFHBGBUwJ:www.biology.ed.ac.uk/research/groups/jdeacon/statistics/tress10.html+Poisson+distribution+for+count+data+student+t&cd=1&hl=en&ct=clnk&gl=uk
		return abs(x - y)/sqrt(x+y);
	};

	point_t CurrentPoint;
	for (size_t C = 10; --C > 0; ) {
		bool Flavour = randombool();
		CreateNewScenario(CurrentPoint, Discrimination, Flavour);
		size_t x, y;
		for (size_t D = 52; D > 47; --D) {
			x = no_in_box(D,
				BlockedTrainingSetOne, CurrentPoint);
			y = no_in_box(D,
				BlockedTrainingSetTwo, CurrentPoint);
			std::cout << Flavour << " "	<< D << " " << x << " " << y <<
				" d="<< (((x+y)>0)? d_stats_test(double(x), double(y)) : -1) << "\n";
		}
		std::cout << "\n";
	}
}

