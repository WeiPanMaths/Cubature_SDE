#include "stdafx.h"
#include "morton.h"

#include <vector>
#include <iostream>
#include <random>
#include <functional>

typedef double point_3d[3];
typedef morton::block<point_3d>::BLOCK BLOCK3D;
typedef morton::block<point_3d>::UNBLOCK UNBLOCK3D; //nb UNBLOCK(BLOCK(UNBLOCK(a))) == UNBLOCK(a) 

typedef double point_2d[2];
typedef morton::block<point_2d>::BLOCK BLOCK2D;
typedef morton::block<point_2d>::UNBLOCK UNBLOCK2D;

typedef double point_1d[1];
typedef morton::block<point_1d>::BLOCK BLOCK1D;
typedef morton::block<point_1d>::UNBLOCK UNBLOCK1D;


void Test3DBox()
{	
	size_t	SetSize = 6;
	std::vector<UNBLOCK3D>  TemporaryNativeDataStore(SetSize);
	
	point_3d myData[6] = {	{0.51, 0.51, 0.51}, 
	{0.575, 0.575, 0.575}, { 0.525, 0.525, 0.525 }, {0.59, 0.59, 0.59} , {0.525, 0.59, 0.525},  {0.59, 0.525, 0.51} };

	auto CreateScenario  = [&myData] (point_3d & pt, size_t position)->void
	{
		size_t index = 0;
		for ( auto it = std::begin(pt); it!=std::end(pt); ++it, ++index)
			*it = myData[position][index];
	};

	size_t index = 0;
	for(auto it = std::begin(TemporaryNativeDataStore); it!=std::end(TemporaryNativeDataStore); ++it, ++index)
		CreateScenario(*it, index);

	std::vector<BLOCK3D> BlockedPointSet(SetSize);
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedPointSet), [](const UNBLOCK3D &arg)->
		BLOCK3D {return BLOCK3D(arg);});

	std::sort(std::begin(BlockedPointSet), std::end(BlockedPointSet));

	std::vector<UNBLOCK3D> UnblockedPointSet(SetSize);
	for(size_t i =0; i < BlockedPointSet.size(); ++i)
		UnblockedPointSet[i] = UNBLOCK3D(BlockedPointSet[i]);

	auto no_in_box = [] (const size_t BoxDepth, const std::vector<BLOCK3D>& Data, const BLOCK3D& Point)->size_t
	{
		morton::block<point_3d>::compare less(BoxDepth);
		auto First = std::begin(Data), Last = std::end(Data);
		
		auto up = std::upper_bound(First, Last, Point, less);
		auto low =std::lower_bound(First, Last, Point, less);

		UNBLOCK3D up_unblock ( *(--up));
		UNBLOCK3D low_unblock( *low);
		
		return	std::upper_bound(First, Last, Point, less) 
				- std::lower_bound(First, Last, Point, less);
	};

	
	point_3d current_pt = {0.52, 0.52, 0.52 };

	size_t depth = 49;
	std::cout << "depth: " << depth << ". number in box: " << no_in_box(depth, BlockedPointSet, current_pt) << std::endl;
}


void TestSubdivsion()
{
	size_t no_of_pts(10000);
	std::mt19937 Engine; // Mersenne twister MT19937
	std::uniform_real_distribution<double> UniformDistribution(0.5, 1.);

	auto LogSpotGenerator = std::bind(UniformDistribution, Engine);
	auto CreateNewPoint = [&LogSpotGenerator](point_2d& currentpoint)->void
	{
		currentpoint[0] = LogSpotGenerator();
		currentpoint[1] = LogSpotGenerator();
	};

	std::vector<point_2d> interpolation_pts(no_of_pts);
	for (auto it = std::begin(interpolation_pts); it != std::end(interpolation_pts); ++it)
		CreateNewPoint(*it);

	std::cout << "done" << std::endl;
}


void Test2DBox(const size_t depth, const point_2d ref_pt)
{	
	std::cout << "Test2DBox()" << std::endl;

	size_t	SetSize = 4;
	std::vector<UNBLOCK2D>  TemporaryNativeDataStore(SetSize);
	
	//point_2d myData[4] = {	{ -0.21, -0.31 }, { -0.15, -0.85 },	{ -0.64, -0.58 }, { -0.88, -0.34 } };
	point_2d myData[4] = {	{ 0.85, 0.85 },	{ 0.55, 0.75 }, { 0.94, 0.52 }, { 0.51, 0.53 } };
	//point_2d myData[8] = {	{ 0.5, 0.5 }, { 0.21, 0.31 }, { 0.64, 0.58 }, { -0.5, -0.5 },	{ -0.5, 0.5 },  { 0.5, -0.5 }, { 0.15, 0.85 },	{ 0.88, 0.34 } };

	auto CreateScenario = [&myData](point_2d & pt, size_t position)->void
	{
		size_t index = 0;
		for (auto it = std::begin(pt); it != std::end(pt); ++it, ++index) {
			*it = myData[position][index];
			//std::cout << myData[position][index] << " ";
		}
		//std::cout << std::endl << std::endl;
	};

	// write myData into TemporaryNativeDataStore
	size_t index = 0;
	for (auto it = std::begin(TemporaryNativeDataStore); it != std::end(TemporaryNativeDataStore); ++it, ++index)
		CreateScenario(*it, index);
	
	std::vector<BLOCK2D> BlockedPointSet(SetSize);
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedPointSet), [](const UNBLOCK2D &arg)->
		BLOCK2D {return BLOCK2D(arg);});

	std::sort(std::begin(BlockedPointSet), std::end(BlockedPointSet));

	std::cout << "in order ..." << std::endl;
	std::vector<UNBLOCK2D> UnblockedPointSetInMortonOrder(SetSize);
	for (size_t i = 0; i < BlockedPointSet.size(); ++i)
	{
		UnblockedPointSetInMortonOrder[i] = UNBLOCK2D(BlockedPointSet[i]);
		point_2d _pt;
		std::copy(std::begin((point_2d &)UnblockedPointSetInMortonOrder[i]), std::end((point_2d &)UnblockedPointSetInMortonOrder[i]), std::begin(_pt));
		std::cout << _pt[0] << " " << _pt[1] << std::endl;
	}
	std::cout << std::endl;

	auto no_in_box = [] (const size_t BoxDepth, const std::vector<BLOCK2D>& Data, const BLOCK2D& Point)->size_t
	{
		morton::block<point_2d>::compare less(BoxDepth);
		auto First = std::begin(Data), Last = std::end(Data);
		
		auto up = std::upper_bound(First, Last, Point, less);
		auto low =std::lower_bound(First, Last, Point, less);

		UNBLOCK2D up_unblock ( *(--up));
		UNBLOCK2D low_unblock( *low);

		point_2d up_pt, low_pt;
		std::copy(std::begin((point_2d &)up_unblock), std::end((point_2d &)up_unblock), std::begin(up_pt));
		std::copy(std::begin((point_2d &)low_unblock), std::end((point_2d &)low_unblock), std::begin(low_pt));

		std::cout << up_pt[0] << " " << up_pt[1] << std::endl;
		std::cout << low_pt[0] << " " << low_pt[1] << std::endl;

		return	std::upper_bound(First, Last, Point, less) - std::lower_bound(First, Last, Point, less);
	};
	
	point_2d current_pt = {ref_pt[0], ref_pt[1]};

	//size_t depth = 63;
	std::cout << "run .." << std::endl;
	std::cout << "depth: " << depth << ". number in box: " << no_in_box(depth, BlockedPointSet, current_pt) << std::endl;
}


// test biggest box that contains two points
void TestBiggestBox()
{
	size_t SetSize = 2;
	double pt1[1] = { 0.511234};
	double pt2[1] = { 2};
	std::vector<UNBLOCK1D> TemporaryNativeDataStore(SetSize);

	TemporaryNativeDataStore[0] = UNBLOCK1D(pt1);
	TemporaryNativeDataStore[1] = UNBLOCK1D(pt2);

	std::vector<BLOCK1D> BlockedPointSet(SetSize);
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedPointSet), [](const UNBLOCK1D &arg)->
		BLOCK1D {return BLOCK1D(arg);});

	size_t BoxDepth = 63;
	for( ; BoxDepth != -1; --BoxDepth)
	{
		morton::block<point_1d>::compare less(BoxDepth);
		if(less(BlockedPointSet[0], BlockedPointSet[1]))
			std::cout << "depth " << BoxDepth << ". not in the same box " << std::endl;
		else
			std::cout << "depth " << BoxDepth << ". in the same box     " << std::endl;
	}
}

// a demonstration of the techniques
void TestBox()
{
	//std::cout << "Test box with points 0.1, 0.125, 0.25, 0.375, 0.5, 0.75, 1.26, 1.33, 1.5, 2.00, 2.343, 3.0" << std::endl;
	
	//std::cout << "Test box with points 0.1, 0.4, 0.75, 1.0" << std::endl;
	const size_t	SetSize = 12;
	
	std::vector<UNBLOCK1D>  TemporaryNativeDataStore(SetSize);
	double myData[SetSize] = {  0.5, 0.75, 0.1, 0.125, 0.25, 0.375,  1.26, 1.33, 1.5, 2.00, 2.343, 3.0 };


	// lambda function definition
	// turns original POD into UNBLOCK
	auto CreateScenario  = [&myData] ( point_1d & pt, size_t position ) -> void
	{
		for ( auto it = std::begin(pt); it!=std::end(pt); ++it)
			*it = myData[position];
	};
	
	// this calls the previous lambda definition and turns myData into UNBLOCK format
	size_t index = 0;
	for(auto it = std::begin(TemporaryNativeDataStore); it!=std::end(TemporaryNativeDataStore); ++it, ++index)
		CreateScenario(*it, index);



	std::vector<BLOCK1D> BlockedPointSet(SetSize);
	std::transform(std::begin(TemporaryNativeDataStore), std::end(TemporaryNativeDataStore), begin(BlockedPointSet), 
		[](const UNBLOCK1D &arg)->	BLOCK1D {return BLOCK1D(arg);}		);

	
	std::sort(std::begin(BlockedPointSet), std::end(BlockedPointSet));


	auto no_in_box = [] (const size_t BoxDepth, const std::vector<BLOCK1D>& Data, const BLOCK1D& Point)->size_t
	{
		morton::block<point_1d>::compare less(BoxDepth);

		auto First = std::begin(Data), Last = std::end(Data);
		
		auto up  = std::upper_bound(First, Last, Point, less);
		auto low = std::lower_bound(First, Last, Point, less);
		
		double * low_pt= UNBLOCK1D(*low);
		double * up_pt = UNBLOCK1D(*up);

		std::cout << std::endl << "low_pt = " << low_pt[0] << std::endl;
		std::cout << "up_pt  = " << *up_pt << std::endl ;

		return	std::upper_bound(First, Last, Point, less) 
				- std::lower_bound(First, Last, Point, less);
	};

	point_1d current_pt = { 0.4 };

	for(size_t depth = 52; depth > 47; --depth)
		std::cout << "depth: " << depth << ". number in box: " << no_in_box(depth, BlockedPointSet, current_pt) << std::endl;
}

//int main()
//{
//	TestBox();
//
//	return 0;
//}
