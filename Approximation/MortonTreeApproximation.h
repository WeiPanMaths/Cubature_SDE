#pragma once

#include <map>
#include <memory>
#include "Utils.h"
#include "MortonTree.h"

//template<typename F>
//struct Func
//{
//	F* pFunc;
//	double operator() (const point_2d & pt)
//	{
//		return (*pFunc)(pt);
//	}
//};

//struct TestValueType
//{
//	TestValueType(double v, double x, double y) : _value(v), _x(x), _y(y) {}
//	double _value;
//	double _x;
//	double _y;
//};

//template<typename Func>
class MortonTreeApproximation
{
	typedef MortonTreeApproximation	Base;
	typedef std::map<BLOCK, double> MapPointsValues;
	//typedef std::map<BLOCK, TestValueType> MapPointsValues;

	MapPointsValues pointsValues;
	MortonTree mortonTree;
	FunctionBase*	pFunc;
	double  err;
	size_t deg;

	// TEST SECTION
	size_t TEST_NUM_EXACT_CALLS;
	size_t TEST_NUM_APPROX_CALLS;
	size_t TEST_NUM_INTERP;
	size_t _test_count_pt;
	// END TEST SECTION
public:
	MortonTreeApproximation(FunctionBase*, size_t, double);
	~MortonTreeApproximation(void) {}

	double operator() (const point_2d&);
	template<typename D>
	double operator() (const point_2d& pt, D& data) { return (*this)(pt); }

	size_t getNoOfExactFunctionCalls() const;
	size_t getNoOfApproxFunctionCalls() const;
	void resetNoOfExactFunctionCalls();
	void resetNoOfApproxFunctionCalls();
	size_t getNoOfInterpolations() const;

	// TEST SECTION
	void TEST_INC_NUM_EXACT_CALLS();
	void TEST_INC_NUM_APPROX_CALLS();
	void TEST_INC_NUM_INTERP();
	void TEST_SHOW_INTERPOLATED_BOXES() const;
	size_t TEST_GET_NUM_OF_PATCHES() const;
	size_t TEST_GET_POINTSVALUES_SIZE() const;
	size_t TEST_GET_TESTCOUNTPT() const;
	void TEST_INC_TESTCOUNTPT();
	// END TEST SECTION

	friend class MortonTreeOccupant;
	friend class MortonTree;

private:
	double getFuncEvaAndStore(const BLOCK&, const point_2d&);
	MapPointsValues::const_iterator getLowerbound(const BLOCK&, COMPARE&) const;
	MapPointsValues::const_iterator getUpperbound(const BLOCK&, COMPARE&) const;
	bool tryInterpolation(MapPointsValues::const_iterator&, MapPointsValues::const_iterator&, MortonTreeOccupant*, const BLOCK& );
	
};	// end MortonTreeApproximation




//	MortonTreeApproximation(Func& f, size_t degree, double e) 
//	: pFunc(f), deg(degree), err(e), pointsValues(),mortonTree(),TEST_NUM_EXACT_CALLS(0), TEST_NUM_APPROX_CALLS(0), TEST_NUM_INTERP(0), _test_count_pt(0)
//{ 
//	mortonTree.pApprox = this; 
//}
//	double operator() (const point_2d& pt)
//{
//	TEST_INC_TESTCOUNTPT();
//	BLOCK pt_block(pt);
//	auto pLeaf = mortonTree.query(pt_block);
//	
//	if (pLeaf->getLocalApprox(pt_block))  // if true ie local approx exists then
//		return pLeaf->getInterpApprox(pt);
//	else {
//		auto value = getFuncEvaAndStore(pt_block, pt);
//		pLeaf->attemptApprox(pt_block);
//		return value;
//	};
//	//else return getFuncEvaAndStore(pt_block, pt);
//}
//	size_t getNoOfExactFunctionCalls() const { return TEST_NUM_EXACT_CALLS;}
//	size_t getNoOfApproxFunctionCalls() const { return TEST_NUM_APPROX_CALLS; }
//	size_t getNoOfInterpolations() const { return TEST_NUM_INTERP;}
//	void TEST_INC_NUM_EXACT_CALLS() { TEST_NUM_EXACT_CALLS++;}
//	void TEST_INC_NUM_APPROX_CALLS() { TEST_NUM_APPROX_CALLS++;}
//	void TEST_INC_NUM_INTERP() { TEST_NUM_INTERP++; }
//	size_t TEST_GET_NUM_OF_PATCHES() const { return mortonTree.size(); }
//	size_t TEST_GET_POINTSVALUES_SIZE() const { return pointsValues.size(); }
//	void TEST_SHOW_INTERPOLATED_BOXES() const
//	{
//		auto display = [&] (const MortonTree::value_type & it)->void
//		{
//			if(it.second->pInterpolation)
//			{
//				std::cout << "id: " << std::distance(mortonTree.begin(), mortonTree.find(it.first)) << ",";
//				// obtain first point in the box
//				auto lower_bound = getLowerbound(it.first,COMPARE(it.second->mDepth));
//				auto upper_bound = --getUpperbound(it.first,COMPARE(it.second->mDepth));
//				auto pod_lower = point_2d_POD(UNBLOCK(lower_bound->first));
//				auto pod_upper = point_2d_POD(UNBLOCK(upper_bound->first));
//				std::cout << "[" << pod_lower[0] << ", " << pod_lower[1] << "],  ";
//				std::cout << "[" << pod_upper[0] << ", " << pod_upper[1] << "],  ";
//				std::cout << "distance: " << std::distance(lower_bound, upper_bound) << "\n";
//			}
//		};
//		std::for_each(mortonTree.begin(), mortonTree.end(),display);
//		std::cout << "size of PointsValues: " << pointsValues.size() << "\n";
//	}
//	size_t TEST_GET_TESTCOUNTPT() const {return _test_count_pt; }
//	void TEST_INC_TESTCOUNTPT()  { _test_count_pt++; }
//
//private:
//	double getFuncEvaAndStore(const BLOCK& pt_block, const point_2d& pt)
//	{
//		//TEST
//		TEST_INC_NUM_EXACT_CALLS();
//		return pointsValues.insert(MapPointsValues::value_type(pt_block,pFunc(pt))).first->second;
//		/*auto pt1 = point_2d_POD(UNBLOCK(pt_block));
//		return pointsValues.insert(MapPointsValues::value_type(pt_block,TestValueType(pFunc(pt),pt1[0],pt1[1]))).first->second._value;*/
//	}
//	MapPointsValues::const_iterator getLowerbound(const BLOCK & pt, COMPARE & compare) const
//	{
//		auto iterator_compare_lowerbound = [&compare] (const MapPointsValues::value_type & val, const BLOCK & ref_point)->bool
//		{
//			return compare(val.first,ref_point);
//		};
//		return std::lower_bound(pointsValues.begin(), pointsValues.end(), pt, iterator_compare_lowerbound);
//	}
//	MapPointsValues::const_iterator getUpperbound(const BLOCK & pt, COMPARE & compare) const
//	{
//		auto iterator_compare_upperbound = [&compare] (const BLOCK & ref_point, const MapPointsValues::value_type & val)->bool
//		{
//			return compare(ref_point, val.first);
//		};
//
//		return std::upper_bound(pointsValues.begin(), pointsValues.end(), pt, iterator_compare_upperbound);
//	}
//	bool tryInterpolation(MapPointsValues::const_iterator& first, 
//		MapPointsValues::const_iterator& last, MortonTreeOccupant* pOc, const BLOCK& pt_block)
//	{
//		const size_t no_of_pts  = std::distance(first,last);
//		std::vector<point_2d>	all_pts(no_of_pts);
//		std::vector<double>		all_values(no_of_pts);
//		std::vector<point_2d>	interpolation_pts(no_of_pts / 2 );
//		std::vector<double>		interpolation_values(no_of_pts / 2 );
//		size_t vec_index = 0;
//		for(auto it=first; it!=last; ++it, ++vec_index)
//		{
//			auto unblocked_pt = UNBLOCK(it->first);
//			std::copy( std::begin(point_2d_POD(unblocked_pt)), std::end(point_2d_POD(unblocked_pt)), std::begin(all_pts[vec_index]));
//			//all_values[vec_index] = it->second._value;
//			all_values[vec_index] = it->second;
//		}
//		for(size_t i = 0; i < no_of_pts /2; ++i)
//		{
//			interpolation_pts[i] = all_pts[i *2];
//			interpolation_values[i] = all_values[i *2];
//		}
//		pOc->pInterpolation = std::make_shared<InterpolationMorton>(deg);
//		//assert(pOc->pInterpolation->interpolate(interpolation_pts, interpolation_values)==true);					// THIS BREAKS IN RELEASE MODE FOR SOME REASON
//		//if(pOc->pInterpolation->isAccurate(all_pts, all_values, err))													
//		if(pOc->pInterpolation->interpolate(interpolation_pts, interpolation_values) && pOc->pInterpolation->isAccurate(all_pts, all_values, err))
//		{
//			TEST_INC_NUM_INTERP();	//TEST
//			//SHOW(pOc->pInterpolation);
//			return true;
//		} else {
//			pOc->pInterpolation.reset();
//			//SHOW(this->_test_count_pt);	// SHOW WHERE WE ARE
//			//auto pod_lower = point_2d_POD(UNBLOCK(first->first));
//			//auto pod_upper = point_2d_POD(UNBLOCK((--last)->first));
//			//auto pod	   = point_2d_POD(UNBLOCK(pt_block));
//			//std::cout << "MT:" << this->TEST_GET_TESTCOUNTPT() << " ("<< pod[0]   << ", " << pod[1] << "), ";
//			//std::cout << "[" << pod_lower[0] << ", " << pod_lower[1] << "], "; //"], value " << first->second << ", ";
//			//std::cout << "[" << pod_upper[0] << ", " << pod_upper[1] << "], "; //"],  value " << last->second  << ", ";
//			//auto pod_next_lower = point_2d_POD(UNBLOCK((++last)->first));
//			////auto pod_next_lower = point_2d_POD(UNBLOCK((last)->first));
//			////std::cout << "[" << pod_next_lower[0] << ", " << pod_next_lower[1] << "],  value " << last->second  << ", ";
//			//std::cout << std::distance(first, last) << " ";
//		 
//			mortonTree.subdivideNode(pOc->pTreeNode);
//			return false;
//		}
//	}
//
//