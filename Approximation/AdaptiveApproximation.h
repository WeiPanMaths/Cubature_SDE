#pragma once

#include <map>
#include <set>
//#include <unordered_set>
//#include <functional>
#include <algorithm>
//#include <cmath>
#include <limits>
#include "Utils.h"

#if dim_ == 2
#include "LeastSquaresInterpolation.h" // 2D?
#endif
#if dim_ == 1
#include "LocalApproximation.h"
#endif

// to access interpolation use BLOCK value and comparison operator resolution
struct ApproxMortonValue
{
	size_t						mCompareDepth;					// comparison operator resolution; for box determination
	double						mExactFunctionEvaluation;
	BLOCK						mInterpolationKey;
	bool						mIsInterpolated;

	ApproxMortonValue(double fneva = 0, size_t depth = DEFAULTDEPTH, bool isinterpolated = false)
		: mCompareDepth(depth)
		, mExactFunctionEvaluation(fneva)
		, mIsInterpolated(isinterpolated) {}
};

#if dim_ == 2
template <typename Function, typename Datatype>
class ApproximationMorton : private std::map< BLOCK, ApproxMortonValue >
{
	typedef		std::map< BLOCK, InterpolationMorton>	MapOfInterpolations;
	typedef		std::map< BLOCK, ApproxMortonValue >	Base;
	typedef		Base::value_type						BaseValueType;
	typedef		Base::iterator							BaseIterator;
	typedef		MapOfInterpolations::value_type			MapOfInterpolationsValueType;
	
	using std::map< BLOCK, ApproxMortonValue >::lower_bound;
	using std::map< BLOCK, ApproxMortonValue >::upper_bound;

	double					myError;						// absolute or relative?
	Function &				myFunction;
	unsigned int			myDegree;
	MapOfInterpolations		myMapOfInterpolations;
	size_t					myMaxPtsPerBox;

	//BaseIterator	mylit;	// used to track box lowerbound upperbound search range
	//BaseIterator	myrit;	// since no box should contain more than myMaxPtsPerBox 
	BaseIterator	mysavedit;

	//////////////
	///	TESTING
	size_t					myNoOfApproxFunctionCalls;
	size_t					myNoOfExactFunctionCalls;
	size_t					myPTCOUNT;
	size_t					myHitBottom;
	/////////////
public:
	
	ApproximationMorton( Function & function, unsigned int degree, double abserror ) 
		: myError(abserror)
		, myFunction(function)
		, myDegree(degree)
#if dim_ == 2
		, myMaxPtsPerBox( (2 + degree)*(1 + degree)/2 * 3 ) // 3 times num of monomials for given degree parameter
#endif
#if dim_ == 1
		, myMaxPtsPerBox(3*(1+degree)) // degree 1
#endif
		//, myrit(end())
		//, mylit(begin())
		, myMapOfInterpolations()
		, myNoOfApproxFunctionCalls(0)
		, myNoOfExactFunctionCalls(0), myPTCOUNT(0), myHitBottom(0)
	{
		//std::cout << DEFAULTDEPTH << " default depth \n";
	}

	~ApproximationMorton(void) {}
	
	size_t  getSizeOfPointCloud() const
	{
		return size();
	}

	size_t  getNoOfInterpolations() const
	{
		return myMapOfInterpolations.size();
	}

	void getInterpolatedPatchesCoeffs(std::vector<double> &patches) const
		/// WARNING: assumes positive domain! and normalised box [0.5, 1)
	{
		for (MapOfInterpolations::const_iterator it = myMapOfInterpolations.begin();
			it != myMapOfInterpolations.end();
			++it)
		{
			BLOCK patch_rep(it->first);  // key
			size_t patch_depth;
			try
			{
				auto pitr = find(patch_rep);

				if (pitr != end())
					patch_depth = pitr->second.mCompareDepth;
				else
					throw myexception();
			}
			catch (myexception& e)
			{
				std::cout << "key not found " << e.what() << std::endl;
			}

			auto unblocked_pt = UNBLOCK(patch_rep);
			point_2d pt;
			std::copy(std::begin((point_2d_POD &)unblocked_pt), std::end((point_2d_POD &)unblocked_pt), std::begin(pt));

			double interval_length = getIntervalLength_(patch_depth);

			//std::cout << "rep " << pt[0] << ", " << pt[1] << std::endl;

			double x_interval_id(std::floor((pt[0] - 0.5) / interval_length));
			double y_interval_id(std::floor((pt[1] - 0.5) / interval_length));

			//std::cout << 0.5 + x_interval_id * interval_length << ", " << 0.5 + (x_interval_id + 1) * interval_length << std::endl;
			//std::cout << 0.5 + y_interval_id * interval_length << ", " << 0.5 + (y_interval_id + 1) * interval_length << std::endl;
			//std::cout << "\n";
			patches.push_back(0.5 + x_interval_id * interval_length);
			patches.push_back(0.5 + (x_interval_id + 1) * interval_length);
			patches.push_back(0.5 + y_interval_id * interval_length);
			patches.push_back(0.5 + (y_interval_id + 1) * interval_length);

			auto coeffs = it->second.getCoeffs();
			for (auto itr = coeffs.begin(); itr != coeffs.end(); ++itr)
				patches.push_back(*itr);
		}
	}

	//void getInterpolatedPatches(std::vector<double> &patches) const
	///// WARNING: assumes positive domain! and normalised box [0.5, 1)
	//{
	//	for (MapOfInterpolations::const_iterator it = myMapOfInterpolations.begin();
	//		it != myMapOfInterpolations.end();
	//		++it)
	//	{
	//		BLOCK patch_rep(it->first);  // key
	//		size_t patch_depth;
	//		try 
	//		{
	//			auto pitr = find(patch_rep);

	//			if (pitr != end())
	//				patch_depth = pitr->second.mCompareDepth;
	//			else
	//				throw myexception();
	//		}
	//		catch (myexception& e)
	//		{
	//			std::cout << "key not found " << e.what() << std::endl;
	//		}
	//		
	//		auto unblocked_pt = UNBLOCK(patch_rep);
	//		point_2d pt;
	//		std::copy(std::begin((point_2d_POD &)unblocked_pt), std::end((point_2d_POD &)unblocked_pt), std::begin(pt));
	//		
	//		double interval_length = getIntervalLength(patch_depth);
	//		
	//		//std::cout << "rep " << pt[0] << ", " << pt[1] << std::endl;
	//		
	//		double x_interval_id(std::floor((pt[0] - 0.5) / interval_length));
	//		double y_interval_id(std::floor((pt[1] - 0.5) / interval_length));

	//		//std::cout << 0.5 + x_interval_id * interval_length << ", " << 0.5 + (x_interval_id + 1) * interval_length << std::endl;
	//		//std::cout << 0.5 + y_interval_id * interval_length << ", " << 0.5 + (y_interval_id + 1) * interval_length << std::endl;
	//		//std::cout << "\n";
	//		patches.push_back(0.5 + x_interval_id * interval_length);
	//		patches.push_back(0.5 + (x_interval_id + 1) * interval_length);
	//		patches.push_back(0.5 + y_interval_id * interval_length);
	//		patches.push_back(0.5 + (y_interval_id + 1) * interval_length);
	//	}
	//}

	size_t	getNoOfApproxFunctionCalls() const
	{
		return myNoOfApproxFunctionCalls;
	}

	size_t	getNoOfExactFunctionCalls() const
	{
		return myNoOfExactFunctionCalls;
	}

	size_t  getmyHitBottom() const{ return myHitBottom;  }
	
	void resetNoOfApproxFunctionCalls()
	{
		myNoOfApproxFunctionCalls = 0;
	}

	void resetNoOfExactFunctionCalls()
	{
		myNoOfExactFunctionCalls = 0;
	}

	/// evaluation
	double operator() (const point_2d & _pt, Datatype & data)
	{
#if _normalisation_
		point_2d transformed_pt = { sigmoid_transform(_pt[0]), sigmoid_transform(_pt[1]) };  // 2d sigmoid transform to [1, 2)

		point_2d_POD pt;
		std::copy(std::begin(transformed_pt), std::end(transformed_pt), std::begin(pt));
#else
		point_2d_POD pt;
		std::copy(std::begin(_pt), std::end(_pt), std::begin(pt));
#endif 
		myPTCOUNT++;

		BLOCK pt_block(pt);

		if (empty())
			return SavePointToCollection(data, pt_block, pt, DEFAULTDEPTH);

		BaseIterator rnbrptr = upper_bound(pt_block);	// O(log) complexity
		BaseIterator lnbrptr(rnbrptr);
		//myrit = rnbrptr;

		if ((rnbrptr != end()) && (rnbrptr != begin()))
		{
			lnbrptr--;
			//std::cout << "in 1" << std::endl;
			return CaseTwoNbrs(lnbrptr, pt_block, pt, rnbrptr, data);
		}
		else if (rnbrptr == begin())
		{
			//std::cout << "in 2" << std::endl;
			return CaseOneNbr(pt_block, pt, rnbrptr, RIGHTNBR, data);			// right nbr
		}
		else
		{
			lnbrptr--;
			//std::cout << "in 3" << std::endl;
			return CaseOneNbr(pt_block, pt, lnbrptr, LEFTNBR, data);			// left nbr
		}
	}

	void TEST_print_subdivide()
	{
		std::set<size_t> depths;
		for (auto itr = begin(); itr != end(); ++itr)
			depths.insert(itr->second.mCompareDepth);
		for (auto itr = depths.begin(); itr != depths.end(); ++itr)
			std::cout << *itr << " ";
		std::cout << std::endl;
	}

private:

	////////////////////////////////////////////////////////////////////////////////////////////
	//// Algorithm methods

	inline void Advance(BaseIterator & it)
	{
		for (auto n= myMaxPtsPerBox; n != 0 && it != end(); n--, it++);
	}

	inline void Retreat(BaseIterator & it)
	{
		for (auto n=myMaxPtsPerBox+2; n != 0 && it != begin(); n--, it--);
	}

	double CaseOneNbr( const BLOCK & pt_block, const point_2d_POD & pt, const BaseIterator & nbrptr, const bool nbrindicator, Datatype & data ) 
	{
		COMPARE less( nbrptr->second.mCompareDepth );
		if ( (nbrindicator==RIGHTNBR)? less( pt_block, nbrptr->first ) : less(nbrptr->first, pt_block) )
			return NotInNbrBox ( nbrptr, pt_block, pt, nbrindicator, data );   // if distinguishable then they are not in the same box
		else 
			return InNbrBox(nbrptr, pt_block, pt, data);
	}

	double CaseTwoNbrs ( const BaseIterator & lhs_it, const BLOCK & pt_block, const point_2d_POD & pt, const BaseIterator & rhs_it, Datatype & data ) 
	{

		COMPARE lless(lhs_it->second.mCompareDepth);
		COMPARE rless(rhs_it->second.mCompareDepth);

		auto getnewlessdepth = [&]()->size_t {
			size_t newlessdepth = std::min(lhs_it->second.mCompareDepth, rhs_it->second.mCompareDepth);
			while (COMPARE(newlessdepth)(lhs_it->first, pt_block) && COMPARE(newlessdepth)(pt_block, rhs_it->first))
				++newlessdepth;
			return newlessdepth-1;
		};

		// check if they have the same tree root
		size_t max_res = sizeof(double)*CHAR_BIT - 1;
		COMPARE lessrootcheck(max_res);

		if( lessrootcheck(lhs_it->first, rhs_it->first) )
		{
			// not the same root 
			if (lessrootcheck(lhs_it->first, pt_block))
				if(lessrootcheck(pt_block, rhs_it->first))
					return SavePointToCollection(data, pt_block, pt, max_res);	// new root, start a forest
				else 
					return CaseOneNbr(pt_block, pt, rhs_it, RIGHTNBR, data);
			else 
				return CaseOneNbr(pt_block, pt, lhs_it, LEFTNBR, data);
		}
		else
		{
			// same root
			if ( lless( lhs_it->first, pt_block) )
				if (rless(pt_block, rhs_it->first))
				{
					//std::cout << "in 1" << std::endl;
					return SavePointToCollection(data, pt_block, pt, getnewlessdepth());	// leaf node
				}
				else {
					//std::cout << "in 2" << std::endl;
					return InNbrBox(rhs_it, pt_block, pt, data);
				}
			else {
				//std::cout << "in 3" << std::endl;
				return InNbrBox(lhs_it, pt_block, pt, data);
			}
		}
	}
	
	double InNbrBox( const BaseIterator & nbritr, const BLOCK & pt_block, const point_2d_POD & pt, Datatype & data )
	{
		if ( nbritr->second.mIsInterpolated == false )								
		{
			// this next line invalidates any previously saved iterators
			// alters mysavedit to point to pt_block's iterator
			double ans = SavePointToCollection(data, pt_block, pt, nbritr->second.mCompareDepth);

			COMPARE less(nbritr->second.mCompareDepth);
			auto const itlow = getLowerbound(mysavedit, less);
			auto const itupper = getUpperbound(mysavedit, less);

			//if ((size_t)std::distance(itlow, itupper) > myMaxPtsPerBox)
			//	std::cout << "error pang pang " << std::endl;
			
			if ((size_t)std::distance(itlow, itupper) >= myMaxPtsPerBox)	// checks the number of data samples in the box is enough to do approximation
			{
				//std::cout << "here" << std::endl;
				if (Interpolation(pt_block, itlow, itupper))	// if interpolation is successful and accurate
					//ModifyBoxElements(itlow, itupper, nbritr->second.mCompareDepth, pt_block, true);
				{
					//std::cout << "in 1" << std::endl;
					std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
					base_value.second.mCompareDepth = nbritr->second.mCompareDepth;
					base_value.second.mIsInterpolated = true;
					base_value.second.mInterpolationKey = pt_block;
				});
				}
				else 
				{
					// find depth such that IsTriggerInterpolation no longer holds
					auto newitlow(itlow);
					auto newitupp(itupper);
					auto newdepth(nbritr->second.mCompareDepth);

					while ((size_t)std::distance(newitlow, newitupp) >= myMaxPtsPerBox){
						COMPARE newless(newdepth-1);
						newitlow = getLowerbound(mysavedit, newless);
						newitupp = getUpperbound(mysavedit, newless);
						newdepth--;
					}
					std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
						base_value.second.mCompareDepth = newdepth;
					});
				}
			}
			return ans;
		} 
		else // Use interpolated function
		{
			point_2d temp = { pt[0], pt[1] };
			++myNoOfApproxFunctionCalls;		/// TEST
			return myMapOfInterpolations[nbritr->second.mInterpolationKey].approximateFunction( temp );		// this is O(log(size()))
		}
	}

	double NotInNbrBox(const BaseIterator & nbr_it, const BLOCK & pt_block, const point_2d_POD & pt, const bool nbrindicator, Datatype & data)
	{
		// need the lowest common ancestor for determining either when nbr can be enlarged
		// or when nbr can't be enlarged but what depth do we assign to pt
		size_t new_depth = (nbrindicator == RIGHTNBR)?  GetLowestCommonAncestor(pt_block, nbr_it->first,nbr_it->second.mCompareDepth ) :
														GetLowestCommonAncestor(nbr_it->first, pt_block, nbr_it->second.mCompareDepth);	

		//double ans = 0;

		if (new_depth == nbr_it->second.mCompareDepth)		// no change means no common ancestor, start a new root
			return SavePointToCollection(data, pt_block, pt, DEFAULTDEPTH);		
		else
			return SavePointToCollection(data, pt_block, pt, new_depth - 1);		// set depth as the max possible that doesn't include nbr

		//// otherwise nbr and pt have the same root so execute the following code
		//if ( nbr_it->second.mIsEnlargeable )				/* assume interpolation always go up to the max domain hence we don't need to worry if the nbr box is interpolated or not */
		//{
		//	COMPARE new_less(new_depth);	// compare at common ancestors resolution
		//	ans = SavePointToCollection(data, pt_block, pt, nbr_it->second.mCompareDepth, ENLARGEABLE );
		//
		//	if ( IsTriggerInterpolation( pt_block, new_less) )
		//	{
		//		//if ( Interpolation(pt_block, less) )			// BUG? should this less be new_less?
		//		if ( Interpolation(pt_block, new_less) )
		//			ModifyBoxElements(pt_block, new_less, new_depth, pt_block, true, NOTENLARGEABLE);
		//		else	// else subdive
		//		{
		//			//TEST_print_subdivide(pt_block, new_less);
		//			ModifyBoxElements(pt_block, new_less, new_depth - 1, NOTENLARGEABLE);
		//		}
		//	} else 
		//		ModifyBoxElements(pt_block, new_less, new_depth, ENLARGEABLE);			// update elements' box size value
		//} else 
		//	ans = SavePointToCollection(data, pt_block, pt, new_depth-1);		// set depth as the max possible that doesn't include nbr
		//
		//return ans;

	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	//// tools
//	
//	bool Interpolation(const BLOCK & pt_block, COMPARE & compare)
//	{
//		// copy values into vector<double> for passing onto interpolation function
//		auto lowerbound			= getLowerbound(pt_block, compare);
//		auto upperbound			= getUpperbound(pt_block, compare);
//		const size_t no_of_pts  = no_in_box(lowerbound, upperbound);
//		
//		if (no_of_pts != myMaxPtsPerBox)
//			std::cout << no_of_pts << " " << myMaxPtsPerBox << std::endl;
//		
//		std::vector<point_2d>	all_pts(no_of_pts);
//		std::vector<double>		all_values(no_of_pts);
//
//		std::vector<point_2d>	interpolation_pts(no_of_pts / 2 );
//		std::vector<double>		interpolation_values(no_of_pts / 2 );
//		
//		size_t vec_index = 0;
//		for(auto it = lowerbound; it!=upperbound; ++it, ++vec_index)
//		{
//			auto unblocked_pt = UNBLOCK(it->first);
//			std::copy( std::begin(point_2d_POD(unblocked_pt)), std::end(point_2d_POD(unblocked_pt)), std::begin(all_pts[vec_index]));
//			//std::for_each(std::begin(interpolation_pts), std::end(interpolation_pts), [](point_2d & arg)->void { std::log(arg[0]); std::log(arg[1]); } );
//			all_values[vec_index] = it->second.mExactFunctionEvaluation;
//		}
//		InterpolationMorton approx(myDegree);
//
//		for(size_t i = 0; i < no_of_pts /2; ++i)
//		{
//			interpolation_pts[i] = all_pts[i *2];
//			interpolation_values[i] = all_values[i *2];
//		}
//
//#ifdef EXPONENTIATE
//		std::vector<point_2d>	interpolation_pts_log_version(no_of_pts);
//		std::transform( std::begin(all_pts), std::end(all_pts), std::begin(interpolation_pts_log_version),
//			[] (const point_2d & arg)->point_2d {  point_2d result; result[0] = std::log(arg[0]); result[1] = std::log(arg[1]); return result;} );
//
//		auto isSuccessful = approx.interpolate(interpolation_pts_log_version,all_values);	// success means successfully solved linear equations
//#else
//		auto isSuccessful = approx.interpolate(interpolation_pts, interpolation_values);	// success means successfully solved linear equations
//#endif
//
//		/// test for accuracy 
//		if (isSuccessful && approx.isAccurate(all_pts, all_values, myError) ) 
//		{
//			myMapOfInterpolations.insert( MapOfInterpolationsValueType( pt_block, approx ) );
//			isSuccessful = true;
//		}
//		else
//			isSuccessful = false;
//
//		return isSuccessful;
//	}

	bool Interpolation(const BLOCK & pt_block, const BaseIterator & itlow, const BaseIterator & itupper)
	{
		// copy values into vector<double> for passing onto interpolation function
		const size_t no_of_pts =  std::distance(itlow, itupper); // myMaxPtsPerBox; //
		
		const size_t no_of_training_pts = (this->myDegree + 2) * (this->myDegree + 1);
		//const size_t no_of_test_pts = no_of_pts - no_of_training_pts;

		if (no_of_pts != myMaxPtsPerBox)
			std::cout << no_of_pts << " " << myMaxPtsPerBox << std::endl;
		//std::cout << "no_of_pts " << no_of_pts << std::endl;

		std::vector<point_2d>	all_pts(no_of_pts);
		std::vector<double>		all_values(no_of_pts);

		std::vector<point_2d>	interpolation_pts(no_of_training_pts);
		std::vector<double>		interpolation_values(no_of_training_pts);

		int vec_index = 0;
		for (auto it = itlow; it != itupper; ++it, ++vec_index)
		{
			//std::cout << vec_index << " ";
			auto unblocked_pt = UNBLOCK(it->first);
			//std::copy(std::begin(point_2d_POD(unblocked_pt)), std::end(point_2d_POD(unblocked_pt)), std::begin(all_pts[vec_index]));
			std::copy(std::begin((point_2d_POD &)unblocked_pt), std::end((point_2d_POD &) unblocked_pt), std::begin(all_pts[vec_index]));
			//std::for_each(std::begin(interpolation_pts), std::end(interpolation_pts), [](point_2d & arg)->void { std::log(arg[0]); std::log(arg[1]); } );
			all_values[vec_index] = it->second.mExactFunctionEvaluation;
		}
		InterpolationMorton approx(myDegree);

		int j(0);
		for (int i = 0; i < no_of_pts && j < no_of_training_pts; j += 2)
		{
			interpolation_pts[j] = all_pts[i];
			interpolation_values[j] = all_values[i];

			interpolation_pts[j+1] = all_pts[i+2];
			interpolation_values[j+1] = all_values[i+2];

			i += 3;
		}

		auto isSuccessful = approx.interpolate_svd(interpolation_pts, interpolation_values);	// success means successfully solved linear equations
		//auto isSuccessful = approx.interpolate_svd(all_pts, all_values);	// success means successfully solved linear equations

		/// test for accuracy 
		bool maxdepthflag = itlow->second.mCompareDepth == DEFAULTMAXDEPTH;

		if (maxdepthflag)
			myHitBottom++;

		if ( (isSuccessful && approx.isAccurate(all_pts, all_values, myError)) || maxdepthflag )
		{
			myMapOfInterpolations.insert(MapOfInterpolationsValueType(pt_block, approx));
			isSuccessful = true;
		}
		else
			isSuccessful = false;

		return isSuccessful;
	}

	size_t GetLowestCommonAncestor(const BLOCK & A, const BLOCK & B, const size_t & res) const
	{
		size_t max_res = sizeof(double)* CHAR_BIT - 1 ;		// for determining whether it's possible to have common ancestor; -1 is to discount sgn bit
		size_t temp_depth = res;
		bool notbigenough(true);

		COMPARE maxless(max_res);
		if ( maxless(A, B))
			notbigenough = false;			// no common ancestors, eg, in 1 D: 0.1 and -0.1 won't be in the same box

		while ( notbigenough )
		{
			COMPARE less(++temp_depth);
			notbigenough = less(A, B);		// true if A < B, hence A and B are not in the same box according to resolution temp_depth
		}
		return temp_depth;
	}

	//inline Base::iterator getLowerbound(const BLOCK & pt, COMPARE & compare)
	//{
	//	auto iterator_compare_lowerbound = [&compare] (const BaseValueType & val, const BLOCK & ref_point)->bool
	//	{
	//		return compare(val.first,ref_point);
	//	};

	//	//return std::lower_bound(begin(), end(), pt, iterator_compare_lowerbound);
	//	return std::lower_bound(mylit, myrit, pt, iterator_compare_lowerbound);  // this is still linear!
	//}

	inline Base::iterator getLowerbound(const Base::iterator & pt_itr, COMPARE & compare)
	//Returns an iterator pointing to the first element in the range [first, last) that is not less than 
	//(i.e. greater or equal to) value, or last if no such element is found.
	{
		BaseIterator litr(pt_itr);  // make a copy of the reference iterator
		if ( !compare(begin()->first, pt_itr->first) )
			return begin();
		while (!compare(litr->first, pt_itr->first))
			--litr;
		return ++litr;
	}

	inline Base::iterator getUpperbound(const Base::iterator & pt_itr, COMPARE & compare)
	//Returns an iterator pointing to the first element in the range [first, last) that is greater than value, 
	//or last if no such element is found.
	{
		BaseIterator ritr(pt_itr);  // make a copy of the reference iterator
		while (ritr != end() && !compare(pt_itr->first, ritr->first))
			++ritr;
		return ritr;
	}

	//inline Base::iterator getUpperbound(const BLOCK & pt, COMPARE & compare)
	//{
	//	auto iterator_compare_upperbound = [&compare] (const BLOCK & ref_point, const BaseValueType & val)->bool
	//	{
	//		return compare(ref_point, val.first);
	//	};

	//	//return std::upper_bound(begin(), end(), pt, iterator_compare_upperbound);
	//	return std::upper_bound(mylit, myrit, pt, iterator_compare_upperbound);
	//}

	//inline size_t no_in_box(const BLOCK & pt, COMPARE & compare) const
	//{
	//	return std::distance(getLowerbound(pt, compare), getUpperbound(pt, compare));
	//}
	
	/*inline size_t no_in_box(const Base::iterator & lower, const Base::iterator & upper)
	{
		return std::distance(lower, upper);
	}*/

	//void ModifyBoxElements( const BLOCK & pt_block, COMPARE & less, const size_t & depth)
	//{
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each( getLowerbound(pt_block, less), getUpperbound(pt_block,less), [&] ( BaseValueType & base_value ) { 
	//								base_value.second.mCompareDepth = depth; 
	//							} );
	//}

	//inline void ModifyBoxElements(const BaseIterator & itlow, const BaseIterator & itupper, const size_t & depth)
	//{
	//	// doesn't change interpolation values
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
	//		base_value.second.mCompareDepth = depth;
	//	});
	//}

	//void ModifyBoxElements( const BLOCK & pt_block, COMPARE & less, const size_t & depth, const BLOCK key, const bool isinterpolated)
	//{
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each( getLowerbound(pt_block, less), getUpperbound(pt_block,less), [&] ( BaseValueType & base_value ) { 
	//								base_value.second.mCompareDepth = depth; 
	//								base_value.second.mIsInterpolated = isinterpolated;
	//								base_value.second.mInterpolationKey = key;
	//							} );
	//}

	//inline ModifyBoxElements(const BaseIterator & itlow, const BaseIterator & itupper, const size_t & depth, const BLOCK key, 
	//						const bool isinterpolated)
	//{
	//	// does change interpolation information
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
	//		base_value.second.mCompareDepth = depth;
	//		base_value.second.mIsInterpolated = isinterpolated;
	//		base_value.second.mInterpolationKey = key;
	//	});
	//}

	//bool IsTriggerInterpolation(const BLOCK & point, COMPARE & less)
	//{
	//	return (no_in_box(point, less) >= myMaxPtsPerBox);		//*****  FOR INTERPOLATION USE HALF OF THE VALUES, AND USE THE REST FOR TESTING
	//}
	//
	////bool IsTriggerInterpolation(const BaseIterator & itlow, const BaseIterator & itupper)
	//{
	//	return (no_in_box(itlow, itupper) >= myMaxPtsPerBox);		//*****  FOR INTERPOLATION USE HALF OF THE VALUES, AND USE THE REST FOR TESTING
	//}

	double SavePointToCollection(Datatype & data, const BLOCK & pt_block, const point_2d_POD & pt, size_t depth = DEFAULTDEPTH)
	{
		double ans = 0;
#if _normalisation_
		point_2d temp = { sigmoid_inv_transform(pt[0]),  sigmoid_inv_transform(pt[1]) };
#else
		point_2d temp = {pt[0], pt[1]};
#endif
		std::pair<BaseIterator, bool> ret;
		ret = insert(BaseValueType(pt_block, ApproxMortonValue(ans, depth)));
		mysavedit = ret.first;	// save the saved location
		if (ret.second == false) // if failed to insert then it already exists
		{
			ans = ret.first->second.mExactFunctionEvaluation;
		}
		else
		{
			ans = myFunction(temp, data);
			ret.first->second.mExactFunctionEvaluation = ans;
			++myNoOfExactFunctionCalls;		// TEST
		}
		return ans;
	}

};

#endif

#if dim_ == 1
template <typename Function, typename Datatype>
class ApproximationMorton : private std::map< BLOCK, ApproxMortonValue >
{
	typedef		std::map< BLOCK, InterpolationMorton>	MapOfInterpolations;
	typedef		std::map< BLOCK, ApproxMortonValue >	Base;
	typedef		Base::value_type						BaseValueType;
	typedef		Base::iterator							BaseIterator;
	typedef		MapOfInterpolations::value_type			MapOfInterpolationsValueType;

	using std::map< BLOCK, ApproxMortonValue >::lower_bound;
	using std::map< BLOCK, ApproxMortonValue >::upper_bound;

	double					myError;						
	Function&				myFunction;
	unsigned int			myDegree;
	MapOfInterpolations		myMapOfInterpolations;
	size_t					myMaxPtsPerBox;

	//BaseIterator	mylit;	// used to track box lowerbound upperbound search range
	//BaseIterator	myrit;	// since no box should contain more than myMaxPtsPerBox 
	BaseIterator	mysavedit;

	//////////////
	///	TESTING
	size_t					myNoOfApproxFunctionCalls;
	size_t					myNoOfExactFunctionCalls;
	size_t					myPTCOUNT;
	size_t					myHitBottom;
	/////////////
public:

	ApproximationMorton(Function& function, unsigned int degree, double abserror)
		: myError(abserror)
		, myFunction(function)
		, myDegree(degree)
		, myMaxPtsPerBox(8 * (1 + degree))  // was 2 for training, 1 for testing;  now 8 it's 4* for training, 4* for testing
		, myMapOfInterpolations()
		, myNoOfApproxFunctionCalls(0)
		, myNoOfExactFunctionCalls(0), myPTCOUNT(0), myHitBottom(0)
	{}

	~ApproximationMorton(void) {}

	size_t  getSizeOfPointCloud() const
	{
		return size();
	}

	size_t  getNoOfInterpolations() const
	{
		return myMapOfInterpolations.size();
	}

	void getInterpolatedPatches() const
	{
		for (MapOfInterpolations::const_iterator it = myMapOfInterpolations.begin();
			it != myMapOfInterpolations.end();
			++it)
		{
			BLOCK patch_rep(it->first);  // key
			size_t patch_depth;
			try
			{
				auto pitr = find(patch_rep);

				if (pitr != end())
					patch_depth = pitr->second.mCompareDepth;
				else
					throw myexception();
			}
			catch (myexception & e)
			{
				std::cout << "key not found " << e.what() << std::endl;
			}

			auto unblocked_pt = UNBLOCK(patch_rep);
			point_POD pt_;
			std::copy(std::begin((point_POD&)unblocked_pt), std::end((point_POD&)unblocked_pt), std::begin(pt_));
			auto interval = getInterval(pt_, patch_depth );

			std::cout << (interval.first) << ", " 
				<< (interval.second) << std::endl;
			
		}
	}

	void getAllPatches() const
	{
		for (auto it = begin();
			it != end();
			++it)
		{
			auto unblocked_pt = UNBLOCK(it->first);
			point_POD pt_;
			std::copy(std::begin((point_POD&)unblocked_pt), std::end((point_POD&)unblocked_pt), std::begin(pt_));
			auto interval = getInterval(pt_, it->second.mCompareDepth);

			std::cout << "pt: " << sigmoid_inv_transform(pt_[0]) << ", sigmoid(pt): " << pt_[0] << ", depth: " << it->second.mCompareDepth << ", ~ ["
				<< (interval.first) << ", " << (interval.second)
				<< "]" << std::endl;
		}
	}

	size_t	getNoOfApproxFunctionCalls() const
	{
		return myNoOfApproxFunctionCalls;
	}

	size_t	getNoOfExactFunctionCalls() const
	{
		return myNoOfExactFunctionCalls;
	}

	size_t  getmyHitBottom() const 
	{ 
		return myHitBottom; 
	}

	void resetNoOfApproxFunctionCalls()
	{
		myNoOfApproxFunctionCalls = 0;
	}

	void resetNoOfExactFunctionCalls()
	{
		myNoOfExactFunctionCalls = 0;
	}

	/// evaluation
	double operator() (const point & _pt, Datatype& data)
	{
		//return myFunction(_pt);

		/* 
		Two sections of the recode require patch renormalisation: 
			1) _pt renormalised to [1, 2] using sigmoid,  so all patches have root node [1, 2]
			2) computing approximation and evaluation -- myApproximateF(...)  and Interpolation(...), need to renormalise to [-1, 1]

		To do 1), we apply sigmoid_transform(_pt), so all patches are 'exp' scale. Thus when calling ExactFunction, we need to map back to original space
		using sigmoid_inv_transform to get f(x)
		
		So this means we need to then get the right scaling to map to [-1, 1];

		*/
		point transformed_pt = { sigmoid_transform(_pt[0]) };   // so we are operating on [1, 2] and its dyadic subintervals
		
		point_POD pt;
		std::copy(std::begin(transformed_pt), std::end(transformed_pt), std::begin(pt));

		myPTCOUNT++;

		BLOCK pt_block(pt);

		if (empty())
			return SavePointToCollection(data, pt_block, pt, DEFAULTDEPTH);

		BaseIterator rnbrptr = upper_bound(pt_block);	// O(log) complexity
		BaseIterator lnbrptr(rnbrptr);
		//myrit = rnbrptr;

		if ((rnbrptr != end()) && (rnbrptr != begin()))
		{
			lnbrptr--;
			//std::cout << "in 1" << std::endl;
			return CaseTwoNbrs(lnbrptr, pt_block, pt, rnbrptr, data);
		}
		else if (rnbrptr == begin())
		{
			//std::cout << "in 2" << std::endl;
			return CaseOneNbr(pt_block, pt, rnbrptr, RIGHTNBR, data);			// right nbr
		}
		else
		{
			lnbrptr--;
			//std::cout << "in 3" << std::endl;
			return CaseOneNbr(pt_block, pt, lnbrptr, LEFTNBR, data);			// left nbr
		}
	}

	void TEST_print_subdivide()
	{
		std::set<size_t> depths;
		for (auto itr = begin(); itr != end(); ++itr)
			depths.insert(itr->second.mCompareDepth);
		for (auto itr = depths.begin(); itr != depths.end(); ++itr)
			std::cout << *itr << " ";
		std::cout << std::endl;
	}

private:

	////////////////////////////////////////////////////////////////////////////////////////////
	//// Algorithm methods

	inline void Advance(BaseIterator& it)
	{
		for (auto n = myMaxPtsPerBox; n != 0 && it != end(); n--, it++);
	}

	inline void Retreat(BaseIterator& it)
	{
		for (auto n = myMaxPtsPerBox + 2; n != 0 && it != begin(); n--, it--);
	}

	double CaseOneNbr(const BLOCK& pt_block, const point_POD& pt, const BaseIterator& nbrptr, const bool nbrindicator, Datatype& data)
	{
		COMPARE less(nbrptr->second.mCompareDepth);
		if ((nbrindicator == RIGHTNBR) ? less(pt_block, nbrptr->first) : less(nbrptr->first, pt_block))
			return NotInNbrBox(nbrptr, pt_block, pt, nbrindicator, data);   // if distinguishable then they are not in the same box
		else
			return InNbrBox(nbrptr, pt_block, pt, data);
	}

	double CaseTwoNbrs(const BaseIterator& lhs_it, const BLOCK& pt_block, const point_POD& pt, const BaseIterator& rhs_it, Datatype& data)
	{

		COMPARE lless(lhs_it->second.mCompareDepth);
		COMPARE rless(rhs_it->second.mCompareDepth);

		auto getnewlessdepth = [&]()->size_t {
			size_t newlessdepth = std::min(lhs_it->second.mCompareDepth, rhs_it->second.mCompareDepth);
			while (COMPARE(newlessdepth)(lhs_it->first, pt_block) && COMPARE(newlessdepth)(pt_block, rhs_it->first))
				++newlessdepth;
			return newlessdepth - 1;
		};

		// check if they have the same tree root
		//size_t max_res = 52; // sizeof(double)* CHAR_BIT - 1;
		COMPARE lessrootcheck(DEFAULT_RES);

		if (lessrootcheck(lhs_it->first, rhs_it->first))
		{
			// not the same root 
			if (lessrootcheck(lhs_it->first, pt_block))
				if (lessrootcheck(pt_block, rhs_it->first))
					return SavePointToCollection(data, pt_block, pt, DEFAULT_RES);	// new root, start a forest
				else
					return CaseOneNbr(pt_block, pt, rhs_it, RIGHTNBR, data);
			else
				return CaseOneNbr(pt_block, pt, lhs_it, LEFTNBR, data);
		}
		else
		{
			// same root
			if (lless(lhs_it->first, pt_block))
				if (rless(pt_block, rhs_it->first))
				{
					//std::cout << "in 1" << std::endl;
					return SavePointToCollection(data, pt_block, pt, getnewlessdepth());	// leaf node
				}
				else {
					//std::cout << "in 2" << std::endl;
					return InNbrBox(rhs_it, pt_block, pt, data);
				}
			else {
				//std::cout << "in 3" << std::endl;
				return InNbrBox(lhs_it, pt_block, pt, data);
			}
		}
	}

	double InNbrBox(const BaseIterator& nbritr, const BLOCK& pt_block, const point_POD& pt, Datatype& data)
	{
		if (nbritr->second.mIsInterpolated == false)
		{
			// this next line invalidates any previously saved iterators
			// alters mysavedit to point to pt_block's iterator
			double ans = SavePointToCollection(data, pt_block, pt, nbritr->second.mCompareDepth);

			COMPARE less(nbritr->second.mCompareDepth);
			auto const itlow   = getLowerbound(mysavedit, less);
			auto const itupper = getUpperbound(mysavedit, less);

			//if ((size_t)std::distance(itlow, itupper) > myMaxPtsPerBox)
			//	std::cout << "error pang pang " << std::endl;

			if ((size_t)std::distance(itlow, itupper) >= myMaxPtsPerBox)	// checks the number of data samples in the box is enough to do approximation
			{
				//std::cout << "here" << std::endl;
				if (Interpolation(pt_block, itlow, itupper, nbritr->second.mCompareDepth, pt))	// if interpolation is successful and accurate
					//ModifyBoxElements(itlow, itupper, nbritr->second.mCompareDepth, pt_block, true);
				{
					//std::cout << "in 1" << std::endl;
					std::for_each(itlow, itupper, [&](BaseValueType& base_value) {
						base_value.second.mCompareDepth = nbritr->second.mCompareDepth;
						base_value.second.mIsInterpolated = true;
						base_value.second.mInterpolationKey = pt_block;
						});
				}
				else
				{
					// find depth such that IsTriggerInterpolation no longer holds
					auto newitlow(itlow);
					auto newitupp(itupper);
					auto newdepth(nbritr->second.mCompareDepth);

					while ((size_t)std::distance(newitlow, newitupp) >= myMaxPtsPerBox) {
						COMPARE newless(newdepth - 1);
						newitlow = getLowerbound(mysavedit, newless);
						newitupp = getUpperbound(mysavedit, newless);
						newdepth--;
					}
					std::for_each(itlow, itupper, [&](BaseValueType& base_value) {
						base_value.second.mCompareDepth = newdepth;
						});
				}
			}
			return ans;
		}
		else // Use interpolated function
		{
			++myNoOfApproxFunctionCalls;		/// TEST
			// need to map data to [-1, 1]
			auto patch1d = getInterval(pt, nbritr->second.mCompareDepth);
			//double lo = patch1d.first;
			//double up = patch1d.second;   // [lo,  up)
			//point temp = {  2 * (pt[0] - lo) / (up - lo) - 1 };
			point temp = { chebyshev_interval_shift(pt[0], patch1d.first, nbritr->second.mCompareDepth) };
			return myMapOfInterpolations[nbritr->second.mInterpolationKey].approximateFunction(temp);		// this is O(log(size()))
		}
	}

	double NotInNbrBox(const BaseIterator& nbr_it, const BLOCK& pt_block, const point_POD& pt, const bool nbrindicator, Datatype& data)
	{
		// need the lowest common ancestor for determining either when nbr can be enlarged
		// or when nbr can't be enlarged but what depth do we assign to pt
		size_t new_depth = (nbrindicator == RIGHTNBR) ? GetLowestCommonAncestor(pt_block, nbr_it->first, nbr_it->second.mCompareDepth) :
			GetLowestCommonAncestor(nbr_it->first, pt_block, nbr_it->second.mCompareDepth);

		//double ans = 0;

		if (new_depth == nbr_it->second.mCompareDepth)		// no change means no common ancestor, start a new root
			return SavePointToCollection(data, pt_block, pt, DEFAULTDEPTH);
		else
			return SavePointToCollection(data, pt_block, pt, new_depth - 1);		// set depth as the max possible that doesn't include nbr

		//// otherwise nbr and pt have the same root so execute the following code
		//if ( nbr_it->second.mIsEnlargeable )				/* assume interpolation always go up to the max domain hence we don't need to worry if the nbr box is interpolated or not */
		//{
		//	COMPARE new_less(new_depth);	// compare at common ancestors resolution
		//	ans = SavePointToCollection(data, pt_block, pt, nbr_it->second.mCompareDepth, ENLARGEABLE );
		//
		//	if ( IsTriggerInterpolation( pt_block, new_less) )
		//	{
		//		//if ( Interpolation(pt_block, less) )			// BUG? should this less be new_less?
		//		if ( Interpolation(pt_block, new_less) )
		//			ModifyBoxElements(pt_block, new_less, new_depth, pt_block, true, NOTENLARGEABLE);
		//		else	// else subdive
		//		{
		//			//TEST_print_subdivide(pt_block, new_less);
		//			ModifyBoxElements(pt_block, new_less, new_depth - 1, NOTENLARGEABLE);
		//		}
		//	} else 
		//		ModifyBoxElements(pt_block, new_less, new_depth, ENLARGEABLE);			// update elements' box size value
		//} else 
		//	ans = SavePointToCollection(data, pt_block, pt, new_depth-1);		// set depth as the max possible that doesn't include nbr
		//
		//return ans;

	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//// tools
//	
//	bool Interpolation(const BLOCK & pt_block, COMPARE & compare)
//	{
//		// copy values into vector<double> for passing onto interpolation function
//		auto lowerbound			= getLowerbound(pt_block, compare);
//		auto upperbound			= getUpperbound(pt_block, compare);
//		const size_t no_of_pts  = no_in_box(lowerbound, upperbound);
//		
//		if (no_of_pts != myMaxPtsPerBox)
//			std::cout << no_of_pts << " " << myMaxPtsPerBox << std::endl;
//		
//		std::vector<point_2d>	all_pts(no_of_pts);
//		std::vector<double>		all_values(no_of_pts);
//
//		std::vector<point_2d>	interpolation_pts(no_of_pts / 2 );
//		std::vector<double>		interpolation_values(no_of_pts / 2 );
//		
//		size_t vec_index = 0;
//		for(auto it = lowerbound; it!=upperbound; ++it, ++vec_index)
//		{
//			auto unblocked_pt = UNBLOCK(it->first);
//			std::copy( std::begin(point_2d_POD(unblocked_pt)), std::end(point_2d_POD(unblocked_pt)), std::begin(all_pts[vec_index]));
//			//std::for_each(std::begin(interpolation_pts), std::end(interpolation_pts), [](point_2d & arg)->void { std::log(arg[0]); std::log(arg[1]); } );
//			all_values[vec_index] = it->second.mExactFunctionEvaluation;
//		}
//		InterpolationMorton approx(myDegree);
//
//		for(size_t i = 0; i < no_of_pts /2; ++i)
//		{
//			interpolation_pts[i] = all_pts[i *2];
//			interpolation_values[i] = all_values[i *2];
//		}
//
//#ifdef EXPONENTIATE
//		std::vector<point_2d>	interpolation_pts_log_version(no_of_pts);
//		std::transform( std::begin(all_pts), std::end(all_pts), std::begin(interpolation_pts_log_version),
//			[] (const point_2d & arg)->point_2d {  point_2d result; result[0] = std::log(arg[0]); result[1] = std::log(arg[1]); return result;} );
//
//		auto isSuccessful = approx.interpolate(interpolation_pts_log_version,all_values);	// success means successfully solved linear equations
//#else
//		auto isSuccessful = approx.interpolate(interpolation_pts, interpolation_values);	// success means successfully solved linear equations
//#endif
//
//		/// test for accuracy 
//		if (isSuccessful && approx.isAccurate(all_pts, all_values, myError) ) 
//		{
//			myMapOfInterpolations.insert( MapOfInterpolationsValueType( pt_block, approx ) );
//			isSuccessful = true;
//		}
//		else
//			isSuccessful = false;
//
//		return isSuccessful;
//	}

	/* Warning, this assumes pt is in the range [1, 2]
	*/
	std::pair<double, double> getInterval(const point_POD& pt, const size_t depth) const
	{
		try {
			if (pt[0] < 1 || pt[0] > 2)
				throw (pt[0]);
		}
		catch (double myNum) {
			std::cout << "getInterval: value outside [1,2]\n" << myNum;
			exit(-1);
		}

		union {
			double value;
			uint64_t value_int;
		};

		value = pt[0];
		value_int = (value_int >> depth) << depth;  // set all masked bits to 0;

		return { value,  value + 1. / double(1 << (DEFAULT_RES - depth)) };  
	}

	bool Interpolation(const BLOCK& pt_block, const BaseIterator& itlow, const BaseIterator& itupper, const size_t depth, const point_POD& pt)
	{
		//std::cout << "interpolation pt: " << pt[0] << ", depth " << depth << "~ [";
		auto patch1d = getInterval(pt, depth);
		//std::cout << patch1d.first << ", " << patch1d.second << "]\n";

		// copy values into vector<double> for passing onto interpolation function
		const size_t no_of_pts = std::distance(itlow, itupper); // myMaxPtsPerBox; //

		const size_t no_of_training_pts = 4 * (this->myDegree + 1);
		//const size_t no_of_test_pts = no_of_pts - no_of_training_pts;

		if (no_of_pts != myMaxPtsPerBox)
			std::cout << no_of_pts << " " << myMaxPtsPerBox << std::endl;
		//std::cout << "no_of_pts " << no_of_pts << std::endl;

		std::vector<point>	all_pts(no_of_pts);
		std::vector<double>	all_values(no_of_pts);

		std::vector<point>	interpolation_pts(no_of_training_pts);
		std::vector<double>	interpolation_values(no_of_training_pts);

		int vec_index = 0;
		for (auto it = itlow; it != itupper; ++it, ++vec_index)
		{
			//std::cout << vec_index << " ";
			auto unblocked_pt = UNBLOCK(it->first);
			std::copy(std::begin((point_POD&)unblocked_pt), std::end((point_POD&)unblocked_pt), std::begin(all_pts[vec_index]));
			//remap [-1,1]
			all_pts[vec_index][0] = chebyshev_interval_shift(all_pts[vec_index][0], patch1d.first, depth);
			//std::cout << "x: " << all_pts[vec_index][0] << ", low: " << patch1d.first << ", depth: " << depth << std::endl
			//	<< "ans: " << all_pts[vec_index][0] << std::endl << std::endl;
				//2 * (all_pts[vec_index][0] - lo) / (up - lo) - 1;
			all_values[vec_index] = it->second.mExactFunctionEvaluation;
		}
		InterpolationMorton approx(myDegree);

		int j(0);
		// j += number of training pts;  i += # training + # testing
		for (int i = 0; i < no_of_pts && j < no_of_training_pts; j += 4)
		{
			interpolation_pts[j] = all_pts[i];
			interpolation_values[j] = all_values[i];

			interpolation_pts[j + 1] = all_pts[i + 2];
			interpolation_values[j + 1] = all_values[i + 2];

			interpolation_pts[j + 2] = all_pts[i + 4];
			interpolation_values[j + 2] = all_values[i + 4];

			interpolation_pts[j + 3] = all_pts[i + 6];
			interpolation_values[j + 3] = all_values[i + 6];

			i += 8;
		}

		// for 3 (2 training, 1 test) it's the following
		/*for (int i = 0; i < no_of_pts && j < no_of_training_pts; j += 2)
		{
			interpolation_pts[j] = all_pts[i];
			interpolation_values[j] = all_values[i];

			interpolation_pts[j + 1] = all_pts[i + 2];
			interpolation_values[j + 1] = all_values[i + 2];

			i += 3;
		}*/

		auto isSuccessful = approx.interpolate_svd(interpolation_pts, interpolation_values);	// success means successfully solved linear equations
		//auto isSuccessful = approx.interpolate_svd(all_pts, all_values);	// success means successfully solved linear equations

		/// test for accuracy 
		bool maxdepthflag = itlow->second.mCompareDepth == DEFAULTMAXDEPTH;

		if (maxdepthflag)
			myHitBottom++;

		if ((isSuccessful && approx.isAccurate(all_pts, all_values, myError)) || maxdepthflag)
		{
			myMapOfInterpolations.insert(MapOfInterpolationsValueType(pt_block, approx));
			isSuccessful = true;
		}
		else
			isSuccessful = false;

		return isSuccessful;
	}

	size_t GetLowestCommonAncestor(const BLOCK& A, const BLOCK& B, const size_t& res) const
	{
		//size_t max_res = sizeof(double) * CHAR_BIT - 1;		// for determining whether it's possible to have common ancestor; -1 is to discount sgn bit
		size_t temp_depth = res;
		bool notbigenough(true);

		COMPARE maxless(DEFAULT_RES);
		if (maxless(A, B))
			notbigenough = false;			// no common ancestors, eg, in 1 D: 0.1 and -0.1 won't be in the same box

		while (notbigenough)
		{
			COMPARE less(++temp_depth);
			notbigenough = less(A, B);		// true if A < B, hence A and B are not in the same box according to resolution temp_depth
		}
		return temp_depth;
	}

	//inline Base::iterator getLowerbound(const BLOCK & pt, COMPARE & compare)
	//{
	//	auto iterator_compare_lowerbound = [&compare] (const BaseValueType & val, const BLOCK & ref_point)->bool
	//	{
	//		return compare(val.first,ref_point);
	//	};

	//	//return std::lower_bound(begin(), end(), pt, iterator_compare_lowerbound);
	//	return std::lower_bound(mylit, myrit, pt, iterator_compare_lowerbound);  // this is still linear!
	//}

	inline Base::iterator getLowerbound(const Base::iterator& pt_itr, COMPARE& compare)
		//Returns an iterator pointing to the first element in the range [first, last) that is not less than 
		//(i.e. greater or equal to) value, or last if no such element is found.
	{
		BaseIterator litr(pt_itr);  // make a copy of the reference iterator
		if (!compare(begin()->first, pt_itr->first))
			return begin();
		while (!compare(litr->first, pt_itr->first))
			--litr;
		return ++litr;
	}

	inline Base::iterator getUpperbound(const Base::iterator& pt_itr, COMPARE& compare)
		//Returns an iterator pointing to the first element in the range [first, last) that is greater than value, 
		//or last if no such element is found.
	{
		BaseIterator ritr(pt_itr);  // make a copy of the reference iterator
		while (ritr != end() && !compare(pt_itr->first, ritr->first))
			++ritr;
		return ritr;
	}

	//inline Base::iterator getUpperbound(const BLOCK & pt, COMPARE & compare)
	//{
	//	auto iterator_compare_upperbound = [&compare] (const BLOCK & ref_point, const BaseValueType & val)->bool
	//	{
	//		return compare(ref_point, val.first);
	//	};

	//	//return std::upper_bound(begin(), end(), pt, iterator_compare_upperbound);
	//	return std::upper_bound(mylit, myrit, pt, iterator_compare_upperbound);
	//}

	//inline size_t no_in_box(const BLOCK & pt, COMPARE & compare) const
	//{
	//	return std::distance(getLowerbound(pt, compare), getUpperbound(pt, compare));
	//}

	/*inline size_t no_in_box(const Base::iterator & lower, const Base::iterator & upper)
	{
		return std::distance(lower, upper);
	}*/

	//void ModifyBoxElements( const BLOCK & pt_block, COMPARE & less, const size_t & depth)
	//{
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each( getLowerbound(pt_block, less), getUpperbound(pt_block,less), [&] ( BaseValueType & base_value ) { 
	//								base_value.second.mCompareDepth = depth; 
	//							} );
	//}

	//inline void ModifyBoxElements(const BaseIterator & itlow, const BaseIterator & itupper, const size_t & depth)
	//{
	//	// doesn't change interpolation values
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
	//		base_value.second.mCompareDepth = depth;
	//	});
	//}

	//void ModifyBoxElements( const BLOCK & pt_block, COMPARE & less, const size_t & depth, const BLOCK key, const bool isinterpolated)
	//{
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each( getLowerbound(pt_block, less), getUpperbound(pt_block,less), [&] ( BaseValueType & base_value ) { 
	//								base_value.second.mCompareDepth = depth; 
	//								base_value.second.mIsInterpolated = isinterpolated;
	//								base_value.second.mInterpolationKey = key;
	//							} );
	//}

	//inline ModifyBoxElements(const BaseIterator & itlow, const BaseIterator & itupper, const size_t & depth, const BLOCK key, 
	//						const bool isinterpolated)
	//{
	//	// does change interpolation information
	//	// range affected by for_each is [lowerbound, upperbound)
	//	std::for_each(itlow, itupper, [&](BaseValueType & base_value) {
	//		base_value.second.mCompareDepth = depth;
	//		base_value.second.mIsInterpolated = isinterpolated;
	//		base_value.second.mInterpolationKey = key;
	//	});
	//}

	//bool IsTriggerInterpolation(const BLOCK & point, COMPARE & less)
	//{
	//	return (no_in_box(point, less) >= myMaxPtsPerBox);		//*****  FOR INTERPOLATION USE HALF OF THE VALUES, AND USE THE REST FOR TESTING
	//}
	//
	////bool IsTriggerInterpolation(const BaseIterator & itlow, const BaseIterator & itupper)
	//{
	//	return (no_in_box(itlow, itupper) >= myMaxPtsPerBox);		//*****  FOR INTERPOLATION USE HALF OF THE VALUES, AND USE THE REST FOR TESTING
	//}

	double SavePointToCollection(Datatype& data, const BLOCK& pt_block, const point_POD& pt, size_t depth = DEFAULTDEPTH)
	{
		double ans = 0;
		point temp = { sigmoid_inv_transform(pt[0]) };	// pt[0] is sigmoid transformed value
		std::pair<BaseIterator, bool> ret;
		ret = insert(BaseValueType(pt_block, ApproxMortonValue(ans, depth)));
		mysavedit = ret.first;	// save the saved location
		if (ret.second == false) // if failed to insert then it already exists
		{
			ans = ret.first->second.mExactFunctionEvaluation;
		}
		else
		{
			ans = myFunction(temp);
			ret.first->second.mExactFunctionEvaluation = ans;
			++myNoOfExactFunctionCalls;		// TEST
		}
		return ans;
	}

};



#endif