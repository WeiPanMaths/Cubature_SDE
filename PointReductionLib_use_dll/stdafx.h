#ifndef stdafx_h__
#define stdafx_h__

#include <iostream>
#include <vector>
#include <valarray>

// TODO: reference additional headers your program requires here
namespace
{
template <class T, class A>
inline	std::ostream& operator <<(std::ostream& out, const std::vector<T, A>& arg)
{
	out << '{';
	for (size_t i = 0; i < arg.size(); ++i)
		out << (const T&)arg[i] << " ";
	out << '}';
	return out;
}

template <class T>
inline	std::ostream& operator <<(std::ostream& out, const std::valarray<T>& arg)
{
	out << '{';
	for (size_t i = 0; i < arg.size(); ++i)
		out << (const T&)arg[i] << " ";
	out << '}';
	return out;
}	

template <class TT>
inline const TT& deref(const void* const arg)
{
	return *(const TT* const)arg;
}

template <class T, class TT>
class dereference
{
	const T& container;

public:
	dereference(const T& arg)
		: container(arg) {};

	friend std::ostream& operator <<(std::ostream& out, const dereference& arg)
	{
		out << '{';
		std::transform(arg.container.begin(), arg.container.end(), std::ostream_iterator<TT>(out, " "), deref<TT>);
		out << '}';
		return out;
	}
};
}
inline double mm(const double& w, const void*& pl)
{
	return w * deref<double>(pl);
}

struct mmm
{
	const double& m_dMean;
	mmm(const double& dMean)
		: m_dMean(dMean) {};
	double operator ()(const double& w, const void*& pl)
	{
		return w * (deref<double>(pl) - m_dMean) * (deref<double>(pl) - m_dMean);
	}
};

template <class T>
T plus(const T lhs, const T rhs)
{
	return lhs + rhs;
}

#include "../test/SHOW.h"
// TODO: reference additional headers your program requires here
#endif // stdafx_h__
