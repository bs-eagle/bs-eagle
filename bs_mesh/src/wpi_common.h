/// @file wpi_common.h
/// @brief Global entities needed by WPI algos
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_COMMON_Y614D09T
#define WPI_COMMON_Y614D09T

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Object.h>

#include <iterator>
#include <cmath>
#include <boost/integer/static_min_max.hpp>

#include "conf.h"

namespace blue_sky { namespace wpi {

// common typedefs
typedef v_float::iterator vf_iterator;
typedef t_float cell_pos[8][3];
using blue_sky::ulong;
using blue_sky::uint;

// CGAL common typedefs
typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// assign for c arrays
// fun with returning reference to array :)
#if defined(_MSC_VER)
// dumb Visual Studio cannot calc min array length
template< class T, class X, int N >
static T (&ca_assign(T (&lhs)[N], const X (&rhs)[N]))[N] {
	std::copy(&rhs[0], &rhs[N], &lhs[0]);
	return lhs;
}
#else
// auto-deduce min array length in compile-time
template< class T, class X, int M, int N >
static T (&ca_assign(T (&lhs)[M], const X (&rhs)[N]))[boost::static_unsigned_min< M, N >::value] {
	std::copy(&rhs[0], &rhs[boost::static_unsigned_min< M, N >::value], &lhs[0]);
	return lhs;
}
#endif

template< class T, class X, int L >
static T (&ca_assign(T (&lhs)[L], const X& v))[L] {
	std::fill(&lhs[0], &lhs[L], v);
	return lhs;
}

// implementation of SGI extension algorithm
template< typename _InputIterator1, typename _InputIterator2 >
int lexicographical_compare_3way(_InputIterator1 __first1,
		_InputIterator1 __last1,
		_InputIterator2 __first2,
		_InputIterator2 __last2)
{
	while (__first1 != __last1 && __first2 != __last2)
	{
		if (*__first1 < *__first2)
			return -1;
		if (*__first2 < *__first1)
			return 1;
		++__first1;
		++__first2;
	}
	if (__first2 == __last2)
		return !(__first1 == __last1);
	else
		return -1;
}

// forward declaraions of strategy traits
struct carray_traits;
struct online_tops_traits;
struct online_tops_traits_bufpool;
struct sgrid_traits;

}} /* { namespace blue_sky::wpi */

#endif /* end of include guard: WPI_COMMON_Y614D09T */

