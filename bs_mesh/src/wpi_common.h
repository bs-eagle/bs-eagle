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

#include "conf.h"
//class bs_mesh_grdecl;

namespace blue_sky { namespace wpi {

// common typedefs
typedef t_ulong ulong;
typedef t_uint uint;
typedef v_float::iterator vf_iterator;

typedef t_float cell_pos[8][3];

// CGAL common typedefs
typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// assign for c arrays
// fun with returning reference to array :)
template< class T, int L >
static T (&ca_assign(T (&lhs)[L], const T (&rhs)[L]))[L] {
	std::copy(&rhs[0], &rhs[L], &lhs[0]);
	return lhs;
}

template< class T, int L >
static T (&ca_assign(T (&lhs)[L], const T& v))[L] {
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

// simple traits for strategies
// assuming that cells vertices and well path are represented as continous C-arrays
struct carray_traits {
	// iterator over raw tops array - simple pointer
	typedef t_float* cell_vertex_iterator;
	typedef t_float* well_traj_iterator;

	// generic converter from iterator to vertex_pos
	// return reference to array to prevent data copying
	template< class pos_t, class iterator_t >
	static pos_t& iter2pos(iterator_t& src) {
		// for simple arrays we can simply use reinterpret_cast
		return reinterpret_cast< pos_t& >(*src);
	}
};

}} /* { namespace blue_sky::wpi */

#endif /* end of include guard: WPI_COMMON_Y614D09T */

