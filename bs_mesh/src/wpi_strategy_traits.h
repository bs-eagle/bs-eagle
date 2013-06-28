/// @file wpi_online_tops_traits.h
/// @brief Online tops WPI strategy traits implementation
/// @author uentity
/// @version 1.0
/// @date 09.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ONLINE_TOPS_TRAITS_U2022LJY
#define WPI_ONLINE_TOPS_TRAITS_U2022LJY

#include "tops_iterator.h"

namespace blue_sky { namespace wpi {

// simple traits for strategies
// assuming that cells vertices and well path are represented as continous C-arrays
// D is dimensionality
// need to introduce this template parameter because iterator can depend on dimesionality
template< uint D_ >
struct carray_traits {
	// iterator over raw tops array - simple pointer
	typedef t_float* cell_vertex_iterator;
	typedef t_float* well_traj_iterator;
	enum { D = D_ };

	// generic converter from iterator to vertex_pos
	// return reference to array to prevent data copying
	template< class pos_t, class iterator_t >
	static pos_t& iter2pos(iterator_t& src) {
		// for simple arrays we can simply use reinterpret_cast
		return reinterpret_cast< pos_t& >(*src);
	}

	static const char* name() {
		return "carray";
	}
};

// cell's vertices coords are calculated on the fly when accessing given cell
template< uint D_ >
struct online_tops_traits {
	typedef tops_iterator< carray_ti_traits, D_ > cell_vertex_iterator;
	typedef t_float*                              well_traj_iterator;
	enum { D = D_ };

	template < typename T, size_t N > char (&array_sz(const T(&)[N]))[N];
	template < typename T, size_t N, size_t M > char (&array_sz(const T(&)[N][M]))[N * M];

	// generic converter from iterator to vertex_pos
	// return reference to array to prevent data copying
	template< class pos_t, class iterator_t >
	static pos_t& iter2pos(iterator_t& src) {
		//assert(
		//	src.fit2data(sizeof(array_sz(pos_t())))
		//);

		// now we can safely use reinterpret_cast
		return reinterpret_cast< pos_t& >(*src);
	}

	static const char* name() {
		return "online_tops";
	}
};

template< uint D_ >
struct online_tops_traits_bufpool : public online_tops_traits< D_ > {
	// add this typedef for tops_iterator::impl::iimpl
	// it needs to know corresponding uncached strategy
	typedef online_tops_traits< D_ > uncached_traits;
	typedef tops_iterator< bufpool_ti_traits, D_ > cell_vertex_iterator;

	static const char* name() {
		return "online_tops_bufpool";
	}
};

template< uint D_ >
struct sgrid_traits : public online_tops_traits< D_ > {
	typedef tops_iterator< sgrid_ti_traits, D_ > cell_vertex_iterator;

	static const char* name() {
		return "sgrid";
	}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ONLINE_TOPS_TRAITS_U2022LJY */

