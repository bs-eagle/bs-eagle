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

// cell's vertices coords are calculated on the fly when accessing given cell
struct online_tops_traits {
	typedef tops_iterator< carray_ti_traits > cell_vertex_iterator;
	typedef t_float*                          well_traj_iterator;

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
};

struct online_tops_traits_bufpool : public online_tops_traits {
	typedef tops_iterator< bufpool_ti_traits > cell_vertex_iterator;
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ONLINE_TOPS_TRAITS_U2022LJY */

