/// @file strategy_2d_common.h
/// @brief Common 2D strategy invariants independent from geometry implementation
/// @author uentity
/// @version 1.0
/// @date 14.02.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef STRATEGY_2D_COMMON_KTKT7HYH
#define STRATEGY_2D_COMMON_KTKT7HYH

#include "wpi_common.h"

namespace blue_sky { namespace wpi {

template< template< uint > class strat_traits >
struct strategy_2d_common {
	// dimens num, inner point id
	enum { D = 2, inner_point_id = 4 };

	typedef t_float vertex_pos[D];
	typedef ulong   vertex_pos_i[D];

	// iterator over source arrays come from traits
	typedef strat_traits< D > traits_t;
	typedef typename traits_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename traits_t::well_traj_iterator   well_traj_iterator;

	// X-Y-Z order!
	static void decode_cell_id(ulong id, vertex_pos_i& res, const vertex_pos_i& m_size) {
		//vertex_pos_i res;
		res[1] = id / m_size[0];
		res[0] = id - m_size[0]*res[1];
	}

	static ulong encode_cell_id(const vertex_pos_i& p, const vertex_pos_i& m_size) {
		return p[0] + m_size[0] * p[1];
	}

	// look into strategy_3d_common.h for description of this function
	typedef int bbox_bnd_offs[2][D];
	static const bbox_bnd_offs& bbox_boundary_offs(const uint dim, const uint bnd_id) {
		static const bbox_bnd_offs t[4] = {
			// X
			{ {0, 1}, {0,  0} },
			{ {0, 0}, {0, -1} },
			// Y
			{ {0, 0}, {-1, 0} },
			{ {1, 0}, { 0, 0} },
		};
		return t[dim*2 + bnd_id];
	}

	// well_data contains specifid MD member for
	// segment distances in 2D
	template< class base_t >
	struct well_data : public base_t {
		t_float md_;

		//empty ctor for map
		well_data() : md_(0) {}
		//std ctor
		well_data(const well_traj_iterator& segment, const well_data* prev = NULL)
			: base_t(segment), md_(0)
		{
			if(prev)
				md_ = prev->md_ + prev->len();
		}

		// MD access
		using base_t::W;
		using base_t::len;

		t_float md() const {
			return md_;
			//return W[3];
		}

		t_float md_finish() const {
			return md() + len();
		}
	};
};

}} // eof blue_sky::wpi

#endif /* end of include guard: STRATEGY_2D_COMMON_KTKT7HYH */

