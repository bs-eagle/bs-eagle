/// @file strategy_3d_common.h
/// @brief Common 3D strategy invariants independent from geometry implementation
/// @author uentity
/// @version 1.0
/// @date 14.02.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef STRATEGY_3D_COMMON_UMT7YBF5
#define STRATEGY_3D_COMMON_UMT7YBF5

#include "wpi_common.h"

namespace blue_sky { namespace wpi {

template< template< uint > class strat_traits >
struct strategy_3d_common {
	// dimens num, inner point id
	enum { D = 3, inner_point_id = 6 };

	typedef t_float vertex_pos[D];
	typedef ulong   vertex_pos_i[D];

	// iterator over source arrays come from traits
	typedef strat_traits< D > traits_t;
	typedef typename traits_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename traits_t::well_traj_iterator   well_traj_iterator;

	// X-Y-Z order!
	static void decode_cell_id(ulong id, vertex_pos_i& res, const vertex_pos_i& m_size) {
		//vertex_pos_i res;
		res[2] = id / (m_size[0] * m_size[1]);
		res[1] = (id - res[2] * m_size[0] * m_size[1]) / m_size[0];
		res[0] = id - m_size[0]*(res[2] * m_size[1] + res[1]);
	}

	static ulong encode_cell_id(const vertex_pos_i& p, const vertex_pos_i& m_size) {
		return p[0] + m_size[0] * (p[1] + p[2] * m_size[1]);
	}

	// This function designed for use in mesh_part::boundary
	// we need to cover mesh_part bbox with non-intersecting mesh_parts,
	// representing original mesh_part boundary
	// if we define boundaries like (example for X dimension):
	// lo_1 = [0, 0, 0], hi_1 = [1, n, n]
	// lo_2 = [n - 1, 0, 0], hi_2 = [n, n, n]
	// then what offsets should we add to each bounday to prevent them
	// from finally intersecting?
	// Purpose of this function is to return needed differencies for each boundary
	// in each dimension.
	// I decided to place it into strategy, because generic dimension-independant algo
	// for calculating boundary isn't obvious to me right now
	typedef int bbox_bnd_offs[2][D];
	static const bbox_bnd_offs& bbox_boundary_offs(const uint dim, const uint bnd_id) {
		static const bbox_bnd_offs t[6] = {
			// X
			{ {0, 1, 0}, {0,  0,  0} },
			{ {0, 0, 0}, {0, -1,  0} },
			// Y
			{ {0, 0, 0}, {-1, 0, 0} },
			{ {1, 0, 0}, { 0, 0, 0} },
			// Z
			{ {1, 1, 0}, {-1, -1, 0} },
			{ {1, 1, 0}, {-1, -1, 0} }
		};
		return t[dim*2 + bnd_id];
	}

	// well_data special part
	template< class base_t >
	struct well_data : public base_t {
		// ctors
		well_data() {}

		well_data(const well_traj_iterator& segment, const well_data* prev = NULL)
			: base_t(segment)
		{}

		// MD access
		using base_t::W;
		t_float md() const {
			return W[3];
		}

		t_float md_finish() const {
			return W[7];
		}
	};
};

}} // eof blue_sky::wpi

#endif /* end of include guard: STRATEGY_3D_COMMON_UMT7YBF5 */

