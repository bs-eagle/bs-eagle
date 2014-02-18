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

	// common cell properties
	struct cell_data {
		enum { n_vertex = 8 };
		enum { n_facets = 6 };
		enum { n_facet_vertex = 4 };
		enum { n_edges = 12 };

		// array to hold facet vertices IDs
		typedef ulong facet_vid_t[n_facet_vertex];

		/* nodes layout
			*                             X
			*                    0+-------+1
			*                    /|     / |
			*                  /  |   /   |
			*               2/-------+3   |
			*              Y |   4+--|----+5
			*                |   /Z  |   /
			*                | /     | /
			*              6 /-------/7
		*/
		/*  facets layout (nodes - plane id)
		*  0-1-2-3 - 0
		*  0-1-4-5 - 1
		*  4-5-6-7 - 2
		*  2-3-6-7 - 3
		*  0-2-4-6 - 4
		*  1-3-5-7 - 5
		*  inside cell - 6
		*/

		// obtain facet idx from given axe ID (dimension)
		// and facet ID (0 - starting from origin, 1 - opposite facet)
		static ulong facet_id(ulong dim, ulong facet) {
			switch(dim) {
				// X
				case 0 :
					return 4 + facet;
				// Y
				case 1 :
					return 1 + facet*2;
				// Z
				case 2 :
					return 0 + facet*2;
			}
			return ulong(-1);
		}

		// obtain IDs of given facet vertices
		static const facet_vid_t& facet_vid(ulong facet) {
			static const facet_vid_t t[6] = {
				{0, 1, 3, 2}, {0, 1, 5, 4}, {4, 5, 7, 6},
				{2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5}
				//{ulong(-1), ulong(-1), ulong(-1), ulong(-1)}
			};
			return t[facet];
		}
		static void facet_vid(ulong dim, ulong facet, facet_vid_t& res) {
			return facet_vid(facet_id(dim, facet, res));
		}
	};

	// fast chack if bbox intersect with segment
	template< class Point >
	static bool bbox_segment_x(
		const Point& b_min, const Point& b_max,
		const Point& s_min, const Point& s_max
		)
	{
		double t_min[D], t_max[D];
		for(uint i = 0; i < D; ++i) {
			double dd = 0;
			if(s_min[i] != s_max[i])
				dd = 1.0 / (s_min[i] - s_max[i]);

			t_min[i] = (b_max[i] - s_min[i]) * dd;
			t_max[i] = (b_min[i] - s_min[i]) * dd;

			// fisrt check in 2D
			if(i == 1) {
				if(t_min[0] > t_max[1] || t_min[1] > t_max[0])
					return false;
				t_min[0] = std::max(t_min[0], t_min[1]);
				t_max[0] = std::min(t_max[0], t_max[1]);
			}
		}

		// last check in 3d
		if(t_min[0] > t_max[2] || t_min[2] > t_max[0])
			return false;
		return true;
	}

	template< class Bbox, class Segment >
	static bool bbox_segment_x(const Bbox& b, const Segment& s) {
		return bbox_segment_x(
			Point(b.min(0), b.min(1), b.min(2)),
			Point(b.max(0), b.max(1), b.max(2)),
			s.min(), s.max()
		);
	}
};

}} // eof blue_sky::wpi

#endif /* end of include guard: STRATEGY_3D_COMMON_UMT7YBF5 */

