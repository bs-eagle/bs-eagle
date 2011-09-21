/// @file wpi_strategy_3d.h
/// @brief 3D strategy for well path identification algorithms
/// @author uentity
/// @version 
/// @date 13.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGY_3D_SJLKT8NL
#define WPI_STRATEGY_3D_SJLKT8NL

#include "wpi_common.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Object.h>
#include <CGAL/Bbox_3.h>
//#include <CGAL/box_intersection_d.h>
//#include <CGAL/intersections.h>

//#include "conf.h"

namespace blue_sky { namespace wpi {

struct wpi_strategy_3d {
	// common typedefs
	//typedef t_ulong ulong;
	//typedef t_uint uint;

	// main typedefs
	//typedef CGAL::Object                                        Object;
	//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
	typedef Kernel::Point_3                                     Point;
	typedef Kernel::Segment_3                                   Segment;
	typedef CGAL::Bbox_3                                        Bbox;
	typedef Kernel::Iso_cuboid_3                                Iso_bbox;

	// 3D specific typedefs
	typedef Kernel::Triangle_3                                  Triangle;
	typedef Kernel::Tetrahedron_3                               Tetrahedron;
	typedef std::vector<Triangle>                               Triangles;
	typedef Triangles::iterator                                 tri_iterator;

	// dimens num, cell vertex num, inner point id
	enum { D = 3, CVN = 8, inner_point_id = 6 };

	typedef t_float vertex_pos[D];
	typedef ulong   vertex_pos_i[D];
	typedef t_float cell_pos[CVN][D];

	// misc helper functions

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

	static Bbox vertex_pos2bbox(const vertex_pos& lo, const vertex_pos& hi) {
		return Bbox(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
	}

	static Point vertex_pos2point(const vertex_pos& p) {
		return Point(p[0], p[1], p[2]);
	}

	// add members nessessary to cover find intersection points and checking
	// whether point in inside cell
	template< class base_t >
	struct cell_data : public base_t {
		typedef std::vector< Tetrahedron > Tetrahedrons;

		// cell facets cover with triangles
		Triangles cover;
		// cell split into tetrahedrons
		Tetrahedrons split;

		// empty ctor for map
		cell_data() {}
		// std ctor
		cell_data(t_float *const cell) : base_t(cell) {}

		using base_t::cpos;
		using base_t::ss;

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

		// implement cell cover with triangles
		// to find exact intersection points
		void cell_tri_cover() {
			if(cover.size() > 0) return;
			cover.resize(12);
			const cell_pos& p = cpos();

			// facet 0-1-2-3
			// 3angle 0-1-3
			cover[0] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[1]),
				vertex_pos2point(p[3])
			);
			// 3angle 0-2-3
			cover[1] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[2]),
				vertex_pos2point(p[3])
			);
			// facet 0-1-4-5
			// 3angle 0-1-5
			cover[2] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[1]),
				vertex_pos2point(p[5])
			);
			// 3angle 0-4-5
			cover[3] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[4]),
				vertex_pos2point(p[5])
			);
			// facet 4-5-6-7
			// 3angle 4-5-7
			cover[4] = Triangle(
				vertex_pos2point(p[4]),
				vertex_pos2point(p[5]),
				vertex_pos2point(p[7])
			);
			// 3angle 4-6-7
			cover[5] = Triangle(
				vertex_pos2point(p[4]),
				vertex_pos2point(p[6]),
				vertex_pos2point(p[7])
			);
			// facet 2-3-6-7
			// 3angle 2-3-7
			cover[6] = Triangle(
				vertex_pos2point(p[2]),
				vertex_pos2point(p[3]),
				vertex_pos2point(p[7])
			);
			// 3angle 2-6-7
			cover[7] = Triangle(
				vertex_pos2point(p[2]),
				vertex_pos2point(p[6]),
				vertex_pos2point(p[7])
			);
			// facet 0-2-4-6
			// 3angle 0-2-6
			cover[8] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[2]),
				vertex_pos2point(p[6])
			);
			// 3angle 0-4-6
			cover[9] = Triangle(
				vertex_pos2point(p[0]),
				vertex_pos2point(p[4]),
				vertex_pos2point(p[6])
			);
			// facet 1-3-5-7
			// 3angle 1-3-7
			cover[10] = Triangle(
				vertex_pos2point(p[1]),
				vertex_pos2point(p[3]),
				vertex_pos2point(p[7])
			);
			// 3angle 1-5-7
			cover[11] = Triangle(
				vertex_pos2point(p[1]),
				vertex_pos2point(p[5]),
				vertex_pos2point(p[7])
			);
		}

		// implement testing whether point is inside cell
		bool contains(const Point& p) {
			// split cell into 5 tetrahedrons
			// and check whether point belongs to any of 'em
			if(!split.size()) {
				split.resize(5);
				// cell A B C D A' B' C' D'
				// ord  0 1 3 2 4  5  7  6
				// A A' B' D'
				// 0 4  5  6
				split[0] = Tetrahedron(ss(0), ss(4), ss(5), ss(6));
				// A B' B C
				// 0 5  1 3
				split[1] = Tetrahedron(ss(0), ss(5), ss(1), ss(3));
				// A C D D'
				// 0 3 2 6
				split[2] = Tetrahedron(ss(0), ss(3), ss(2), ss(6));
				// B' C' D' C
				// 5  7  6  3
				split[3] = Tetrahedron(ss(5), ss(7), ss(6), ss(3));
				// A C B' D'
				// 0 3 5  6
				split[4] = Tetrahedron(ss(0), ss(3), ss(5), ss(6));
			}

			// check each tetrahedron
			for(uint i = 0; i < split.size(); ++i) {
				if(!split[i].has_on_unbounded_side(p))
					return true;
			}
			return false;
		}
	};

	// well_data special part
	template< class base_t >
	struct well_data : public base_t {
		// ctors
		well_data() {}

		well_data(t_float *const segment, const well_data* prev = NULL) : base_t(segment) {}

		// MD access
		using base_t::W;
		t_float md() const {
			return W[3];
		}

		t_float md_finish() const {
			return W[7];
		}
	};

	typedef std::list< std::pair< Point, uint > > xpoints_list;

	// action taken on well & mesh boxes intersection
	template< class cell_data_t >
	static xpoints_list precise_intersection(cell_data_t& c, const Segment& well_seg) {
		// cover cell with triangles
		if(!c.cover.size())
			c.cell_tri_cover();

		xpoints_list res;
		// check that each triangle really intersects with given well segment
		//Segment s = w.segment();
		uint tri_count = 0;
		for(tri_iterator tri = c.cover.begin(), end = c.cover.end(); tri != end; ++tri, ++tri_count) {
			// test intersection
			//if(!CGAL::do_intersect(well_seg, *tri)) continue;
			// really do intersection
			Object xres = CGAL::intersection(well_seg, *tri);
			// in 99% of cases we should get a point of intersection
			if(const Point* xpoint = CGAL::object_cast< Point >(&xres))
				res.push_back(std::make_pair(*xpoint, tri_count >> 1));
			else if(const Segment* xseg = CGAL::object_cast< Segment >(&xres)) {
				// in rare 1% of segment lying on the facet, add begin and end of segment
				res.push_back(std::make_pair(xseg->source(), tri_count >> 1));
				res.push_back(std::make_pair(xseg->target(), tri_count >> 1));
			}
		}
		return res;
	}
};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_3D_SJLKT8NL */

