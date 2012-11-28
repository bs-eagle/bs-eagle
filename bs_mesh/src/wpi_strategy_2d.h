/// @file wpi_strategy_3d.h
/// @brief 3D strategy for well path identification algorithms
/// @author uentity
/// @version 
/// @date 13.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGY_2D_C5EHIKCY
#define WPI_STRATEGY_2D_C5EHIKCY

#include "wpi_common.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Object.h>
//#include <CGAL/box_intersection_d.h>
//#include <CGAL/intersections.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>

//#include "conf.h"

namespace blue_sky { namespace wpi {

template< class strat_traits >
struct strategy_2d_ex {
	// main typedefs
	typedef Kernel::Point_2                                     Point;
	typedef Kernel::Segment_2                                   Segment;
	typedef CGAL::Bbox_2                                        Bbox;
	typedef Kernel::Iso_rectangle_2                             Iso_bbox;

	// 2D specific typedefs
	typedef CGAL::Polygon_2< Kernel >                           Polygon_2;

	// dimens num, inner point id
	enum { D = 2, inner_point_id = 4 };

	typedef t_float vertex_pos[D];
	typedef ulong   vertex_pos_i[D];

	// iterator over source arrays come from traits
	typedef strat_traits traits_t;
	typedef typename traits_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename traits_t::well_traj_iterator   well_traj_iterator;

	// misc helper functions

	// X-Y-Z order!
	static void decode_cell_id(ulong id, vertex_pos_i& res, const vertex_pos_i& m_size) {
		//vertex_pos_i res;
		res[1] = id / m_size[0];
		res[0] = id - m_size[0]*res[1];
	}

	static ulong encode_cell_id(const vertex_pos_i& p, const vertex_pos_i& m_size) {
		return p[0] + m_size[0] * p[1];
	}

	static Bbox vertex_pos2bbox(const vertex_pos& lo, const vertex_pos& hi) {
		return Bbox(lo[0], lo[1], hi[0], hi[1]);
	}

	static Point vertex_pos2point(const vertex_pos& p) {
		return Point(p[0], p[1]);
	}

	// cell_data specific in 2D
	template< class base_t >
	struct cell_data : public base_t {
		enum { n_vertex = 4 };
		enum { n_facets = 4 };
		enum { n_facet_vertex = 2 };
		enum { n_edges = 4 };
		// array to hold facet vertices IDs
		typedef ulong facet_vid_t[n_facet_vertex];

		// empty ctor for map
		cell_data() {}
		// std ctor
		cell_data(const cell_vertex_iterator& cell) : base_t(cell) {}

		using base_t::V;

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
		/* working with XY upper plane of cell (0-1-3-2)
		*  facets (and edges) layout (nodes - facet id)
		*  0-1 - 0
		*  1-3 - 1
		*  3-2 - 2
		*  2-0 - 3
		*  inside plane - 4
		*/

		static ulong facet_id(ulong dim, ulong facet) {
			if(dim == 0) {
				// X
				return 3 - facet*2;
			}
			else if(dim == 1) {
				// Y
				return facet*2;
			}
			return ulong(-1);
		}
		// obtain IDs of given facet vertices
		static void facet_vid(ulong facet, facet_vid_t& res) {
			switch(facet) {
				case 0 : {
					facet_vid_t t = {0, 1};
					ca_assign(res, t); }
					break;
				case 1 : {
					facet_vid_t t = {1, 3};
					ca_assign(res, t); }
					break;
				case 2 : {
					facet_vid_t t = {3, 2};
					ca_assign(res, t); }
					break;
				case 3 : {
					facet_vid_t t = {2, 0};
					ca_assign(res, t); }
					break;
				default : {
					facet_vid_t t = {ulong(-1), ulong(-1)};
					ca_assign(res, t); }
			}
		}
		static void facet_vid(ulong dim, ulong facet, facet_vid_t& res) {
			return facet_vid(facet_id(dim, facet, res));
		}

		Polygon_2 polygon() {
			// 2D upper plane of cell
			Point points[] = {
				Point(V[0], V[1]),  Point(V[3], V[4]),
				Point(V[9], V[10]), Point(V[6], V[7])
			};
			return Polygon_2(points, points + 4);
		}

		bool contains(const Point& p) {
			return !(polygon().has_on_unbounded_side(p));
		}
	};

	// well_data contains specifid MD member for
	// segment distances in 2D

	template< class base_t >
	struct well_data : public base_t {
		double md_;

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

	typedef std::list< std::pair< Point, uint > > xpoints_list;

	// action taken on well & mesh boxes intersection
	template< class cell_data_t >
	static xpoints_list precise_intersection(cell_data_t& c, const Segment& well_seg) {
		// obtain well segment
		//const Segment& s = w.segment();
		// and cell polygon
		const Polygon_2& p = c.polygon();

		xpoints_list res;
		for(ulong i = 0; i < 4; ++i) {
			Object xres = CGAL::intersection(p.edge(i), well_seg);
			// in 99% of cases we should get a point of intersection
			if(const Point* xpoint = CGAL::object_cast< Point >(&xres))
				res.push_back(std::make_pair(*xpoint, i));
			else if(const Segment* xseg = CGAL::object_cast< Segment >(&xres)) {
				// in rare 1% of segment lying on the facet, add begin and end of segment
				res.push_back(std::make_pair(xseg->source(), i));
				res.push_back(std::make_pair(xseg->target(), i));
			}
		}
		return res;
	}

	// fast chack if bbox intersect with segment
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
		}

		if(t_min[0] > t_max[1] || t_min[1] > t_max[0])
			return false;
		return true;
	}

	static bool bbox_segment_x(const Bbox& b, const Segment& s) {
		return bbox_segment_x(
			Point(b.min(0), b.min(1)),
			Point(b.max(0), b.max(1)),
			s.min(), s.max()
		);
	}

	// look into wpi_strategy_3d.h for description of this function
	typedef int bbox_bnd_offs[2][D];
	static const bbox_bnd_offs& bbox_boundary_offs(const uint dim, const uint bnd_id) {
		static const bbox_bnd_offs t[4] = {
			// X
			{ {0, 1}, {0, -1} },
			{ {0, 0}, {0, -1} },
			// Y
			{ {0, 0}, {-1, 0} },
			{ {0, 0}, { 0, 0} },
		};
		return t[dim*2 + bnd_id];
	}
};

// shortcoming typedef
typedef strategy_2d_ex< online_tops_traits > strategy_2d;
//typedef strategy_2d_ex< carray_traits > strategy_2d;

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_2D_C5EHIKCY */

