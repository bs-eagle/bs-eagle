/// @file wpi_strategy_3d.h
/// @brief 3D strategy for well path identification algorithms
/// @author uentity
/// @version 
/// @date 13.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGY_2D_C5EHIKCY
#define WPI_STRATEGY_2D_C5EHIKCY

#include "strategy_2d_common.h"
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>

//#include "conf.h"

namespace blue_sky { namespace wpi {

template< template< uint > class strat_traits >
struct strategy_2d_ex : public strategy_2d_common< strat_traits > {
	// main typedefs
	typedef Kernel::Point_2                                     Point;
	typedef Kernel::Segment_2                                   Segment;
	typedef CGAL::Bbox_2                                        Bbox;
	typedef Kernel::Iso_rectangle_2                             Iso_bbox;

	// 2D specific typedefs
	typedef CGAL::Polygon_2< Kernel >                           Polygon_2;

	// import common typedefs
	typedef strategy_2d_common< strat_traits >          strat_common;
	typedef typename strat_common::vertex_pos           vertex_pos;
	typedef typename strat_common::vertex_pos_i         vertex_pos_i;
	typedef typename strat_common::traits_t             traits_t;
	typedef typename strat_common::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_common::well_traj_iterator   well_traj_iterator;

	using strat_common::D;

	// misc helper functions
	static const char* name() {
		static std::string name_ = std::string("2D:") + traits_t::name();
		return name_.c_str();
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
		static facet_vid_t& facet_vid(ulong facet) {
			static const facet_vid_t t[4] = {
				{0, 1}, {1, 3}, {3, 2}, {2, 0}
			};
			return t[facet];
		}
		static void facet_vid(ulong dim, ulong facet, facet_vid_t& res) {
			return facet_vid(facet_id(dim, facet, res));
		}

		Polygon_2 polygon() const {
			// 2D upper plane of cell
			Point points[] = {
				Point(V[0], V[1]),  Point(V[3], V[4]),
				Point(V[9], V[10]), Point(V[6], V[7])
			};
			return Polygon_2(points, points + 4);
		}

		bool contains(const Point& p) const {
			return !(polygon().has_on_unbounded_side(p));
		}
	};

	typedef std::list< std::pair< Point, uint > > xpoints_list;

	// action taken on well & mesh boxes intersection
	template< class cell_data_t >
	static xpoints_list precise_intersection(const cell_data_t& c, const Segment& well_seg) {
		// obtain well segment
		//const Segment& s = w.segment();
		// and cell polygon
		const Polygon_2& p = c.polygon();

		xpoints_list res;
		for(ulong i = 0; i < 4; ++i) {
			try {
				Object xres = CGAL::intersection(p.edge(i), well_seg);
				// in 99% of cases we should get a point of intersection
				if(const Point* xpoint = CGAL::object_cast< Point >(&xres))
					res.push_back(std::make_pair(*xpoint, i));
				else if(const Segment* xseg = CGAL::object_cast< Segment >(&xres)) {
					// in rare 1% of segment lying on the facet, add begin and end of segment
					// update: add only first point, two points confuse different algorithms
					res.push_back(std::make_pair(xseg->source(), i));
					//res.push_back(std::make_pair(xseg->target(), i));
				}
			}
			catch(std::exception& e) {
				BOSERR << "WARNING! wpi::strategy_2d::precise_intersection: " << e.what() << bs_end;
			}
			catch(...) {
				BOSERR << "WARNING! wpi::strategy_2d::precise_intersection: unknown error" << bs_end;
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
};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_2D_C5EHIKCY */

