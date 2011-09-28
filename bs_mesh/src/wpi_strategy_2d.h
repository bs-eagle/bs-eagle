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

struct wpi_strategy_2d {
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
		// empty ctor for map
		cell_data() {}
		// std ctor
		cell_data(t_float *const cell) : base_t(cell) {}

		using base_t::V;

		Polygon_2 polygon() const {
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
		well_data(t_float *const segment, const well_data* prev = NULL)
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

};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_2D_C5EHIKCY */

