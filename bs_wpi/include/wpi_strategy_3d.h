/// @file wpi_strategy_3d.h
/// @brief 3D strategy for well path identification algorithms
/// @author uentity
/// @version 
/// @date 13.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGY_3D_SJLKT8NL
#define WPI_STRATEGY_3D_SJLKT8NL

#include "strategy_3d_common.h"
#include <CGAL/Bbox_3.h>

namespace blue_sky { namespace wpi {

template< template< uint > class strat_traits >
struct strategy_3d_ex : public strategy_3d_common< strat_traits > {
	// main typedefs
	typedef typename Kernel::Point_3                                     Point;
	typedef typename Kernel::Segment_3                                   Segment;
	typedef typename CGAL::Bbox_3                                        Bbox;
	typedef typename Kernel::Iso_cuboid_3                                Iso_bbox;

	// 3D specific typedefs
	typedef typename Kernel::Triangle_3                                  Triangle;
	typedef typename Kernel::Tetrahedron_3                               Tetrahedron;
	typedef typename std::vector<Triangle>                               Triangles;
	typedef typename Triangles::iterator                                 tri_iterator;

	// import common typedefs
	typedef strategy_3d_common< strat_traits >          strat_common;
	typedef typename strat_common::vertex_pos           vertex_pos;
	typedef typename strat_common::vertex_pos_i         vertex_pos_i;
	typedef typename strat_common::traits_t             traits_t;
	typedef typename strat_common::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_common::well_traj_iterator   well_traj_iterator;

	using strat_common::D;

	// misc helper functions
	static const char* name() {
		static std::string name_ = std::string("3D:") + traits_t::name();
		return name_.c_str();
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
	struct cell_data : public base_t, public strat_common::cell_data {
		typedef std::vector< Tetrahedron > Tetrahedrons;

		// cell facets cover with triangles
		//Triangles cover;
		// cell split into tetrahedrons
		//Tetrahedrons split;

		// empty ctor for map
		cell_data() {}
		// std ctor
		cell_data(const cell_vertex_iterator& cell) : base_t(cell) {}

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

		// implement testing whether point is inside cell
		bool contains(const Point& p) const {
			// split cell into 5 tetrahedrons
			// and check whether point belongs to any of 'em
			//if(!split.size()) {
			//	split.resize(5);
			const Tetrahedron split[] = {
				// cell A B C D A' B' C' D'
				// ord  0 1 3 2 4  5  7  6
				// A A' B' D'
				// 0 4  5  6
				Tetrahedron(ss(0), ss(4), ss(5), ss(6)),
				// A B' B C
				// 0 5  1 3
				Tetrahedron(ss(0), ss(5), ss(1), ss(3)),
				// A C D D'
				// 0 3 2 6
				Tetrahedron(ss(0), ss(3), ss(2), ss(6)),
				// B' C' D' C
				// 5  7  6  3
				Tetrahedron(ss(5), ss(7), ss(6), ss(3)),
				// A C B' D'
				// 0 3 5  6
				Tetrahedron(ss(0), ss(3), ss(5), ss(6))
			};

			// check each tetrahedron
			for(uint i = 0; i < 5; ++i) {
				if(!split[i].is_degenerate() && !split[i].has_on_unbounded_side(p))
					return true;
			}
			return false;
		}
	};

	typedef std::list< std::pair< Point, uint > > xpoints_list;

	// action taken on well & mesh boxes intersection
	template< class cell_data_t >
	static xpoints_list precise_intersection(const cell_data_t& c, const Segment& well_seg) {
		// cover cell with triangles
		// to find exact intersection points
		const Triangle cover[] = {
			// facet 0-1-2-3
			// 3angle 0-1-3
			Triangle(c.ss(0), c.ss(1), c.ss(3)),
			// 3angle 0-2-3
			Triangle(c.ss(0), c.ss(2), c.ss(3)),
			// facet 0-1-4-5
			// 3angle 0-1-5
			Triangle(c.ss(0), c.ss(1), c.ss(5)),
			// 3angle 0-4-5
			Triangle(c.ss(0), c.ss(4), c.ss(5)),
			// facet 4-5-6-7
			// 3angle 4-5-7
			Triangle(c.ss(4), c.ss(5), c.ss(7)),
			// 3angle 4-6-7
			Triangle(c.ss(4), c.ss(6), c.ss(7)),
			// facet 2-3-6-7
			// 3angle 2-3-7
			Triangle(c.ss(2), c.ss(3), c.ss(7)),
			// 3angle 2-6-7
			Triangle(c.ss(2), c.ss(6), c.ss(7)),
			// facet 0-2-4-6
			// 3angle 0-2-6
			Triangle(c.ss(0), c.ss(2), c.ss(6)),
			// 3angle 0-4-6
			Triangle(c.ss(0), c.ss(4), c.ss(6)),
			// facet 1-3-5-7
			// 3angle 1-3-7
			Triangle(c.ss(1), c.ss(3), c.ss(7)),
			// 3angle 1-5-7
			Triangle(c.ss(1), c.ss(5), c.ss(7))
		};

		xpoints_list res;
		// check that each triangle really intersects with given well segment
		for(uint tri_count = 0; tri_count < 12; ++tri_count) {
			try {
				// sanity
				// disabled: is_degenerate() checking is very expensive
				// so use try-catch approach to just skip invalid elements
				//if(tri->is_degenerate()) continue;
				//if(!CGAL::do_intersect(well_seg, *tri)) continue;

				// really do intersection
				Object xres = CGAL::intersection(well_seg, cover[tri_count]);
				// in 99% of cases we should get a point of intersection
				if(const Point* xpoint = CGAL::object_cast< Point >(&xres))
					res.push_back(std::make_pair(*xpoint, tri_count >> 1));
				else if(const Segment* xseg = CGAL::object_cast< Segment >(&xres)) {
					// in rare 1% of segment lying on the facet, add begin and end of segment
					// update: add only first point, two points confuse different algorithms
					res.push_back(std::make_pair(xseg->source(), tri_count >> 1));
					//res.push_back(std::make_pair(xseg->target(), tri_count >> 1));
				}
			}
			catch(const std::exception& e) {
				BOSERR << "WARNING! wpi::strategy_3d::precise_intersection: " << e.what() << bs_end;
			}
			catch(...) {
				BOSERR << "WARNING! wpi::strategy_3d::precise_intersection: unknown error" << bs_end;
				//throw;
			}
		}
		return res;
	}
};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_3D_SJLKT8NL */

