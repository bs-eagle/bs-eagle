/// @file wpi_strategy_3d_bg.h
/// @brief 3D strategy for WPI algorithms using boost::geometry
/// @author uentity
/// @version 1.0
/// @date 17.02.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGY_3D_BG_HG2O3QNM
#define WPI_STRATEGY_3D_BG_HG2O3QNM

#include "strategy_3d_common.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/algorithms/intersection.hpp> 

namespace blue_sky { namespace wpi {

namespace bg = boost::geometry;

template< template< uint > class strat_traits >
struct strategy_3d_bg : public strategy_3d_common< strat_traits > {
	// import common typedefs
	typedef strategy_3d_common< strat_traits >          strat_common;
	typedef typename strat_common::vertex_pos           vertex_pos;
	typedef typename strat_common::vertex_pos_i         vertex_pos_i;
	typedef typename strat_common::traits_t             traits_t;
	typedef typename strat_common::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_common::well_traj_iterator   well_traj_iterator;

	using strat_common::D;

	// add some useful functions to Point class
	//struct Point : public bg::model::point< t_float, D, bg::cs::cartesian > {
	//	typedef bg::model::point< t_float, D, bg::cs::cartesian > base_t;

	//	// ctor from castable values (int, flaot, double, etc)
	//	template< class value_t >
	//	Point(const value_t& x, const value_t& y, const value_t& z)
	//		: base_t(x, y, z)
	//	{}

	//	// elements access
	//	const t_float& operator[](const ulong& i) const {
	//		
	//	}
	//};

	typedef bg::model::point< t_float, D, bg::cs::cartesian > Point;
	typedef bg::model::segment< t_float, D, bg::cs::cartesian > Segment;
	typedef bg::model::box< Point > Bbox;

	// misc helper functions
	static const char* name() {
		static std::string name_ = std::string("bg3D:") + traits_t::name();
		return name_.c_str();
	}

	static Bbox vertex_pos2bbox(const vertex_pos& lo, const vertex_pos& hi) {
		return Bbox(Point(lo[0], lo[1], lo[2]), Point(hi[0], hi[1], hi[2]));
	}

	static Point vertex_pos2point(const vertex_pos& p) {
		return Point(p[0], p[1], p[2]);
	}

	/*-----------------------------------------------------------------
	 * cell_data implementation
	 *----------------------------------------------------------------*/
	template< class base_t >
	struct cell_data : public base_t, public strat_common::cell_data {
		// empty ctor for map
		cell_data() {}
		// std ctor
		cell_data(const cell_vertex_iterator& cell) : base_t(cell) {}
	};

	/*-----------------------------------------------------------------
	 * action taken on well & mesh boxes intersection
	 *----------------------------------------------------------------*/
	// assume that cell_data_t is cell_data
	template< class cell_data_t >
	static xpoints_list precise_intersection(cell_data_t& c, const Segment& well_seg) {
	}
};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_STRATEGY_3D_BG_HG2O3QNM */

