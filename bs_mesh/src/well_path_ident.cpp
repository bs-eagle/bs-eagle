/// @file well_path_ident.cpp
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 05.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"

#include "well_path_ident.h"
#include "wpi_iface.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo.h"

#include "export_python_wrapper.h"

// profiling
#include <google/profiler.h>

namespace blue_sky {
// alias
namespace bp = boost::python;

namespace {
template< class strat_t >
spv_uint where_is_points_impl(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
	//typedef wpi::strategy_3d strat_t;
	typedef wpi::pods< strat_t > pods_t;
	typedef wpi::mesh_tools< strat_t > mesh_tools_t;
	typedef wpi::algo< strat_t > algo;

	typedef wpi::ulong ulong;
	typedef typename pods_t::trimesh trimesh;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Point Point;

	//enum { D = strat_t::D };

	// convert coord & zcorn to tops
	trimesh M(nx, ny, coord, zcorn);
	//vertex_pos_i mesh_size;
	//spv_float tops = algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);

	// convert plain array of coords to array of Point objects
	// points are always specified in 3D
	ulong pnum = points->size() / 3;
	std::vector< Point > P(pnum);
	v_float::const_iterator p = points->begin();
	for(ulong i = 0; i < pnum; ++i) {
		P[i] = pods_t::rawptr2point(&*p);
		p += 3;
	}

	// real action
	const std::vector< ulong >& hit_idx = mesh_tools_t::where_is_point(M, P);

	// return result
	spv_uint res = BS_KERNEL.create_object(v_uint::bs_type());
	res->resize(hit_idx.size());
  if (hit_idx.size()) {
	  std::copy(hit_idx.begin(), hit_idx.end(), res->begin());
  }
	return res;
}

template< class strat_t >
t_uint where_is_point_impl(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
	//typedef wpi::strategy_3d strat_t;
	typedef wpi::pods< strat_t > pods_t;
	typedef wpi::mesh_tools< strat_t > mesh_tools_t;
	typedef wpi::algo< strat_t > algo;

	typedef typename pods_t::trimesh trimesh;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Point Point;

	//enum { D = strat_t::D };

	// convert coord & zcorn to tops
	trimesh M(nx, ny, coord, zcorn);
	//vertex_pos_i mesh_size;
	//spv_float tops = algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);

	// convert plain array of coords to Point object
	Point P = pods_t::rawptr2point(&*point->begin());

	// real action
	return mesh_tools_t::where_is_point(M, P);
}

} /* hidden implementation */

// specialization for 3D
spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	typedef wpi::algo< wpi::strategy_3d_ex< wpi::online_tops_traits > > wpi_algo_t;
	ProfilerStart("/home/uentity/my_projects/blue-sky.git/plugins/bs-eagle/examples/well_path_ident.prof");
	spv_float res = wpi_algo_t::well_path_ident_d< true >(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
	ProfilerStop();
	return res;
}

// specialization for 2D
spv_float well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return wpi::algo< wpi::strategy_2d >::well_path_ident_d< true >(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
}

// 3D
spv_uint where_is_points(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
	return where_is_points_impl< wpi::strategy_3d >(nx, ny, coord, zcorn, points);
}
// 2D
spv_uint where_is_points_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
	return where_is_points_impl< wpi::strategy_2d >(nx, ny, coord, zcorn, points);
}

// 3D
t_uint where_is_point(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
	return where_is_point_impl< wpi::strategy_3d >(nx, ny, coord, zcorn, point);
}
// 2D
t_uint where_is_point_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
	return where_is_point_impl< wpi::strategy_2d >(nx, ny, coord, zcorn, point);
}

/*-----------------------------------------------------------------
 * C++ iface
 *----------------------------------------------------------------*/
namespace wpi {

std::vector< well_hit_cell_3d > well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return algo< strategy_3d >::well_path_ident_d< false >(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
}

std::vector< well_hit_cell_2d > well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return algo< strategy_2d >::well_path_ident_d< false >(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
}

} /* wpi */

/*-----------------------------------------------------------------
 * Python bindings
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl, well_path_ident, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl_2d, well_path_ident_2d, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl_2d_old, well_path_ident_2d_old, 5, 6)

namespace python {

void py_export_wpi() {
	bp::def("well_path_ident", &well_path_ident, well_path_ident_overl());
	bp::def("well_path_ident_2d", &well_path_ident_2d, well_path_ident_overl_2d());
	bp::def("well_path_ident_2d_old", &well_path_ident_2d_old, well_path_ident_overl_2d_old());
	bp::def("where_is_point", &where_is_point);
	bp::def("where_is_points", &where_is_points);
	bp::def("where_is_point_2d", &where_is_point_2d);
	bp::def("where_is_points_2d", &where_is_points_2d);
}

}

}	// eof blue-sky namespace

