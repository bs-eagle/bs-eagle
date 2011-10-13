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
#include "wpi_algo.h"

#include "export_python_wrapper.h"

namespace blue_sky {
// alias
namespace bp = boost::python;

// specialization for 3D
spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return wpi::algo< wpi::strategy_3d >::well_path_ident_d< true >(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
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

spv_uint where_is_points(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
	typedef wpi::strategy_3d strat_t;
	typedef wpi::pods< strat_t > pods_t;
	typedef wpi::mesh_tools< strat_t > mesh_tools_t;
	typedef wpi::algo< strat_t > algo;

	typedef wpi::ulong ulong;
	typedef pods_t::trimesh trimesh;
	typedef strat_t::vertex_pos_i vertex_pos_i;
	typedef strat_t::Point Point;

	enum { D = strat_t::D };

	// convert coord & zcorn to tops
	trimesh M;
	vertex_pos_i mesh_size;
	spv_float tops = algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);

	// convert plain array of coords to array of Point objects
	ulong pnum = points->size() / D;
	std::vector< Point > P(pnum);
	v_float::const_iterator p = points->begin();
	for(ulong i = 0; i < pnum; ++i) {
		P[i] = pods_t::rawptr2point(&*p);
		p += D;
	}

	// real action
	const std::vector< ulong >& hit_idx = mesh_tools_t::where_is_point(M, mesh_size, P);

	// return result
	spv_uint res = BS_KERNEL.create_object(v_uint::bs_type());
	res->resize(hit_idx.size());
	std::copy(hit_idx.begin(), hit_idx.end(), res->begin());
	return res;
}

t_uint where_is_point(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
	typedef wpi::strategy_3d strat_t;
	typedef wpi::pods< strat_t > pods_t;
	typedef wpi::mesh_tools< strat_t > mesh_tools_t;
	typedef wpi::algo< strat_t > algo;

	typedef pods_t::trimesh trimesh;
	typedef strat_t::vertex_pos_i vertex_pos_i;
	typedef strat_t::Point Point;

	enum { D = strat_t::D };

	// convert coord & zcorn to tops
	trimesh M;
	vertex_pos_i mesh_size;
	spv_float tops = algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);

	// convert plain array of coords to Point object
	Point P = pods_t::rawptr2point(&*point->begin());

	// real action
	return mesh_tools_t::where_is_point(M, mesh_size, P);
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
}

}

}	// eof blue-sky namespace

