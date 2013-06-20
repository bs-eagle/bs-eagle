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
#include "wpi_algo_vtk.h"

#include "export_python_wrapper.h"

#include <string.h>

// profiling
//#include <google/profiler.h>

namespace blue_sky {
// alias
namespace bp = boost::python;

namespace {

typedef wpi::strategy_3d_ex< wpi::carray_traits >              carray_3d;
typedef wpi::strategy_3d_ex< wpi::online_tops_traits >         onlinett_3d;
typedef wpi::strategy_3d_ex< wpi::online_tops_traits_bufpool > onlinett_bp_3d;
typedef wpi::strategy_3d_ex< wpi::sgrid_traits >               sgrid_3d;

typedef wpi::strategy_2d_ex< wpi::carray_traits >              carray_2d;
typedef wpi::strategy_2d_ex< wpi::online_tops_traits >         onlinett_2d;
typedef wpi::strategy_2d_ex< wpi::online_tops_traits_bufpool > onlinett_bp_2d;
typedef wpi::strategy_2d_ex< wpi::sgrid_traits >               sgrid_2d;

template< class strat_t >
spv_ulong where_is_points_impl(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
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
	spv_ulong res = BS_KERNEL.create_object(v_ulong::bs_type());
	res->resize(hit_idx.size());
	if (hit_idx.size()) {
		std::copy(hit_idx.begin(), hit_idx.end(), res->begin());
	}
	return res;
}

template< class strat_t >
t_ulong where_is_point_impl(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
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

/*-----------------------------------------------------------------
 * VTK-related algorithms wrappers for Python
 *----------------------------------------------------------------*/
#ifdef BSPY_EXPORTING_PLUGIN

//typedef wpi::strategy_3d_ex< wpi::sgrid_traits > vtk_strat_t;
//typedef wpi::algo_vtk< vtk_strat_t > wpi_algo_vtk;

template< class vtk_strat_t >
bp::tuple enum_border_facets_vtk_impl(t_long nx, t_long ny, sp_obj trim_backend,
	spv_int mask, int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	typedef wpi::algo_vtk< vtk_strat_t > wpi_algo_vtk;

	spv_long cell_idx = BS_KERNEL.create_object(v_long::bs_type());
	spv_float points = BS_KERNEL.create_object(v_float::bs_type());
	//ProfilerStart("/home/uentity/my_projects/blue-sky.git/gui/enum_border_facets_vtk.prof");
	spv_long res = wpi_algo_vtk::enum_border_facets_vtk(nx, ny, trim_backend, mask, cell_idx,
		points, slice_dim, slice_idx, min_split_threshold, facet_filter);
	//ProfilerStop();
	return bp::make_tuple(res, cell_idx, points);
}

template< class vtk_strat_t >
bp::tuple enum_border_edges_vtk_impl(t_long nx, t_long ny, sp_obj trim_backend,
	spv_int mask, int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	typedef wpi::algo_vtk< vtk_strat_t > wpi_algo_vtk;

	spv_long cell_idx = BS_KERNEL.create_object(v_long::bs_type());
	spv_float points = BS_KERNEL.create_object(v_float::bs_type());
	//ProfilerStart("/home/uentity/my_projects/blue-sky.git/gui/enum_border_edges_vtk.prof");
	//HeapProfilerStart("/home/uentity/my_projects/blue-sky.git/gui/enum_border_edges_vtk");
	spv_long res = wpi_algo_vtk::enum_border_edges_vtk(nx, ny, trim_backend, mask, cell_idx,
		points, slice_dim, slice_idx, min_split_threshold, facet_filter);
	//ProfilerStop();
	return bp::make_tuple(res, cell_idx, points);
}

#endif

} /* hidden implementation */

// specialization for 3D
spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes, const char* strat_traits)
{
	//ProfilerStart("/home/uentity/my_projects/blue-sky.git/plugins/bs-eagle/examples/well_path_ident.prof");
	if(strcmp(strat_traits, "online_tops") == 0) {
		return  wpi::algo< onlinett_3d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return  wpi::algo< onlinett_bp_3d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return  wpi::algo< sgrid_3d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else {
		// fallback to carray traits
		return  wpi::algo< carray_3d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	//ProfilerStop();
}

spv_float well_path_ident(t_long nx, t_long ny, sp_obj trimesh_backend,
	spv_float well_info, bool include_well_nodes, const char* strat_traits)
{
	//ProfilerStart("/home/uentity/my_projects/blue-sky.git/plugins/bs-eagle/examples/well_path_ident.prof");
	if(strcmp(strat_traits, "online_tops") == 0) {
		return  wpi::algo< onlinett_3d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return  wpi::algo< onlinett_bp_3d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return  wpi::algo< sgrid_3d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else {
		// fallback to carray traits
		return  wpi::algo< carray_3d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	//ProfilerStop();
}

// specialization for 2D
spv_float well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes, const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return  wpi::algo< onlinett_2d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return  wpi::algo< onlinett_bp_2d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return  wpi::algo< sgrid_2d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
	else {
		// fallback to carray traits
		return  wpi::algo< carray_2d >::well_path_ident_d< true >(
			nx, ny, coord, zcorn, well_info, include_well_nodes
		);
	}
}

spv_float well_path_ident_2d(t_long nx, t_long ny, sp_obj trimesh_backend,
	spv_float well_info, bool include_well_nodes, const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return  wpi::algo< onlinett_2d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return  wpi::algo< onlinett_bp_2d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return  wpi::algo< sgrid_2d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
	else {
		// fallback to carray traits
		return  wpi::algo< carray_2d >::well_path_ident_d< true >(
			nx, ny, trimesh_backend, well_info, include_well_nodes
		);
	}
}

// 3D
spv_ulong where_is_points(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points,
	const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return where_is_points_impl< onlinett_3d >(nx, ny, coord, zcorn, points);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return where_is_points_impl< onlinett_bp_3d >(nx, ny, coord, zcorn, points);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return where_is_points_impl< sgrid_3d >(nx, ny, coord, zcorn, points);
	}
	else {
		// fallback to carray traits
		return where_is_points_impl< carray_3d >(nx, ny, coord, zcorn, points);
	}
}
// 2D
spv_ulong where_is_points_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points,
	const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return where_is_points_impl< onlinett_2d >(nx, ny, coord, zcorn, points);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return where_is_points_impl< onlinett_bp_2d >(nx, ny, coord, zcorn, points);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return where_is_points_impl< sgrid_2d >(nx, ny, coord, zcorn, points);
	}
	else {
		// fallback to carray traits
		return where_is_points_impl< carray_2d >(nx, ny, coord, zcorn, points);
	}
}

// 3D
t_ulong where_is_point(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point,
	const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return where_is_point_impl< onlinett_3d >(nx, ny, coord, zcorn, point);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return where_is_point_impl< onlinett_bp_3d >(nx, ny, coord, zcorn, point);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return where_is_point_impl< sgrid_3d >(nx, ny, coord, zcorn, point);
	}
	else {
		// fallback to carray traits
		return where_is_point_impl< carray_3d >(nx, ny, coord, zcorn, point);
	}
}
// 2D
t_ulong where_is_point_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point,
	const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return where_is_point_impl< onlinett_2d >(nx, ny, coord, zcorn, point);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return where_is_point_impl< onlinett_bp_2d >(nx, ny, coord, zcorn, point);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return where_is_point_impl< sgrid_2d >(nx, ny, coord, zcorn, point);
	}
	else {
		// fallback to carray traits
		return where_is_point_impl< carray_2d >(nx, ny, coord, zcorn, point);
	}
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

#ifdef BSPY_EXPORTING_PLUGIN

bp::tuple enum_border_facets_vtk(t_long nx, t_long ny, sp_obj trim_backend,
	spv_int mask, int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = MIN_SPLIT_THRESHOLD, const int facet_filter = -1,
	const char* strat_traits = "sgrid")
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return enum_border_facets_vtk_impl< onlinett_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return enum_border_facets_vtk_impl< onlinett_bp_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return enum_border_facets_vtk_impl< sgrid_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else {
		// fallback to carray traits
		return enum_border_facets_vtk_impl< carray_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
}

bp::tuple enum_border_edges_vtk(t_long nx, t_long ny, sp_obj trim_backend,
	spv_int mask, int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = MIN_SPLIT_THRESHOLD, const int facet_filter = -1,
	const char* strat_traits = "sgrid")
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return enum_border_edges_vtk_impl< onlinett_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return enum_border_edges_vtk_impl< onlinett_bp_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return enum_border_edges_vtk_impl< sgrid_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
	else {
		// fallback to carray traits
		return enum_border_edges_vtk_impl< carray_3d >(
			nx, ny, trim_backend, mask, slice_dim, slice_idx,
			min_split_threshold, facet_filter
		);
	}
}

sp_obj make_trimesh_backend(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	const char* strat_traits = "sgrid")
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return wpi::pods< onlinett_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return wpi::pods< onlinett_bp_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return wpi::pods< sgrid_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else {
		return wpi::pods< carray_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
}

/*-----------------------------------------------------------------
 * Python bindings
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl1, ::blue_sky::well_path_ident, 5, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl2, ::blue_sky::well_path_ident, 4, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_2d_overl1, ::blue_sky::well_path_ident_2d, 5, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_2d_overl2, ::blue_sky::well_path_ident_2d, 4, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(where_is_point_overl, where_is_point, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(where_is_points_overl, where_is_points, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(where_is_point_2d_overl, where_is_point_2d, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(where_is_points_2d_overl, where_is_points_2d, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(enumb_facets_overl, enum_border_facets_vtk, 4, 9)
BOOST_PYTHON_FUNCTION_OVERLOADS(enumb_edges_overl, enum_border_edges_vtk, 4, 9)
BOOST_PYTHON_FUNCTION_OVERLOADS(make_trimbe_overl, make_trimesh_backend, 4, 5)

namespace python {

void py_export_wpi() {
	spv_float (*wpi_3d_1)(t_long, t_long, spv_float, spv_float, spv_float, bool, const char*) =
		&::blue_sky::well_path_ident;
	spv_float (*wpi_3d_2)(t_long, t_long, sp_obj, spv_float, bool, const char*) =
		&::blue_sky::well_path_ident;
	spv_float (*wpi_2d_1)(t_long, t_long, spv_float, spv_float, spv_float, bool, const char*) =
		&::blue_sky::well_path_ident_2d;
	spv_float (*wpi_2d_2)(t_long, t_long, sp_obj, spv_float, bool, const char*) =
		&::blue_sky::well_path_ident_2d;

	bp::def("well_path_ident", wpi_3d_1, well_path_ident_overl1());
	bp::def("well_path_ident", wpi_3d_2, well_path_ident_overl2());
	bp::def("well_path_ident_2d", wpi_2d_1, well_path_ident_2d_overl1());
	bp::def("well_path_ident_2d", wpi_2d_2, well_path_ident_2d_overl2());
	bp::def("where_is_point", &where_is_point, where_is_point_overl());
	bp::def("where_is_points", &where_is_points, where_is_points_overl());
	bp::def("where_is_point_2d", &where_is_point_2d, where_is_point_2d_overl());
	bp::def("where_is_points_2d", &where_is_points_2d, where_is_points_2d_overl());

	// export trimesh backend creation fcn
	bp::def("make_trimesh_backend", &make_trimesh_backend, make_trimbe_overl());

	bp::def("enum_border_facets_vtk", &enum_border_facets_vtk, enumb_facets_overl());
	bp::def("enum_border_edges_vtk", &enum_border_edges_vtk, enumb_edges_overl());
}

}	// eof python namespace
#endif

}	// eof blue-sky namespace

