/// @file well_path_ident_vtk.cpp
/// @brief VTK-related functions implementation
/// @author uentity
/// @version 1.0
/// @date 29.11.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

// VTK-related algos are only exported to Python now

#include "bs_common.h"
#include "well_path_ident_vtk.h"
#include "wpi_strategies.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo_vtk.h"

#include <string.h>

namespace blue_sky {

/*-----------------------------------------------------------------
 * implementation details
 *----------------------------------------------------------------*/
namespace { // hiddent implementation details

template< int prim_id >
spv_ulong enum_border_vtk_impl(
	t_ulong nx, t_ulong ny, sp_obj trim_backend, spv_int mask,
	spv_ulong cell_idx, spv_float points, const char* strat_traits,
	int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter
) {
	if(strcmp(strat_traits, "online_tops") == 0) {
		return wpi::algo_vtk< wpi::onlinett_3d >::enum_border_vtk(
			nx, ny, trim_backend, mask, cell_idx, points, Loki::Int2Type< prim_id >(),
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return wpi::algo_vtk< wpi::onlinett_bp_3d >::enum_border_vtk(
			nx, ny, trim_backend, mask, cell_idx, points, Loki::Int2Type< prim_id >(),
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0) {
		return wpi::algo_vtk< wpi::sgrid_3d >::enum_border_vtk(
			nx, ny, trim_backend, mask, cell_idx, points, Loki::Int2Type< prim_id >(),
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}
	else if(strcmp(strat_traits, "rgrid") == 0) {
		return wpi::algo_vtk< wpi::rgrid_3d >::enum_border_vtk(
			nx, ny, trim_backend, mask, cell_idx, points, Loki::Int2Type< prim_id >(),
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}
	else {
		// fallback to carray traits
		return wpi::algo_vtk< wpi::carray_3d >::enum_border_vtk(
			nx, ny, trim_backend, mask, cell_idx, points, Loki::Int2Type< prim_id >(),
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}
}

#ifdef BSPY_EXPORTING_PLUGIN
// alias
namespace bp = boost::python;

template< int prim_id >
bp::tuple py_enum_border_vtk_impl(t_ulong nx, t_ulong ny, sp_obj trim_backend,
	spv_int mask, const char* strat_traits = "online_tops",
	int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = MIN_SPLIT_THRESHOLD, const int facet_filter = -1)
{
	spv_ulong cell_idx = BS_KERNEL.create_object(v_ulong::bs_type());
	spv_float points = BS_KERNEL.create_object(v_float::bs_type());
	//ProfilerStart("/home/uentity/my_projects/blue-sky.git/gui/enum_border_facets_vtk.prof");
	spv_ulong res = enum_border_vtk_impl< prim_id >(
		nx, ny, trim_backend, mask, cell_idx, points, strat_traits,
		slice_dim, slice_idx, min_split_threshold, facet_filter
	);
	//ProfilerStop();
	return bp::make_tuple(res, cell_idx, points);
}

#endif
} //eof hidden namespace

/*-----------------------------------------------------------------
 * implementation of declared main functions
 *----------------------------------------------------------------*/
sp_obj make_trimesh_backend(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	const char* strat_traits)
{
	if(strcmp(strat_traits, "online_tops") == 0) {
		return wpi::pods< wpi::onlinett_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else if(strcmp(strat_traits, "online_tops_bufpool") == 0) {
		return wpi::pods< wpi::onlinett_bp_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else if(strcmp(strat_traits, "sgrid") == 0 || strcmp(strat_traits, "rgrid") == 0) {
		return wpi::pods< wpi::sgrid_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
	else {
		return wpi::pods< wpi::carray_3d >::trimesh::create_backend(
			nx, ny, coord, zcorn
		);
	}
}

spv_ulong enum_border_facets_vtk(
	t_ulong nx, t_ulong ny, sp_obj trim_backend, spv_int mask,
	spv_ulong cell_idx, spv_float points, const char* strat_traits,
	int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	return enum_border_vtk_impl< 0 >(
		nx, ny, trim_backend, mask, cell_idx, points, strat_traits, slice_dim, slice_idx,
		min_split_threshold, facet_filter
	);
}

spv_ulong enum_border_facets_vtk(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_int mask,
	spv_ulong cell_idx, spv_float points, const char* strat_traits,
	int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	return enum_border_vtk_impl< 0 >(
		nx, ny, make_trimesh_backend(nx, ny, coord, zcorn, strat_traits),
		mask, cell_idx, points, strat_traits, slice_dim, slice_idx,
		min_split_threshold, facet_filter
	);
}

spv_ulong enum_border_edges_vtk(
	t_ulong nx, t_ulong ny, sp_obj trim_backend, spv_int mask,
	spv_ulong cell_idx, spv_float points, const char* strat_traits,
	int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	return enum_border_vtk_impl< 1 >(
		nx, ny, trim_backend, mask, cell_idx, points, strat_traits, slice_dim, slice_idx,
		min_split_threshold, facet_filter
	);
}

spv_ulong enum_border_edges_vtk(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_int mask,
	spv_ulong cell_idx, spv_float points, const char* strat_traits,
	int slice_dim, ulong slice_idx,
	const ulong min_split_threshold, const int facet_filter)
{
	return enum_border_vtk_impl< 1 >(
		nx, ny, make_trimesh_backend(nx, ny, coord, zcorn, strat_traits),
		mask, cell_idx, points, strat_traits, slice_dim, slice_idx,
		min_split_threshold, facet_filter
	);
}

#ifdef BSPY_EXPORTING_PLUGIN
/*-----------------------------------------------------------------
 * Python bindings
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(make_trimbe_overl, make_trimesh_backend, 4, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(enumb_facets_overl, py_enum_border_vtk_impl< 0 >, 4, 9)
BOOST_PYTHON_FUNCTION_OVERLOADS(enumb_edges_overl, py_enum_border_vtk_impl< 1 >, 4, 9)

namespace python {

void py_export_wpi_vtk() {
	// export trimesh backend creation fcn
	bp::def("make_trimesh_backend", &make_trimesh_backend, make_trimbe_overl());

	bp::def("enum_border_facets_vtk", &py_enum_border_vtk_impl< 0 >, enumb_facets_overl());
	bp::def("enum_border_edges_vtk", &py_enum_border_vtk_impl< 1 >, enumb_edges_overl());
}

} // eof python namespace
#endif

} // eof blue_sky namespace

