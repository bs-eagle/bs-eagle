/// @file well_path_ident_vtk.h
/// @brief VTK-related well path identification functions
/// @author uentity
/// @version 1.0
/// @date 29.11.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WELL_PATH_IDENT_VTK_34TI6RCA
#define WELL_PATH_IDENT_VTK_34TI6RCA

#include "conf.h"

namespace blue_sky {
/*-----------------------------------------------------------------
 * Well path intersection with given mesh
 * params:
 * nx, ny -- number of cells in X & Y directions
 * coord, zcorn -- COORD, ZCORN arrays that describes mesh
 * trimesh_backend -- alternate way to describe mesh using corresponding optimizations
 * include_well_nodes -- include well path nodes in resulting intersections set?
 * strat_traits -- strategy defining the way of internal mesh interpretation (abstract virtual layer)
 * hit_idx -- indexes of cells that contains well path nodes returned here (if != NULL)
 *----------------------------------------------------------------*/

BS_API_PLUGIN spv_ulong enum_border_facets_vtk(
	ulong nx, ulong ny, sp_obj trim_backend,
	spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
	const char* strat_traits = "online_tops",
	int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = 10000, int facet_filter = -1
);

BS_API_PLUGIN spv_ulong enum_border_facets_vtk(
	ulong nx, ulong ny, spv_float coord, spv_float zcorn,
	spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
	const char* strat_traits = "online_tops",
	int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = 10000, int facet_filter = -1
);

BS_API_PLUGIN spv_ulong enum_border_edges_vtk(
	ulong nx, ulong ny, sp_obj trim_backend,
	spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
	const char* strat_traits = "online_tops",
	int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = 10000, int facet_filter = -1
);

BS_API_PLUGIN spv_ulong enum_border_edges_vtk(
	ulong nx, ulong ny, spv_float coord, spv_float zcorn,
	spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
	const char* strat_traits = "online_tops",
	int slice_dim = -1, ulong slice_idx = 0,
	const ulong min_split_threshold = 10000, int facet_filter = -1
);

BS_API_PLUGIN sp_obj make_trimesh_backend(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	const char* strat_traits = "sgrid"
);

} // eof blue_sky namespace

#endif /* end of include guard: WELL_PATH_IDENT_VTK_34TI6RCA */

