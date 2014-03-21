/// @file well_path_ident.h
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 11.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WELL_PATH_IDENT_F1N56F3I
#define WELL_PATH_IDENT_F1N56F3I

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

BS_API_PLUGIN spv_float well_path_ident(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
	const ulong min_split_threshold = 0, spv_ulong hit_idx = NULL
);
BS_API_PLUGIN spv_float well_path_ident(
	t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
	spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
	const ulong min_split_threshold = 0, spv_ulong hit_idx = NULL
);

BS_API_PLUGIN spv_float well_path_ident_2d(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
	const ulong min_split_threshold = 0, spv_ulong hit_idx = NULL
);
BS_API_PLUGIN spv_float well_path_ident_2d(
	t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
	spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
	const ulong min_split_threshold = 0, spv_ulong hit_idx = NULL
);

BS_API_PLUGIN spv_ulong where_is_points(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
	const char* strat_traits = "online_tops"
);
BS_API_PLUGIN spv_ulong where_is_points_2d(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
	const char* strat_traits = "online_tops"
);

BS_API_PLUGIN t_ulong where_is_point(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
	const char* strat_traits = "online_tops"
);
BS_API_PLUGIN t_ulong where_is_point_2d(
	t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
	const char* strat_traits = "online_tops"
);

} 	// eof blue_sky

#endif /* end of include guard: WELL_PATH_IDENT_F1N56F3I */

