/// @file well_path_ident.h
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 11.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
#include "conf.h"

namespace blue_sky {

spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true);

spv_float well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true);

} 	// eof blue_sky

