/// @file wpi_iface.h
/// @brief C++ interface for acessing well path ident algos
/// @author uentity
/// @version 
/// @date 29.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_IFACE_NSM6SQCA
#define WPI_IFACE_NSM6SQCA

#include "wpi_algo_pod.h"
#include "wpi_strategy_3d.h"
#include "wpi_strategy_2d.h"

namespace blue_sky { namespace wpi {
// handy typedefs
typedef pods< strategy_3d >::well_hit_cell well_hit_cell_3d;
typedef pods< strategy_2d >::well_hit_cell well_hit_cell_2d;

std::vector< well_hit_cell_3d > well_path_ident(
	t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true);

std::vector< well_hit_cell_2d > well_path_ident_2d(
	t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes = true);

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_IFACE_NSM6SQCA */

