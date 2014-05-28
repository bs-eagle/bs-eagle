/// @file wpi_iface.h
/// @brief C++ interface for acessing well path ident algos
/// @author uentity
/// @version 1.0
/// @date 29.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_IFACE_NSM6SQCA
#define WPI_IFACE_NSM6SQCA

#include "bs_kernel.h"
#include "conf.h"

namespace blue_sky {

// interface to allow access WPI functions without linking to library
class BS_API_PLUGIN wpi_iface : public objbase {
public:
	/*-----------------------------------------------------------------
	 * Generic well path ident functions
	 *----------------------------------------------------------------*/
	virtual spv_float well_path_ident(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "sgrid",
		spv_ulong hit_idx = NULL
	) = 0;
	virtual spv_float well_path_ident(
		t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
		spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "sgrid",
		spv_ulong hit_idx = NULL
	) = 0;

	virtual spv_float well_path_ident_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
		spv_ulong hit_idx = NULL
	) = 0;
	virtual spv_float well_path_ident_2d(
		t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
		spv_float well_info, bool include_well_nodes = true, const char* strat_traits = "online_tops",
		spv_ulong hit_idx = NULL
	) = 0;

	virtual spv_ulong where_is_points(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
		const char* strat_traits = "online_tops"
	) = 0;
	virtual spv_ulong where_is_points_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
		const char* strat_traits = "online_tops"
	) = 0;

	virtual t_ulong where_is_point(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
		const char* strat_traits = "online_tops"
	) = 0;
	virtual t_ulong where_is_point_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
		const char* strat_traits = "online_tops"
	) = 0;

	/*-----------------------------------------------------------------
	 * VTK-related WPI functions
	 *----------------------------------------------------------------*/
	virtual spv_ulong enum_border_facets_vtk(
		ulong nx, ulong ny, sp_obj trim_backend,
		spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
		const char* strat_traits = "online_tops",
		int slice_dim = -1, ulong slice_idx = 0,
		const ulong min_split_threshold = 10000, int facet_filter = -1
	) = 0;

	virtual spv_ulong enum_border_facets_vtk(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
		const char* strat_traits = "online_tops",
		int slice_dim = -1, ulong slice_idx = 0,
		const ulong min_split_threshold = 10000, int facet_filter = -1
	) = 0;

	virtual spv_ulong enum_border_edges_vtk(
		ulong nx, ulong ny, sp_obj trim_backend,
		spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
		const char* strat_traits = "online_tops",
		int slice_dim = -1, ulong slice_idx = 0,
		const ulong min_split_threshold = 10000, int facet_filter = -1
	) = 0;

	virtual spv_ulong enum_border_edges_vtk(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_ulong cell_idx, spv_float points, spv_int mask = NULL,
		const char* strat_traits = "online_tops",
		int slice_dim = -1, ulong slice_idx = 0,
		const ulong min_split_threshold = 10000, int facet_filter = -1
	) = 0;

	virtual sp_obj make_trimesh_backend(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		const char* strat_traits = "sgrid"
	) = 0;

	BLUE_SKY_TYPE_DECL_IFACE(wpi_iface)
};
typedef smart_ptr< wpi_iface > sp_iwpi;

} /* eof blue_sky */

#endif /* end of include guard: WPI_IFACE_NSM6SQCA */

