/// @file wpi_iface.cpp
/// @brief Implementation of wpi_iface
/// @author uentity
/// @version 1.0
/// @date 12.12.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "wpi_iface.h"
#include "well_path_ident.h"
#include "well_path_ident_vtk.h"
#include "wpi_algo.h"

namespace blue_sky {

namespace {
/*-----------------------------------------------------------------
 * Forward interface calls to corresponding functions
 *----------------------------------------------------------------*/
class wpi_object : public wpi_iface {
public:
	wpi_object(bs_type_ctor_param = NULL) {}

	/*-----------------------------------------------------------------
	 * Generic well path ident functions
	 *----------------------------------------------------------------*/
	spv_float well_path_ident(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes, const char* strat_traits,
		spv_ulong hit_idx
	) {
		return blue_sky::well_path_ident(
			nx, ny, coord, zcorn, well_info, include_well_nodes, strat_traits, hit_idx
		);
	}

	spv_float well_path_ident(
		t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
		spv_float well_info, bool include_well_nodes, const char* strat_traits,
		spv_ulong hit_idx
	) {
		return blue_sky::well_path_ident(
			nx, ny, trimesh_backend, well_info, include_well_nodes, strat_traits, hit_idx
		);
	}

	spv_float well_path_ident_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes, const char* strat_traits,
		spv_ulong hit_idx
	) {
		return blue_sky::well_path_ident_2d(
			nx, ny, coord, zcorn, well_info, include_well_nodes, strat_traits, hit_idx
		);
	}

	spv_float well_path_ident_2d(
		t_ulong nx, t_ulong ny, sp_obj trimesh_backend,
		spv_float well_info, bool include_well_nodes, const char* strat_traits,
		spv_ulong hit_idx
	) {
		return blue_sky::well_path_ident_2d(
			nx, ny, trimesh_backend, well_info, include_well_nodes, strat_traits, hit_idx
		);
	}

	spv_ulong where_is_points(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
		const char* strat_traits
	) {
		return blue_sky::where_is_points(
			nx, ny, coord, zcorn, points, strat_traits
		);
	}

	spv_ulong where_is_points_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float points,
		const char* strat_traits
	) {
		return blue_sky::where_is_points_2d(
			nx, ny, coord, zcorn, points, strat_traits
		);
	}

	t_ulong where_is_point(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
		const char* strat_traits
	) {
		return blue_sky::where_is_point(
			nx, ny, coord, zcorn, point, strat_traits
		);
	}

	t_ulong where_is_point_2d(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn, spv_float point,
		const char* strat_traits
	) {
		return blue_sky::where_is_point_2d(
			nx, ny, coord, zcorn, point, strat_traits
		);
	}

	/*-----------------------------------------------------------------
	 * VTK-related WPI functions
	 *----------------------------------------------------------------*/
	virtual spv_ulong enum_border_facets_vtk(
		ulong nx, ulong ny, sp_obj trim_backend,
		spv_ulong cell_idx, spv_float points, spv_int mask,
		const char* strat_traits,
		int slice_dim, ulong slice_idx,
		const ulong min_split_threshold, int facet_filter
	) {
		return blue_sky::enum_border_facets_vtk(
			nx, ny, trim_backend, cell_idx, points, mask, strat_traits,
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}

	virtual spv_ulong enum_border_facets_vtk(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_ulong cell_idx, spv_float points, spv_int mask,
		const char* strat_traits,
		int slice_dim, ulong slice_idx,
		const ulong min_split_threshold, int facet_filter
	) {
		return blue_sky::enum_border_facets_vtk(
			nx, ny, coord, zcorn, cell_idx, points, mask, strat_traits,
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}

	virtual spv_ulong enum_border_edges_vtk(
		ulong nx, ulong ny, sp_obj trim_backend,
		spv_ulong cell_idx, spv_float points, spv_int mask,
		const char* strat_traits,
		int slice_dim, ulong slice_idx,
		const ulong min_split_threshold, int facet_filter
	) {
		return blue_sky::enum_border_edges_vtk(
			nx, ny, trim_backend, cell_idx, points, mask, strat_traits,
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}

	virtual spv_ulong enum_border_edges_vtk(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_ulong cell_idx, spv_float points, spv_int mask,
		const char* strat_traits,
		int slice_dim, ulong slice_idx,
		const ulong min_split_threshold, int facet_filter
	) {
		return blue_sky::enum_border_edges_vtk(
			nx, ny, coord, zcorn, cell_idx, points, mask, strat_traits,
			slice_dim, slice_idx, min_split_threshold, facet_filter
		);
	}

	sp_obj make_trimesh_backend(
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		const char* strat_traits
	) {
		return blue_sky::make_trimesh_backend(
			nx, ny, coord, zcorn, strat_traits
		);
	}

	BS_RESOLVE_TYPE_IMPL_MEM
};

} // eof hidden namespace

/*-----------------------------------------------------------------
 * BS-related implementation stuff
 *----------------------------------------------------------------*/
BLUE_SKY_TYPE_IMPL(wpi_iface, objbase, "wpi_iface", "Interface to WPI functions", "")
BLUE_SKY_TYPE_STD_CREATE_IFACE(wpi_iface, wpi_object)
BLUE_SKY_TYPE_STD_COPY_IFACE(wpi_iface, wpi_object)

bool register_wpi_iface(const plugin_descriptor& pd) {
	return BS_KERNEL.register_type(pd, wpi_iface::bs_type());
}

} /* namespace blue_sky */

