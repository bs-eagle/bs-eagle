/// @file czt_iface.cpp
/// @brief BS interface for CZT functions implementation
/// @author uentity
/// @version 1.0
/// @date 12.12.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "czt_iface.h"

namespace blue_sky {
//namespace czt = coord_zcorn_tools;

namespace {

class czt_object : public czt_iface {
public:
	// empty ctor
	czt_object(bs_type_ctor_param = NULL) {}

	/*-----------------------------------------------------------------
	* Main processing functions
	*----------------------------------------------------------------*/
	spfp_storarr_t gen_coord(
		int_t nx, int_t ny, spfp_storarr_t dx, spfp_storarr_t dy, fp_t x0, fp_t y0
	) {
		return czt::gen_coord(
			nx, ny, dx, dy, x0, y0
		);
	}

	coord_zcorn_pair gen_coord_zcorn(
		int_t nx, int_t ny, int_t nz,
		spv_float dx, spv_float dy, spv_float dz,
		fp_stor_t x0, fp_stor_t y0, fp_stor_t z0
	) {
		return czt::gen_coord_zcorn(
			nx, ny, nz, dx, dy, dz, x0, y0, z0
		);
	}

	spfp_storarr_t gen_coord2(spfp_storarr_t x, spfp_storarr_t y) {
		return czt::gen_coord2(x, y);
	}

	coord_zcorn_pair refine_mesh_deltas(
		int_t& nx, int_t& ny, spfp_storarr_t coord,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx
	) {
		return czt::refine_mesh_deltas(
			nx, ny, coord, points, cell_merge_thresh, band_thresh, hit_idx
		);
	}

	coord_zcorn_pair refine_mesh(
		int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx
	) {
		return czt::refine_mesh(
			nx, ny, coord, zcorn, points, cell_merge_thresh, band_thresh, hit_idx
		);
	}

	/*-----------------------------------------------------------------
	* convert points given in (i, j) format to absolute fp coordinates
	*----------------------------------------------------------------*/
	spfp_storarr_t point_index2coord(
		int_t nx, int_t ny, spfp_storarr_t coord, spi_arr_t points_pos, spfp_storarr_t points_param
	) {
		return czt::point_index2coord(
			nx, ny, coord, points_pos, points_param
		);
	}
	/*-----------------------------------------------------------------
	* refine_mesh_deltas for points given in (i,j) format
	*----------------------------------------------------------------*/
	coord_zcorn_pair refine_mesh_deltas(
		int_t& nx, int_t& ny, spfp_storarr_t coord,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx
	) {
		return czt::refine_mesh_deltas(
			nx, ny, coord, points_pos, points_param, cell_merge_thresh, band_thresh, hit_idx
		);
	};

	/*-----------------------------------------------------------------
	* refine_mesh for points given in (i,j) format
	*----------------------------------------------------------------*/
	coord_zcorn_pair refine_mesh(
		int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx
	) {
		return czt::refine_mesh(
			nx, ny, coord, zcorn, points_pos, points_param, cell_merge_thresh, band_thresh, hit_idx
		);
	}

	/*-----------------------------------------------------------------
	* points-based deltas generation algorithm
	*----------------------------------------------------------------*/
	coord_zcorn_pair wave_mesh_deltas_s1(
		fp_stor_t max_dx, fp_stor_t max_dy, fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param
	) {
		return czt::wave_mesh_deltas_s1(
			max_dx, max_dy, len_x, len_y, points_pos, points_param
		);
	}

	void wave_mesh_deltas_s2(
		fp_stor_t cell_dx, fp_stor_t cell_dy,
		fp_stor_t min_dx, fp_stor_t min_dy,
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_t max_sz_tol, bool strict_max_sz
	) {
		return czt::wave_mesh_deltas_s2(
			cell_dx, cell_dy, min_dx, min_dy, dx, dy, max_sz_tol, strict_max_sz
		);
	}

	void wave_mesh_deltas_s2(
		fp_stor_t cell_dx, fp_stor_t cell_dy,
		spfp_storarr_t points_param,
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_t max_sz_tol, bool strict_max_sz
	) {
		return czt::wave_mesh_deltas_s2(
			cell_dx, cell_dy, points_param, dx, dy, max_sz_tol, strict_max_sz
		);
	}

	coord_zcorn_pair wave_mesh_deltas(
		fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param,
		spi_arr_t hit_idx
	) {
		return czt::wave_mesh_deltas(
			max_dx, max_dy, len_x, len_y, points_pos, points_param, hit_idx
		);
	}

	/*-----------------------------------------------------------------
	* points-based mesh generation algorithm
	*----------------------------------------------------------------*/
	coord_zcorn_pair wave_mesh(
		int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param,
		int_t nz, spfp_storarr_t dz,
		fp_stor_t x0, fp_stor_t y0, fp_stor_t z0
	) {
		return czt::wave_mesh(
			nx, ny, max_dx, max_dy, len_x, len_y, points_pos, points_param,
			nz, dz, x0, y0, z0
		);
	}

	// for points in (i,j) format
	coord_zcorn_pair wave_mesh(
		spfp_storarr_t coord,
		int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
		spi_arr_t points_pos, spfp_storarr_t points_param,
		int_t nz, spfp_storarr_t dz,
		fp_stor_t x0, fp_stor_t y0, fp_stor_t z0
	) {
		return czt::wave_mesh(
			coord, nx, ny, max_dx, max_dy, points_pos, points_param,
			nz, dz, x0, y0, z0
		);
	}

	/*-----------------------------------------------------------------
	* find hit index based on given dx, dy
	*----------------------------------------------------------------*/
	spi_arr_t find_hit_idx(
		spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t points_pos,
		fp_stor_t x0, fp_stor_t y0
	) {
		return czt::find_hit_idx(
			dx, dy, points_pos, x0, y0
		);
	}

	spi_arr_t find_hit_idx(
		uint_t nx, uint_t ny, spfp_storarr_t coord,
		spfp_storarr_t points_pos
	) {
		return czt::find_hit_idx(nx, ny, coord, points_pos);
	}

	// mesh refine with wave algorithm
	// returns refined coord & zcorn
	coord_zcorn_pair refine_wave_mesh(
		int_t& nx, int_t& ny,
		spfp_storarr_t coord, spfp_storarr_t zcorn,
		fp_stor_t max_dx, fp_stor_t max_dy,
		spfp_storarr_t points_pos, spfp_storarr_t points_param
	) {
		return czt::refine_wave_mesh(
			nx, ny, coord, zcorn, max_dx, max_dy,
			points_pos, points_param
		);
	}

	coord_zcorn_pair refine_wave_mesh_deltas(
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_stor_t max_dx, fp_stor_t max_dy,
		spfp_storarr_t points_pos, spfp_storarr_t points_param
	) {
		return czt::refine_wave_mesh_deltas(
			dx, dy, max_dx, max_dy,
			points_pos, points_param
		);
	}

	// convert array of cell's tops coord to structured grid representation
	spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float tops) {
		return czt::tops2struct_grid(nx, ny, tops);
	}
	// memory-efficient implementation directly from COORD and ZCORN using tops_iterator
	spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float coord, spv_float zcorn) {
		return czt::tops2struct_grid(nx, ny, coord, zcorn);
	}
	// same as above using provided rs_smesh_iface object
	spv_float tops2struct_grid(smart_ptr< rs_smesh_iface > mesh) {
		return czt::tops2struct_grid(mesh);
	}

	spv_float gen_sgrid(
		t_ulong nx, t_ulong ny, t_ulong nz,
		spv_float dx, spv_float dy, spv_float dz,
		t_float x0, t_float y0, t_float z0,
		bool zyx_order
	) {
		return czt::gen_sgrid(
			nx, ny, nz, dx, dy, dz, x0, y0, z0, zyx_order
		);
	}

	spv_float gen_sgrid_2d(
		t_ulong nx, t_ulong ny,
		spv_float dx, spv_float dy,
		t_float x0, t_float y0,
		bool yx_order
	) {
		return czt::gen_sgrid_2d(
			nx, ny, dx, dy, x0, y0, yx_order
		);
	}

	BS_RESOLVE_TYPE_IMPL_MEM
};

} // eof hidden namespace

/*-----------------------------------------------------------------
 * BS-related implementation stuff
 *----------------------------------------------------------------*/
BLUE_SKY_TYPE_IMPL(czt_iface, objbase, "czt_iface", "Interface to CZT functions", "")
BLUE_SKY_TYPE_STD_CREATE_IFACE(czt_iface, czt_object)
BLUE_SKY_TYPE_STD_COPY_IFACE(czt_iface, czt_object)

bool register_czt_iface(const plugin_descriptor& pd) {
	return BS_KERNEL.register_type(pd, czt_iface::bs_type());
}

} /* namespace blue_sky */

