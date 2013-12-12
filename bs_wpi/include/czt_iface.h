/// @file czt_iface.h
/// @brief BS interface for CZT functions
/// @author uentity
/// @version 1.0
/// @date 12.12.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef CZT_IFACE_1510IJ8B
#define CZT_IFACE_1510IJ8B

#include "bs_kernel.h"
#include "coord_zcorn_tools.h"

namespace blue_sky {
namespace czt = coord_zcorn_tools;

class BS_API_PLUGIN czt_iface : public objbase {
public:
	typedef czt::int_t int_t;
	typedef czt::uint_t uint_t;
	typedef czt::fp_t fp_t;
	typedef czt::fp_stor_t fp_stor_t;

	typedef czt::fp_storarr_t fp_storarr_t;
	typedef czt::spfp_storarr_t spfp_storarr_t;
	typedef czt::spi_arr_t spi_arr_t;
	typedef czt::coord_zcorn_pair coord_zcorn_pair;

	/*-----------------------------------------------------------------
	* Main processing functions
	*----------------------------------------------------------------*/
	virtual spfp_storarr_t gen_coord(
		int_t nx, int_t ny, spfp_storarr_t dx, spfp_storarr_t dy, fp_t x0, fp_t y0
	) = 0;

	virtual coord_zcorn_pair gen_coord_zcorn(
		int_t nx, int_t ny, int_t nz,
		spv_float dx, spv_float dy, spv_float dz,
		fp_stor_t x0, fp_stor_t y0, fp_stor_t z0
	) = 0;

	virtual spfp_storarr_t gen_coord2(spfp_storarr_t x, spfp_storarr_t y) = 0;

	virtual coord_zcorn_pair refine_mesh_deltas(
		int_t& nx, int_t& ny, spfp_storarr_t coord,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx = NULL
	) = 0;

	virtual coord_zcorn_pair refine_mesh(
		int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx = NULL
	) = 0;

	/*-----------------------------------------------------------------
	* convert points given in (i, j) format to absolute fp coordinates
	*----------------------------------------------------------------*/
	virtual spfp_storarr_t point_index2coord(
		int_t nx, int_t ny, spfp_storarr_t coord, spi_arr_t points_pos, spfp_storarr_t points_param
	) = 0;
	/*-----------------------------------------------------------------
	* refine_mesh_deltas for points given in (i,j) format
	*----------------------------------------------------------------*/
	virtual coord_zcorn_pair refine_mesh_deltas(
		int_t& nx, int_t& ny, spfp_storarr_t coord,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx = NULL
	) = 0;

	/*-----------------------------------------------------------------
	* refine_mesh for points given in (i,j) format
	*----------------------------------------------------------------*/
	virtual coord_zcorn_pair refine_mesh(
		int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx
	) = 0;

	/*-----------------------------------------------------------------
	* points-based deltas generation algorithm
	*----------------------------------------------------------------*/
	virtual coord_zcorn_pair wave_mesh_deltas_s1(
		fp_stor_t max_dx, fp_stor_t max_dy, fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param) = 0;

	virtual void wave_mesh_deltas_s2(
		fp_stor_t cell_dx, fp_stor_t cell_dy,
		fp_stor_t min_dx, fp_stor_t min_dy,
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_t max_sz_tol = 0.3, bool strict_max_sz = false) = 0;

	virtual void wave_mesh_deltas_s2(
		fp_stor_t cell_dx, fp_stor_t cell_dy,
		spfp_storarr_t points_param,
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_t max_sz_tol = 0.3, bool strict_max_sz = false) = 0;

	virtual coord_zcorn_pair wave_mesh_deltas(
		fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param,
		spi_arr_t hit_idx = NULL) = 0;

	/*-----------------------------------------------------------------
	* points-based mesh generation algorithm
	*----------------------------------------------------------------*/
	virtual coord_zcorn_pair wave_mesh(
		int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param,
		int_t nz, spfp_storarr_t dz,
		fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0
	) = 0;

	// for points in (i,j) format
	virtual coord_zcorn_pair wave_mesh(
		spfp_storarr_t coord,
		int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
		spi_arr_t points_pos, spfp_storarr_t points_param,
		int_t nz, spfp_storarr_t dz,
		fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0
	) = 0;

	/*-----------------------------------------------------------------
	* find hit index based on given dx, dy
	*----------------------------------------------------------------*/
	virtual spi_arr_t find_hit_idx(
		spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t points_pos,
		fp_stor_t x0 = 0, fp_stor_t y0 = 0
	) = 0;

	virtual spi_arr_t find_hit_idx(
		uint_t nx, uint_t ny, spfp_storarr_t coord,
		spfp_storarr_t points_pos
	) = 0;

	// mesh refine with wave algorithm
	// returns refined coord & zcorn
	virtual coord_zcorn_pair refine_wave_mesh(
		int_t& nx, int_t& ny,
		spfp_storarr_t coord, spfp_storarr_t zcorn,
		fp_stor_t max_dx, fp_stor_t max_dy,
		spfp_storarr_t points_pos, spfp_storarr_t points_param
	) = 0;

	virtual coord_zcorn_pair refine_wave_mesh_deltas(
		spfp_storarr_t dx, spfp_storarr_t dy,
		fp_stor_t max_dx, fp_stor_t max_dy,
		spfp_storarr_t points_pos, spfp_storarr_t points_param
	) = 0;

	// convert array of cell's tops coord to structured grid representation
	virtual spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float tops) = 0;
	// memory-efficient implementation directly from COORD and ZCORN using tops_iterator
	virtual spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float coord, spv_float zcorn) = 0;
	// same as above using provided rs_smesh_iface object
	virtual spv_float tops2struct_grid(smart_ptr< rs_smesh_iface > mesh) = 0;

	virtual spv_float gen_sgrid(
		t_ulong nx, t_ulong ny, t_ulong nz,
		spv_float dx, spv_float dy, spv_float dz,
		t_float x0 = 0, t_float y0 = 0, t_float z0 = 0,
		bool zyx_order = false
	) = 0;

	virtual spv_float gen_sgrid_2d(
		t_ulong nx, t_ulong ny,
		spv_float dx, spv_float dy,
		t_float x0 = 0, t_float y0 = 0,
		bool yx_order = false
	) = 0;

	BLUE_SKY_TYPE_DECL(czt_iface)
};

} /* namespace blue_sky */

#endif /* end of include guard: CZT_IFACE_1510IJ8B */

