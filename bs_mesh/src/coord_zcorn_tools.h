/// @file coord_zcorn_tools.h
/// @brief COORD & ZCORN generation and refinement tools
/// @author uentity
/// @date 10.05.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef COORD_ZCORN_TOOLS_1M4YIRKV
#define COORD_ZCORN_TOOLS_1M4YIRKV

#include "bs_mesh_stdafx.h"
#include "rs_smesh_iface.h"
#include <iterator>
#include <cmath>

/*-----------------------------------------------------------------
 * helpers for gen_coord_zcorn & refine_mesh
 *----------------------------------------------------------------*/
namespace blue_sky { namespace coord_zcorn_tools {

// shorter aliases
typedef t_long int_t;
typedef t_ulong uint_t;
typedef t_double fp_t;
typedef t_float fp_stor_t;

// other typedefs
typedef v_float fp_storarr_t;
typedef spv_float spfp_storarr_t;
typedef spv_long spi_arr_t;
typedef std::pair< spfp_storarr_t, spfp_storarr_t > coord_zcorn_pair;

/*-----------------------------------------------------------------
 * Main processing functions
 *----------------------------------------------------------------*/
BS_API_PLUGIN spfp_storarr_t gen_coord(int_t nx, int_t ny, spfp_storarr_t dx, spfp_storarr_t dy, fp_t x0, fp_t y0);

BS_API_PLUGIN coord_zcorn_pair refine_mesh_deltas(int_t& nx, int_t& ny, spfp_storarr_t coord,
	spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
	spi_arr_t hit_idx = NULL);

BS_API_PLUGIN coord_zcorn_pair refine_mesh(int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx = NULL);

/*-----------------------------------------------------------------
 * convert points given in (i, j) format to absolute fp coordinates
 *----------------------------------------------------------------*/
spfp_storarr_t point_index2coord(int_t nx, int_t ny, spfp_storarr_t coord, spi_arr_t points_pos, spfp_storarr_t points_param);
/*-----------------------------------------------------------------
 * refine_mesh_deltas for points given in (i,j) format
 *----------------------------------------------------------------*/
BS_API_PLUGIN coord_zcorn_pair refine_mesh_deltas(int_t& nx, int_t& ny, spfp_storarr_t coord,
	spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
	spi_arr_t hit_idx = NULL);

/*-----------------------------------------------------------------
 * refine_mesh for points given in (i,j) format
 *----------------------------------------------------------------*/
BS_API_PLUGIN coord_zcorn_pair refine_mesh(int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx);

/*-----------------------------------------------------------------
 * points-based deltas generation algorithm
 *----------------------------------------------------------------*/
BS_API_PLUGIN coord_zcorn_pair wave_mesh_deltas_s1(
	fp_stor_t max_dx, fp_stor_t max_dy, fp_stor_t len_x, fp_stor_t len_y,
	spfp_storarr_t points_pos, spfp_storarr_t points_param);

BS_API_PLUGIN void wave_mesh_deltas_s2(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	fp_stor_t min_dx, fp_stor_t min_dy,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol = 0.3, bool strict_max_sz = false);

BS_API_PLUGIN void wave_mesh_deltas_s2(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	spfp_storarr_t points_param,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol = 0.3, bool strict_max_sz = false);

BS_API_PLUGIN coord_zcorn_pair wave_mesh_deltas(
	fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y,
	spfp_storarr_t points_pos, spfp_storarr_t points_param,
	spi_arr_t hit_idx = NULL);

/*-----------------------------------------------------------------
 * points-based mesh generation algorithm
 *----------------------------------------------------------------*/
BS_API_PLUGIN coord_zcorn_pair wave_mesh(
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0);

// for points in (i,j) format
BS_API_PLUGIN coord_zcorn_pair wave_mesh(
	spfp_storarr_t coord,
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	spi_arr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0);

/*-----------------------------------------------------------------
 * find hit index based on given dx, dy
 *----------------------------------------------------------------*/
BS_API_PLUGIN spi_arr_t find_hit_idx(
	spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t points_pos,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0);

BS_API_PLUGIN spi_arr_t find_hit_idx(
	uint_t nx, uint_t ny, spfp_storarr_t coord,
	spfp_storarr_t points_pos);

// mesh refine with wave algorithm
// returns refined coord & zcorn
BS_API_PLUGIN coord_zcorn_pair refine_wave_mesh(
	int_t& nx, int_t& ny,
	spfp_storarr_t coord, spfp_storarr_t zcorn,
	fp_stor_t max_dx, fp_stor_t max_dy,
	spfp_storarr_t points_pos, spfp_storarr_t points_param);

BS_API_PLUGIN coord_zcorn_pair refine_wave_mesh_deltas(
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_stor_t max_dx, fp_stor_t max_dy,
	spfp_storarr_t points_pos, spfp_storarr_t points_param);

// convert array of cell's tops coord to structured grid representation
spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float tops);
// memory-efficient implementation directly from COORD and ZCORN using tops_iterator
spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float coord, spv_float zcorn);
// same as above using provided rs_smesh_iface object
spv_float tops2struct_grid(smart_ptr< rs_smesh_iface > mesh);

spv_float gen_sgrid(
	t_ulong nx, t_ulong ny, t_ulong nz,
	spv_float dx, spv_float dy, spv_float dz,
	t_float x0 = 0, t_float y0 = 0, t_float z0 = 0,
	bool zyx_order = false
);

spv_float gen_sgrid_2d(
	t_ulong nx, t_ulong ny,
	spv_float dx, spv_float dy,
	t_float x0 = 0, t_float y0 = 0,
	bool yx_order = false
);

}}  // eof namespace blue_sky::coord_zcorn_tools

#endif /* end of include guard: COORD_ZCORN_TOOLS_1M4YIRKV */

