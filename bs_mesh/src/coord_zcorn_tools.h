/// @file coord_zcorn_tools.h
/// @brief COORD & ZCORN generation and refinement tools
/// @author uentity
/// @date 10.05.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
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


BS_API_PLUGIN coord_zcorn_pair wave_mesh_deltas(
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param);

}}  // eof namespace blue_sky::coord_zcorn_tools

