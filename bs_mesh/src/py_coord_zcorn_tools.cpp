/// @file py_coord_zcorn_tools.cpp
/// @brief Python bindnings for COORD & ZCORN tools
/// @author uentity
/// @date 10.05.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
#include "coord_zcorn_tools.h"
#include "export_python_wrapper.h"
#include "py_pair_converter.h"
#include <boost/python/tuple.hpp>

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;
namespace czt = blue_sky::coord_zcorn_tools;

BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_overl, wave_mesh, 8, 11)
BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_ij_overl, wave_mesh_ij, 9, 12)

namespace {

typedef t_long int_t;
typedef t_double fp_t;
typedef t_float fp_stor_t;
typedef spv_float spfp_storarr_t;
typedef spv_long spi_arr_t;
typedef spi_arr_t::pure_pointed_t int_arr_t;
typedef std::pair< spv_float, spv_float > coord_zcorn_pair;

// wave_mesh_deltas with overloads
tuple wave_mesh_deltas(fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	int_t nx, ny;
	coord_zcorn_pair r = czt::wave_mesh_deltas(
		nx, ny, max_dx, max_dy, len_x, len_y, points_pos, points_param
	);
	return make_tuple(r.first, r.second, nx, ny);
}

tuple wave_mesh(
	fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0)
{
	int_t nx, ny;
	coord_zcorn_pair r = czt::wave_mesh(
		nx, ny, max_dx, max_dy, len_x, len_y, points_pos, points_param,
		nz, dz, x0, y0, z0
	);
	return make_tuple(r.first, r.second, nx, ny);
}

tuple wave_mesh_ij(
	spfp_storarr_t coord,
	int_t nx, int_t ny,
	fp_stor_t max_dx, fp_stor_t max_dy,
	spi_arr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0, fp_stor_t z0 = 0)
{
	coord_zcorn_pair r = czt::wave_mesh(
		coord,
		nx, ny, max_dx, max_dy, points_pos, points_param,
		nz, dz, x0, y0, z0
	);
	return make_tuple(r.first, r.second, nx, ny);
}

}

namespace blue_sky { namespace python {

void py_export_czt() {
	def("wave_mesh_deltas", &wave_mesh_deltas);
	def("wave_mesh", &wave_mesh, wave_mesh_overl());
	def("wave_meshj", &wave_mesh_ij, wave_mesh_ij_overl());
}

}} 	// eof blue_sky::python

#endif

