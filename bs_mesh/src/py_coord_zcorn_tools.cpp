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

namespace {

typedef t_long int_t;
typedef czt::uint_t uint_t;
typedef t_double fp_t;
typedef t_float fp_stor_t;
typedef v_float fp_storarr_t;
typedef spv_float spfp_storarr_t;
typedef spv_long spi_arr_t;
typedef spi_arr_t::pure_pointed_t int_arr_t;
typedef std::pair< spv_float, spv_float > coord_zcorn_pair;

tuple wave_mesh_deltas_s2_o1(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	fp_stor_t min_dx, fp_stor_t min_dy,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol = 0.3, bool strict_max_sz = false)
{
	// Make copies of delta arrays as thay can be referenced by Python
	spfp_storarr_t dx_cpy = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	dx_cpy->assign(*dx);
	spfp_storarr_t dy_cpy = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	dy_cpy->assign(*dy);
	// invoke stage 2
	czt::wave_mesh_deltas_s2(
		cell_dx, cell_dy, min_dx, min_dy,
		dx_cpy, dy_cpy, max_sz_tol, strict_max_sz
	);
	return make_tuple(dx_cpy, dy_cpy);
}

tuple wave_mesh_deltas_s2_o2(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	spfp_storarr_t points_param,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol = 0.3, bool strict_max_sz = false)
{
	// Make copies of delta arrays as thay can be referenced by Python
	spfp_storarr_t dx_cpy = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	dx_cpy->assign(*dx);
	spfp_storarr_t dy_cpy = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	dy_cpy->assign(*dy);
	// invoke stage 2
	czt::wave_mesh_deltas_s2(
		cell_dx, cell_dy, points_param,
		dx, dy, max_sz_tol, strict_max_sz
	);
	return make_tuple(dx_cpy, dy_cpy);
}

// wave_mesh_deltas with overloads
tuple wave_mesh_deltas(fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
	coord_zcorn_pair r = czt::wave_mesh_deltas(
		max_dx, max_dy, len_x, len_y, points_pos, points_param,
		hit_idx
	);
	return make_tuple(r.first, r.second, hit_idx);
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

spi_arr_t find_hit_idx1(
	spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t points_pos,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0)
{
	return czt::find_hit_idx(dx, dy, points_pos, x0, y0);
}

spi_arr_t find_hit_idx2(
	uint_t nx, uint_t ny, spfp_storarr_t coord,
	spfp_storarr_t points_pos)
{
	return czt::find_hit_idx(nx, ny, coord, points_pos);
}

}  // eof hidden namespace

BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_deltas_s2_overl1, wave_mesh_deltas_s2_o1, 6, 8)
BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_deltas_s2_overl2, wave_mesh_deltas_s2_o2, 5, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_overl, wave_mesh, 8, 11)
BOOST_PYTHON_FUNCTION_OVERLOADS(wave_mesh_ij_overl, wave_mesh_ij, 9, 12)
BOOST_PYTHON_FUNCTION_OVERLOADS(find_hit_idx_overl, find_hit_idx1, 3, 5)

namespace blue_sky { namespace python {

void py_export_czt() {
	def("wave_mesh_deltas_s1", &czt::wave_mesh_deltas_s1);
	def("wave_mesh_deltas_s2", &wave_mesh_deltas_s2_o1, wave_mesh_deltas_s2_overl1());
	def("wave_mesh_deltas_s2", &wave_mesh_deltas_s2_o2, wave_mesh_deltas_s2_overl2());
	def("wave_mesh_deltas", &wave_mesh_deltas);
	def("wave_mesh", &wave_mesh, wave_mesh_overl());
	def("wave_mesh", &wave_mesh_ij, wave_mesh_ij_overl());
	def("find_hit_idx", &find_hit_idx1, find_hit_idx_overl());
	def("find_hit_idx", &find_hit_idx2);
	def("refine_wave_mesh", &czt::refine_wave_mesh);
	def("refine_wave_mesh_deltas", &czt::refine_wave_mesh_deltas);
}

}} 	// eof blue_sky::python

#endif

