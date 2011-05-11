/// @file py_coord_zcorn_tools.cpp
/// @brief Python bindnings for COORD & ZCORN tools
/// @author uentity
/// @date 10.05.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "coord_zcorn_tools.h"
#include "export_python_wrapper.h"
#include "py_pair_converter.h"
#include <boost/python/tuple.hpp>

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;
namespace czt = blue_sky::coord_zcorn_tools;

namespace {

typedef t_long int_t;
typedef t_double fp_t;
typedef t_float fp_stor_t;
typedef spv_float spfp_storarr_t;
typedef spv_long spi_arr_t;
typedef typename spi_arr_t::pure_pointed_t int_arr_t;
typedef std::pair< spv_float, spv_float > coord_zcorn_pair;

// refine_mesh_deltas with overloads
static tuple refine_mesh_deltas_s(int_t nx, int_t ny, fp_stor_t max_dx, fp_stor_t max_dy,
		fp_stor_t len_x, fp_stor_t len_y,
		spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	std::pair< spfp_storarr_t, spfp_storarr_t > r = czt::refine_mesh_deltas_s(nx, ny, max_dx, max_dy, len_x, len_y, points_pos, points_param);
	return make_tuple(r.first, r.second, nx, ny);
}

}

namespace blue_sky { namespace python {

void py_export_czt() {
	def("refine_mesh_deltas_s", &refine_mesh_deltas_s);
}

}} 	// eof blue_sky::python

#endif

