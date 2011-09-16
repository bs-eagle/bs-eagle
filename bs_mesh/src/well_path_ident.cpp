/// @file well_path_ident.cpp
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 05.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "well_path_ident.h"
#include "wpi_strategy_3d.h"
#include "wpi_strategy_2d.h"
#include "wpi_algo.h"

#include "bs_mesh_stdafx.h"
#include "export_python_wrapper.h"

namespace blue_sky {
// alias
namespace bp = boost::python;

// specialization for 3D
spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return wpi::wpi_algo< wpi::wpi_strategy_3d >::well_path_ident_d(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
}

// specialization for 2D
spv_float well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	//return well_path_ident_(nx, ny, coord, zcorn, well_info, include_well_nodes);
	return wpi::wpi_algo< wpi::wpi_strategy_2d >::well_path_ident_d(
		nx, ny, coord, zcorn, well_info, include_well_nodes
	);
}

/*-----------------------------------------------------------------
 * Python bindings
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl, well_path_ident, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl_2d, well_path_ident_2d, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl_2d_old, well_path_ident_2d_old, 5, 6)

namespace python {

void py_export_wpi() {
	bp::def("well_path_ident", &well_path_ident, well_path_ident_overl());
	bp::def("well_path_ident_2d", &well_path_ident_2d, well_path_ident_overl_2d());
	bp::def("well_path_ident_2d_old", &well_path_ident_2d_old, well_path_ident_overl_2d_old());
}

}

}	// eof blue-sky namespace

