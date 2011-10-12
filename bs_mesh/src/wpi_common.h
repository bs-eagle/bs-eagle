/// @file wpi_common.h
/// @brief Global entities needed by WPI algos
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_COMMON_Y614D09T
#define WPI_COMMON_Y614D09T

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Object.h>

#include <iterator>
#include <cmath>

#include "conf.h"
//class bs_mesh_grdecl;

namespace blue_sky { namespace wpi {

// common typedefs
typedef t_ulong ulong;
typedef t_uint uint;
typedef v_float::iterator vf_iterator;

typedef t_float cell_pos[8][3];

// CGAL common typedefs
typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

}} /* { namespace blue_sky::wpi */

#endif /* end of include guard: WPI_COMMON_Y614D09T */

