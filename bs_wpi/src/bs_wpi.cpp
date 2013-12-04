/// @file bs_wpi.cpp
/// @brief Entry point for bs_wpi shared library
/// @author uentity
/// @version 1.0
/// @date 28.11.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_kernel.h"
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python/module.hpp>
#endif

BLUE_SKY_PLUGIN_DESCRIPTOR_EXT("bs_wpi", "1.0.0", "Well path identification", "", "wpi")

namespace blue_sky {
namespace {

bool register_types(const plugin_descriptor& pd) {
	// no types for a while
	return true;
}

} /* eof hidden namespace */

BLUE_SKY_REGISTER_PLUGIN_FUN {
	return register_types(*bs_init.pd_);
}

#ifdef BSPY_EXPORTING_PLUGIN
/*-----------------------------------------------------------------
 * forward declaration of exporting functions
 *----------------------------------------------------------------*/
namespace python {

void py_export_wpi();
void py_export_wpi_vtk();

}

/*-----------------------------------------------------------------
 * callback that make Python bindings
 *----------------------------------------------------------------*/
namespace {

void init_py_subsystem() {
	python::py_export_wpi();
	python::py_export_wpi_vtk();
}

} /* eof hidden namespace */

BLUE_SKY_INIT_PY_FUN {
	init_py_subsystem();
}

#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_wpi_d)
#else
BOOST_PYTHON_MODULE (bs_wpi)
#endif
{
	init_py_subsystem ();
	std::cout << &BS_KERNEL << std::endl;
	bool res = blue_sky::register_types (*bs_get_plugin_descriptor());
	if (!res)
		throw "Can't register bs_wpi types";
}

#endif
} /* namespace blue_sky */

