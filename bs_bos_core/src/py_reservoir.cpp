/**
 *       \file  py_reservoir.cpp
 *      \brief  Export python wrappers for reservoir,
 *              see reservoir.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"

#include "py_reservoir.h"
#include "reservoir.h"
#include "well_connection.h"
#include "py_facility_manager.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

  PY_EXPORTER (reservoir_exporter, default_exporter)
    .def ("add_filter_well",                        &T::add_filter_well)
    .add_property ("facility_list",                 &T::get_facility_list)
    .add_property ("well_factory",                  make_function (&T::get_well_factory),                 make_function (&T::set_well_factory))
    .add_property ("well_controller_factory",       make_function (&T::get_well_controller_factory),      make_function (&T::set_well_controller_factory))
    .add_property ("well_limit_operation_factory",  make_function (&T::get_well_limit_operation_factory), make_function (&T::set_well_limit_operation_factory))
  PY_EXPORTER_END;

  void py_export_reservoir()
  {
    base_exporter <reservoir, reservoir_exporter>::export_class ("reservoir");
  }

} // namespace python
} // namespace blue_sky
#endif
