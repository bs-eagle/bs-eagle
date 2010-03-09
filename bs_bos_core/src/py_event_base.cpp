/**
 *       \file  py_event_base.cpp
 *      \brief  Python wrappers for event_base, python events iterator,
 *              for event_base see event_base.h
 *     \author  Nikonov Max
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */

#include "stdafx.h"

#include "calc_model.h"
#include "reservoir.h"

#include "py_event_base.h"

//#include "py_reservoir.h"
//#include "py_calc_model.h"

//#include "facility_manager.h"
//#include "calc_model_types.h"
//#include "well_connection.h"

//#include BS_FORCE_PLUGIN_IMPORT ()
//#include "py_rs_mesh.h"
//#include BS_STOP_PLUGIN_IMPORT ()

#ifdef BSPY_EXPORTING_PLUGIN

using namespace boost::python;

namespace blue_sky  {
namespace python    {

  PY_EXPORTER (event_base_exporter, default_exporter)
    .def ("apply", &T::apply)
  PY_EXPORTER_END;

  void py_export_events ()
  {
    strategy_exporter::export_base <event_base, event_base_exporter> ("event_base");
    strategy_exporter::export_class <py_event_base, event_base, event_base_exporter> ("py_event_base");
  }

} // namespace python
} // namespace blue_sky
#endif