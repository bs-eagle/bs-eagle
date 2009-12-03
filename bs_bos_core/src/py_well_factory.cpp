/**
 *       \file  py_well_factory.cpp
 *      \brief  Python wrappers for well factories (well, well_controller,
 *              well_rate_controller, well_limit_operation factories),
 *              see calc_well.h and related files
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "py_well_factory.h"
#include "well_connection.h"

namespace blue_sky {
namespace python {

  PY_EXPORTER (well_factory_exporter, default_exporter)
    .def ("create_well", &T::create_well)
    .def ("create_connection", &T::create_connection)
  PY_EXPORTER_END;

  PY_EXPORTER (well_controller_factory_exporter, default_exporter)
    .def ("create_controller", &T::create_controller)
  PY_EXPORTER_END;

  void 
  py_export_well_factories ()
  {
    strategy_exporter::export_base <well_factory, well_factory_exporter> ("well_factory");
    strategy_exporter::export_base <well_controller_factory, well_controller_factory_exporter> ("well_controller_factory");

    strategy_exporter::export_class <py_well_factory, well_factory, well_factory_exporter> ("py_well_factory");
    strategy_exporter::export_class <py_well_controller_factory, well_controller_factory, well_controller_factory_exporter> ("py_well_controller_factory");
  }

} // namespace python
} // namespace blue_sky

