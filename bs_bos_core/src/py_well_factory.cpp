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
    .def ("create_control", &T::create_control)
    .def ("set_rate_control_factory", &T::set_rate_control_factory)
  PY_EXPORTER_END;

  PY_EXPORTER (well_rate_control_factory_exporter, default_exporter)
    .def ("create_control", &T::create_control)
  PY_EXPORTER_END;

  //PYTHON_EXPORT_WRAPPER (py_export_well_factory, well_factory, py_well_factory, 2, (create_well, create_connection));
  //PYTHON_EXPORT_WRAPPER (py_export_well_controller_factory, well_controller_factory, py_well_controller_factory, 3, (create_controller, create_control, set_rate_control_factory));
  ////PYTHON_EXPORT_WRAPPER (py_export_well_limit_operation_factory, well_limit_operation_factory, py_well_limit_operation_factory, 1, (create_limit));

  //PYTHON_EXPORT_WRAPPER (py_export_well_rate_control_factory, well_rate_control_factory, py_well_rate_control_factory, 1, (create_control));

  void 
  py_export_well_factories ()
  {
    strategy_exporter::export_base <well_factory, well_factory_exporter> ("well_factory");
    strategy_exporter::export_base <well_controller_factory, well_controller_factory_exporter> ("well_controller_factory");
    strategy_exporter::export_base <well_rate_control_factory, well_rate_control_factory_exporter> ("well_rate_control_factory");

    strategy_exporter::export_class <py_well_factory, well_factory, well_factory_exporter> ("py_well_factory");
    strategy_exporter::export_class <py_well_controller_factory, well_controller_factory, well_controller_factory_exporter> ("py_well_controller_factory");
    strategy_exporter::export_class <py_well_rate_control_factory, well_rate_control_factory, well_rate_control_factory_exporter> ("py_well_rate_control_factory");
  }

} // namespace python
} // namespace blue_sky

