/**
 *       \file  py_reservoir_simulator.cpp
 *      \brief  Exports python wrappers for reservoir_simulator,
 *              see reservoir_simulator.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"

#include "py_reservoir_simulator.h"
#include "reservoir_simulator.h"

#include "calc_model.h"
#include "jacobian.h"
#include "reservoir.h"
#include "facility_manager.h"
#include "data_storage_interface.h"

#include "py_reservoir.h"
#include "py_facility_manager.h"
#include "keyword_manager.h"
#include "well_connection.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "export_python_wrapper.h"
#include <boost/python/call_method.hpp>

using namespace boost::python;

namespace blue_sky {
namespace python {

  BS_SP (bcsr_matrix_iface)
  get_matrix (reservoir_simulator *t, std::string const &name)
  {
    return t->get_jacobian ()->get_matrix (name);
  }

  reservoir::sp_facility_manager_t
  get_facility_list (reservoir_simulator *t)
  {
    return t->get_reservoir ()->get_facility_list ();
  }

  bool
  subscribe (reservoir_simulator *t, int signal_code, const python_slot &slot)
  {
    bool result = t->subscribe (signal_code, slot.spslot);
    result &=     t->subscribe (objbase::on_delete, new blue_sky::tools::py_object_handler (boost::python::detail::wrapper_base_::get_owner (slot)));

    return result;
  }

  PY_EXPORTER (reservoir_simulator_exporter, default_exporter)
    .def ("init",                     &T::read_keyword_file_and_init)
    .def ("simulate",                 &T::main_loop)
    .def ("simulate",                 &T::simulate)
    .def ("subscribe",                subscribe)
    .add_property ("facility_list",   make_function (get_facility_list))
    .def ("matrix",                   make_function (get_matrix))
    .add_property ("calc_model",      &T::get_calc_model)
    .add_property ("reservoir",       &T::get_reservoir)
    .add_property ("event_manager",   &T::get_event_manager)
    .add_property ("hydrodynamic_model",    &T::get_hydrodynamic_model)
    .add_property ("jacobian",        &T::get_jacobian)
    ;
  reservoir_simulator::python_exporter (class__); // hmm, what changed that I should add this line to old model_loader.py works?
  PY_EXPORTER_END;

  void 
  py_export_reservoir_simulator()
  {
      base_exporter <reservoir_simulator, reservoir_simulator_exporter>::export_class ("reservoir_simulator");
  }

} // namespace python
} // namespace blue_sky
#endif
