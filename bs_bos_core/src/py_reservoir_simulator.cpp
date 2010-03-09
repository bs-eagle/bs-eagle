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

  template <typename T>
  typename T::sp_jmatrix_t 
  get_jmatrix (T *t)
  {
    return t->jacobian_->get_jmatrix ();
  }

  template <typename T>
  typename T::reservoir_t::sp_facility_manager_t
  get_facility_list (T *t)
  {
    return t->reservoir_->get_facility_list ();
  }

  template <typename T>
  bool
  subscribe (T *t, int signal_code, const python_slot &slot)
  {
    bool result = t->subscribe (signal_code, slot.spslot);
    result &=     t->subscribe (objbase::on_delete, new blue_sky::tools::py_object_handler (boost::python::detail::wrapper_base_::get_owner (slot)));

    return result;
  }

  PY_EXPORTER (reservoir_simulator_exporter, default_exporter)
    .def ("init",                     &T::read_keyword_file_and_init)
    .def ("simulate",                 &T::main_loop)
    .def ("simulate",                 &T::simulate)
    .def ("subscribe",                subscribe <T>)
    .add_property ("facility_list",   make_function (get_facility_list <T>))
    .add_property ("jmatrix",         make_function (get_jmatrix <T>))
    .add_property ("calc_model",      make_getter (&T::cm, return_value_policy <copy_non_const_reference> ()))
    .add_property ("reservoir",       make_getter (&T::reservoir_, return_value_policy <copy_non_const_reference> ()))
    .add_property ("mesh",            make_getter (&T::mesh, return_value_policy <copy_non_const_reference> ()))
    .add_property ("event_manager",   make_getter (&T::em, return_value_policy <copy_non_const_reference> ()))
    .add_property ("data_manager",    make_getter (&T::dm, return_value_policy <copy_non_const_reference> ()))
    .add_property ("jacobian",        make_getter (&T::jacobian_, return_value_policy <copy_non_const_reference> ()))
    .add_property ("keyword_manager", make_getter (&T::keyword_manager_, return_value_policy <copy_non_const_reference> ()))
  PY_EXPORTER_END;

  void 
  py_export_reservoir_simulator()
  {
    strategy_exporter::export_base <reservoir_simulator, reservoir_simulator_exporter> ("reservoir_simulator");
  }

} // namespace python
} // namespace blue_sky
#endif