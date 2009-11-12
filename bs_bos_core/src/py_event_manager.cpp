/**
 *       \file  py_event_manager.cpp
 *      \brief  Python wrapper for event_manager, 
 *              see event_manager.h
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "event_manager.h"
#include "py_event_manager.h"

#include "reservoir.h"
#include "facility_manager.h"
#include "well_connection.h"
#include "calc_model.h"

using namespace boost;
using namespace boost::posix_time;

namespace blue_sky {
namespace python {

  template <typename strategy_t>
  py_event_manager<strategy_t>::~py_event_manager ()
  {

  }

  template <typename strategy_t>
  py_event_manager<strategy_t>::py_event_manager()
      : py_objbase(wrapped_t::bs_type())
  {
  }

  template <typename strategy_t>
  py_event_manager<strategy_t>::py_event_manager (const sp_em_t &src)
  : py_objbase (src)
  {
  }

  template <typename strategy_t>
  py_event_manager<strategy_t>::py_event_manager(const this_t &src)
  : py_objbase (src.get_spx <wrapped_t> ())
  {}

  template <typename strategy_t>
  typename py_event_manager<strategy_t>::py_event_base_t *
  py_event_manager<strategy_t>::create_event (const std::string &date, const std::string &event_name, const std::string &event_params)
  {
    return new py_event_base_t (this->template get_spx <wrapped_t> ()->create_event (boost::posix_time::ptime (time_from_string (date)), event_name, event_params));
  }

  template <typename strategy_t>
  void py_event_manager<strategy_t>::factory_register(const std::string &/*name*/, const event_creater_base_t &/*creater*/)
  {
    BS_ASSERT (false && "NOT_IMPL_YET");
  }

  template <typename strategy_t>
  event_creater_base<strategy_t>::event_creater_base (boost::python::object &pobj)
      : obj (pobj)
  {}

  template <typename strategy_t>
  typename py_event_manager<strategy_t>::event_list_iterator_t
  py_event_manager<strategy_t>::el_begin ()
  {
    return event_list_iterator_t (this->get_spx <wrapped_t> ());
  }

  //event_list_iterator_t el_end ();

  template <typename strategy_t>
  typename event_creater_base<strategy_t>::py_event_base_t *
  event_creater_base<strategy_t>::create ()
  {
    // CALL PYTHON IMPL (pure virtual)
  }


  template <typename strategy_t>
  void py_export_event_manager_ (const char *name)
  {
    using namespace boost::python;
    class_<py_event_manager <strategy_t>, bases<py_objbase> >(name)
      .def ("create_event", &py_event_manager <strategy_t>::create_event, return_value_policy <manage_new_object> ())
      .def ("el_begin", &py_event_manager <strategy_t>::el_begin)
      ;
  }

  void
  py_export_event_manager ()
  {
    py_export_event_manager_ <base_strategy_fi> ("event_manager_seq_fi");
    py_export_event_manager_ <base_strategy_di> ("event_manager_seq_di");
    py_export_event_manager_ <base_strategy_mixi> ("event_manager_seq_mixi");
  }

  template class py_event_manager <base_strategy_fi>;
  template class py_event_manager <base_strategy_di>;
  template class py_event_manager <base_strategy_mixi>;
} // namespace python
} // namespace blue_sky

