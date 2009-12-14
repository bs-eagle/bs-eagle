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

#include "export_python_wrapper.h"

#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <datetime.h>

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

  PY_EXPORTER (event_manager_exporter, default_exporter)
    .add_property ("event_list", &T::event_list)
  PY_EXPORTER_END;

  template <typename T>
  void
  list_begin (T *t)
  {
  }

  template <typename T>
  void
  export_event_list (const char *name)
  {
    using namespace boost::python;

    class_ <T> (name)
      .def ("begin", list_begin <T>)
      ;
  }

  template <typename T>
  void
  export_event_map (const char *name)
  {
    using namespace boost::python;

    class_ <T> (name)
      .def (map_indexing_suite <T> ())
      ;
  }

  struct boost_ptime_to_python_datetime
  {
    static PyObject *
    convert (const boost::posix_time::ptime &pt)
    {
      boost::gregorian::date date = pt.date ();
      boost::posix_time::time_duration td = pt.time_of_day ();

      int year    = date.year ();
      int month   = date.month ();
      int day     = date.day ();
      int hour    = td.hours ();
      int minute  = td.minutes ();
      int second  = td.seconds ();

      return PyDateTime_FromDateAndTime (year, month, day, hour, minute, second, 0);
    }
  };

  struct boost_ptime_from_python_datetime
  {
    boost_ptime_from_python_datetime ()
    {
      boost::python::converter::registry::push_back (
        &convertible, &construct,
        boost::python::type_id <boost::posix_time::ptime> ());
    }

    static void *
    convertible (PyObject *obj_ptr)
    {
      return PyDateTime_Check (obj_ptr) ? obj_ptr : 0;
    }

    static void
    construct (PyObject *obj_ptr, 
      boost::python::converter::rvalue_from_python_stage1_data *data)
    {
      const PyDateTime_DateTime *pydate = reinterpret_cast <PyDateTime_DateTime *> (obj_ptr);

      int year      = PyDateTime_GET_YEAR (pydate);
      int month     = PyDateTime_GET_MONTH (pydate);
      int day       = PyDateTime_GET_DAY (pydate);
      int hour      = PyDateTime_DATE_GET_HOUR (pydate);
      int minute    = PyDateTime_DATE_GET_MINUTE (pydate);
      int second    = PyDateTime_DATE_GET_SECOND (pydate);
      int mcsecond  = PyDateTime_DATE_GET_MICROSECOND (pydate);

      boost::gregorian::date d (year, month, day);
      boost::posix_time::time_duration td (hour, minute, second, 0);
      td += boost::posix_time::microseconds (mcsecond);

      using namespace boost::python;
      using namespace boost::posix_time;

      void *storage = ((converter::rvalue_from_python_storage <ptime>*)data)->storage.bytes;

      new (storage) boost::posix_time::ptime (d, td);

      data->convertible = storage;
    }
  };


  void
  py_export_event_manager ()
  {
    using namespace boost::python;

    PyDateTime_IMPORT;

    boost_ptime_from_python_datetime ();
    to_python_converter <boost::posix_time::ptime, boost_ptime_to_python_datetime> ();

    typedef event_base <base_strategy_di> event_base_di_t;
    typedef std::list <smart_ptr <event_base_di_t> > event_list_di_t;

    export_event_list <event_list_di_t> ("event_list_di");
    export_event_map <std::map <boost::posix_time::ptime, event_list_di_t> > ("event_map_di");

    strategy_exporter::export_base <event_manager, event_manager_exporter> ("event_manager");
  }

  template class py_event_manager <base_strategy_fi>;
  template class py_event_manager <base_strategy_di>;
  template class py_event_manager <base_strategy_mixi>;
} // namespace python
} // namespace blue_sky

